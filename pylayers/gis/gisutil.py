# -*- coding: utf-8 -*-
import numpy as np
import pdb
import simplejson
import urllib
import requests
import json
"""

.. currentmodule:: pylayers.gis.gisutil

.. autosummary::

"""
def enctile(lon,lat):
    """ encode tile prefix from (lon,lat)

    Parameters
    ----------

    lon : float
        longitude (degrees)
    lat : float
        latitude (degrees)

    Returns
    -------

    prefix : string
        srtm prefix filename

    Examples
    --------

        >>> from pylayers.gis.gisutil import *
        >>> assert enctile(-1.5,48.5)=='N48W002'
        >>> assert enctile(0.5,48.5)=='N48E000'

    """
    if lon>0:
        slon='E'
    else:
        slon='W'

    if lat>0:
        slat='N'
    else:
        slat='S'
    # rounding and integer to string  conversion
    if lat>0:
        clat = str(np.floor(abs(lat)).astype('int'))
    else:
        clat = str(np.ceil(abs(lat)).astype('int'))
    if lon>0:
        clon = str(np.floor(abs(lon)).astype('int'))
    else:
        clon = str(np.ceil(abs(lon)).astype('int'))
    # 2 char formating
    if len(clon)<2:
        clon ='00'+clon
    else:
        clon='0'+clon
    if len(clat)<2:
        clat ='0'+clat

    tilename = slat+clat+slon+clon
    return tilename

def dectile(prefix='N48W002'):
    """ decode tile name

    Parameters
    ----------

    prefix : string

    Returns
    -------

    extent : tuple

        (lonmin,lonmax,latmin,latmax)

    Examples
    --------

    >>> from pylayers.gis.gisutil import *
    >>> dectile('N48W002')
    (-2.0, -1.0, 48,49)

    """
    if prefix[0] == 'N':
        latmin = int(prefix[1:3])
        latmax = latmin+1
    if prefix[0] == 'S':
        latmin = -int(prefix[1:3])
        latmax = latmin+1

    if prefix[3] == 'W':
        st = prefix[4:]
        if st[0] == '0':
            lonmin = -int(st[1:])
        else:
            lonmin = -int(st)
        lonmax = lonmin+1

    if prefix[3] == 'E':
        st = prefix[4:]
        if st[0] == '0':
            lonmin = int(st[1:])
        else:
            lonmin = int(st)

        lonmax = lonmin+1
    return (lonmin, lonmax, latmin, latmax)


def ent(lL):
    """ encode lon Lat in natural integer

    Parameters
    ----------

    lL : nd.array (N,2)
        longitude latitude array


    Returns
    -------

    'iA-B' : array of string (2xN)

    Longitude and latitude offset (a,b) with a<1 and b<1
    are multiplied by 256 and converted in unsigned int8.
    This means that the 1 degree  tile is subdivided into 256*256 subtiles.
    'i0-0' is the lower corner subtile
    'i255-255' is the upper corner subtile

    Examples
    --------

    >>> from pylayers.gis.gisutil import *
    >>> lL= np.array([[-2,48],[-1,49]])
    >>> ent(lL)

    """
    # offset from the lower left corner
    if len(lL.shape) != 2:
        lL = lL[None, :]
    d = lL-np.floor(lL)
    dui8 = np.floor(d*255).astype('uint8')
    lab = ['i'+str(x[0])+'-'+str(x[1]) for x in dui8]
    return(lab)


def eqt(lL):
    """ encode lon Lat in quad tree integer

    Parameters
    ----------

    lL : nd.array (2xN)
        longitude Latitude

    Returns
    -------
    ud16 : array of uint16

    """
    N = lL.shape[1]
    # offset from the lower left corner
    d = lL-np.floor(lL)
    dui8 = np.floor(d*256).astype('uint8')
    ndui8 = np.unpackbits(dui8).reshape(2, N, 8)
    d16 = np.empty((N, 16)).astype('int')
    d16[:, 1::2] = ndui8[0, :, :]
    d16[:, 0::2] = ndui8[1, :, :]
    ud8 = np.packbits(d16, axis=1)
    ud16 = (ud8[:, 0]*256+ud8[:, 1]).astype('uint16')
    return(ud16)


def dqt(ud16, lL0):
    """ decode quad tree integer to lon Lat

    Parameters
    ----------

    ud16 : nd.array (2xN)
        longitude Latitude
    lL0 : nd.array (,2)
        lower left corner of the 1degree tile

    Returns
    -------

    lL : longitude, Latitude

    """

    N = len(ud16)
    # offset from the lower left corner
    #d = lL-lL0[:,None]
    #dui8 = np.floor(d*256).astype('uint8')
    uh8 = ud16/256
    ul8 = ud16-uh8*256
    ud8 = (np.vstack((uh8,ul8)).T).astype('uint8')
    ud16 = np.unpackbits(ud8).reshape(N,16)
    ndu8 = np.empty((2,N,8)).astype('int')
    ndu8[0,:,:]=ud16[:,1::2]
    ndu8[1,:,:]=ud16[:,0::2]
    du8 = np.packbits(ndu8).reshape(2,N)/256.
    lL = lL0[:,None]+du8

    return(lL)


def lab2ext(lab, lL0):
    """ label to extent

    Parameters
    ----------

    lab : string
        'i0-0' to 'i255-255'
    lL0 : tuple
        (lon0,lat0)

    Returns
    -------

    ext : extent
        (lonmin,lonmax,latmin,Latmax)

    """
    lab=lab.replace('i','')
    ilon = int(lab.split('-')[0])
    ilat = int(lab.split('-')[1])
    lonmin = lL0[0] + ilon/256
    lonmax = lL0[0] + (ilon+1)/256
    latmin = lL0[1] + ilat/256
    latmax = lL0[1] + (ilat+1)/256

    ext = (lonmin, lonmax, latmin, latmax)
    return ext

def lL2ext(lL):
    """ convert a lonLat into its qt extent

    """
    label = ent(lL)
    ext = [lab2ext(x , np.floor(lL)) for x in label]
    return ext


def ext2qt(extent=np.array([-1.8,-1.7,48.4,48.5]),lL0=np.array([-2,48])):
    """ convert an extent region into a list of qt regions

    Parameters
    ----------

    extent : np.array
        (lonmin,lonmax,latmin,latmax)
    lL0 : np.array
        tile origin (lon,lat) of lower left corner

    Returns
    -------

    ltiles : list of sub tiles

    Examples
    --------

    >>> import numpy as np
    >>> from pylayers.gis.gisutil import *
    >>> lL0 = np.array([-2,48])
    >>> extent = np.array([-1.8,-1.7,48.4,48.5])
    >>> ltile = ext2qt(extent,lL0)

    """

    lm = extent[0]
    lM = extent[1]
    Lm = extent[2]
    LM = extent[3]

    lL = np.array([[lm,Lm],[lM,LM]])
    uf8 = ent(lL)
    ill = uf8[0].replace('i','').split('-')
    iur = uf8[1].replace('i','').split('-')
    il  = np.arange(eval(ill[0]),eval(iur[0]))
    iL  = np.arange(eval(ill[1]),eval(iur[1]))
    ltile = []

    for l in il:
        for L in iL:
            ltile.append('i'+str(l)+'-'+str(L))

    return(ltile)


def arr2lp(arr):
    """ convert zeros separated array to list of array

    Returns
    -------

    lp : list or arrays

    """
    lp=[]
    ustart = 0
    sarr = np.sum(arr,axis=1)
    uz = np.where(sarr==0)[0]

    for k in range(len(uz)):
        p = arr[ustart:uz[k]]
        lp.append(p)
        ustart = uz[k]+1
    return(lp)


def ctrad2qt(extent):
    """ convert center,radius into a list of qt regions
    """
    pass


if __name__=='__main__':
    lL0 = np.array([-2,48])
    lL = lL0[:,None] + np.random.rand(2,5)
    print(lL)
    ud16 = eqt(lL)
    un8  = ent(lL)
    lLb  = dqt(ud16,lL0)


def minsec2dec(old):
    """
    convert latlon from DMS (minute second) to DD (decimal)

    Parameters
    ----------

    old : string
        format : DMS format

    Example
    -------

        >>> from pylayers.gis.gisutil import minsec2dec
        >>> minsec2dec('50 03 59 N')
         -50.06638888888889


    """
    direction = {'N':-1, 'S':1, 'E': -1, 'W':1}
    new = old.replace(u'',' ').replace('\'',' ').replace('"',' ')
    new = new.split()
    new_dir = new.pop()
    return (int(new[0])+int(new[1])/60.0+int(new[2])/3600.0) * direction[new_dir]


def distance_on_earth(lat1, long1, lat2, long2):
    """
    Compute great circle distance (the shortest distance over the earths surface)
    between 2 points on earth: A(lat1,lon1) and B(lat2,lon2) 

    (This is a John Cook code, thanks to him !)

    Parameters
    ---------

    lat1 : float for  DDformat | str for DMS format
        latitude first point 
    lat2 : float for  DDformat | str for DMS format
        latitude second point 
    lon1 : float for  DDformat | str for DMS format
        longitude first point 
    lon2 : float for  DDformat | str for DMS format
        longitude second point 

    Return
    ------

    arc : float
        length of the arc on a spherical earth


    Reference
    ---------

    http://www.johndcook.com/blog/python_longitude_latitude/
    http://www.johndcook.com/lat_long_details.html


    """

    if isinstance(lat1,str):
        lat1 = minsec2dec(lat1)
    if isinstance(lat2,str):
        lat2 = minsec2dec(lat2)
    if isinstance(long1,str):
        long1 = minsec2dec(long1)
    if isinstance(long2,str):
        long2 = minsec2dec(long2)

    R = 6371000 # earth radius in meter

    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = np.pi/180.0

    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians

    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians

    # Compute spherical distance from spherical coordinates.

    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length

    cos = (np.sin(phi1)*np.sin(phi2)*np.cos(theta1 - theta2) + 
           np.cos(phi1)*np.cos(phi2))
    arc = np.arccos( cos )

    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return arc*R




# def haversine(lat1, lon1, lat2, lon2):


#     R = 6371000 # earth radius in meter
#     phi1 = lat1 * np.pi/180.
#     phi2 = lat2 * np.pi/180.
#     deltaphi = (lat2-lat1) * np.pi/180.
#     deltalamda = (lon2-lon1) * np.pi/180.

#     a = np.sin(deltaphi/2.) * np.sin(deltaphi/2.) + np.cos(phi1) * np.cos(phi2) * np.sin(deltalamda/2.) * np.sin(deltalamda/2.)
#     c = 2. * np.arctan2(np.sqrt(a), np.sqrt(1.-a))

#     return R * c




def haversine(lat1, lon1, lat2, lon2,mode = 'normal'):
    """
    lat 1 (Na)
    lon 1 (Na)
    lat 2 (Nb)
    lon 2 (Nb)


    mode = 'normal' | 'comb'

        if normal : Na=Nb shape of lat1 & lon1 == lat2 & lon2.
            compute distnace elementwise
        if comb : return a matrix (Na,Nb) of distance. make all combinations

    result(Na,Nb)
    """

    if len(lat1.shape) == 0:
        lat1 = np.array([lat1])
    if len(lon1.shape) == 0:
        lon1 = np.array([lon1])

    if len(lat2.shape) == 0:
        lat2 = np.array([lat2])
    if len(lon2.shape) == 0:
        lon2 = np.array([lon2])

    R = 6371000 # earth radius in meter
    if mode == 'comb':
        phi1 = (lat1 * np.pi/180.)[:,None]
        phi2 = (lat2 * np.pi/180.)[None,:]
        deltaphi = (lat2[None,:]-lat1[:,None]) * np.pi/180.
        deltalamda = (lon2[None,:]-lon1[:,None]) * np.pi/180.
    elif mode == 'normal':
        phi1 = (lat1 * np.pi/180.)
        phi2 = (lat2 * np.pi/180.)
        deltaphi = (lat2-lat1) * np.pi/180.
        deltalamda = (lon2-lon1) * np.pi/180.

    a = np.sin(deltaphi/2.) * np.sin(deltaphi/2.) + np.cos(phi1) * np.cos(phi2) * np.sin(deltalamda/2.) * np.sin(deltalamda/2.)
    c = 2. * np.arctan2(np.sqrt(a), np.sqrt(1.-a))

    return R * c

def get_google_elev_profile(node0,node1,nb_samples=10):
    """
        return elevation profile between 2 nodes
        using google elevation API data

        Parameters
        ----------

            node0 = [lon,lat]
            node1 = [lon,lat]
            nb_samples = int

        Returns
        -------

        profile : np.ndarray (nb_samples)
            elevation for the nb_sambples between node0 and node1



    """
    ELEVATION_BASE_URL = 'https://maps.googleapis.com/maps/api/elevation/json?'
    
    n0  = str(node0[0]) + ',' + str(node0[1])
    n1  = str(node1[0]) + ',' + str(node1[1])
    url = ELEVATION_BASE_URL + 'path=' + n0 + '|' + n1 +'&samples=' + str(nb_samples)
    response = simplejson.load(urllib.request.urlopen(url))
    print(response['status'])
    if response['status'] =='OK':
        profile = []
        for k in range(nb_samples):
            profile.append(response['results'][k]['elevation'])
        return np.array(profile)
    else:
        print ('issue in getting elevation from google')

def get_osm_elev_profile(node0,node1,nb_samples=10):
    """
        return elevation profile between 2 nodes
        using google elevation API data

        Parameters
        ----------

            node0 = [lon,lat]
            node1 = [lon,lat]
            nb_samples = int

        Return
        ------

        profile : np.ndarray (nb_samples)
            elevation for the nb_sambples between node0 and node1



    """


    n0 = np.array(node0)
    n1 = np.array(node1)

    locations = (n1-n0)*np.linspace(0,1,nb_samples)[:,None] + n0
    loc = [{"latitude":k[0],"longitude":k[1]} for k in locations]
    data = {"locations":loc}
    dataj= json.dumps(data)

    URL = 'https://api.open-elevation.com/api/v1/lookup'
    headers = {'Content-type': 'application/json', 'Accept': 'text/plain'}
    r = requests.post(URL, data=dataj,headers=headers)

    if r.status_code == 200:
        pass
    elif r.status_code == 400:
        print('invalid json? ')
        return np.array([])
    else:
        print('issue getting elevation info')
        return np.array([])
    dres=json.loads(r.text)['results']

    # res = np.array([[x['latitude'],x['longitude'],x['elevation']] for x in delev])
    elev = np.array([x['elevation'] for x in dres])
    return elev
