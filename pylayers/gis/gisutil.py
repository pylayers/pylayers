# -*- coding:Utf-8 -*-

import numpy as np
import pdb

def ent(lL,lL0):
    """ encode lon Lat in natural integer

    Parameters
    ----------

    lL : nd.array (2xN)
        longitude latitude array
    lL0 : nd.array (,2)
        lower left corner of the 1 degree tile


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
    >>> lL0 = np.array([-2,48])
    >>> ent(lL,lL0)

    """
    # offset from the lower left corner
    d = lL-lL0
    dui8 = np.floor(d*255).astype('uint8')
    lab = map(lambda x:'i'+str(x[0])+'-'+str(x[1]),dui8)
    return(lab)


def eqt(lL,lL0):
    """ encode lon Lat in quad tree integer

    Parameters
    ----------

    lL : nd.array (2xN)
        longitude Latitude
    lL0 : nd.array (,2)
        lower left corner of the 1degree tile

    """
    N = lL.shape[1]
    # offset from the lower left corner
    d = lL-lL0[:,None]
    dui8 = np.floor(d*256).astype('uint8')
    ndui8 = np.unpackbits(dui8).reshape(2,N,8)
    d16 = np.empty((N,16)).astype('int')
    d16[:,1::2]=ndui8[0,:,:]
    d16[:,0::2]=ndui8[1,:,:]
    ud8 = np.packbits(d16,axis=1)
    ud16 = (ud8[:,0]*256+ud8[:,1]).astype('uint16')
    return(ud16 )


def dqt(ud16,lL0):
    """ decode quad tree integer to lon Lat

    Parameters
    ----------

    lL : nd.array (2xN)
        longitude Latitude
    lL0 : nd.array (,2)
        lower left corner of the 1degree tile

    """

    N = len(ud16)
    # offset from the lower left corner
    #d = lL-lL0[:,None]
    #dui8 = np.floor(d*256).astype('uint8')
    uh8  = ud16/256
    ul8  = ud16-uh8*256
    ud8 =  (np.vstack((uh8,ul8)).T).astype('uint8')
    ud16 = np.unpackbits(ud8).reshape(N,16)
    ndu8 = np.empty((2,N,8)).astype('int')
    ndu8[0,:,:]=ud16[:,1::2]
    ndu8[1,:,:]=ud16[:,0::2]
    du8 = np.packbits(ndu8).reshape(2,N)/256.
    lL = lL0[:,None]+du8

    return(lL)


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
    uf8 = ent(lL,lL0)
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
    print lL
    ud16 = eqt(lL,lL0)
    un8  = ent(lL,lL0)
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




def haversine(lat1, lon1, lat2, lon2):


    R = 6371000 # earth radius in meter
    phi1 = lat1 * np.pi/180.
    phi2 = lat2 * np.pi/180.
    deltaphi = (lat2-lat1) * np.pi/180.
    deltalamda = (lon2-lon1) * np.pi/180.

    a = np.sin(deltaphi/2) * np.sin(deltaphi/2) + np.cos(phi1) * np.cos(phi2) * np.sin(deltalamda/2) * np.sin(deltalamda/2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))

    return R * c
