# -*- coding:Utf-8 -*-
import os
import numpy as np 
import scipy as sp 
import matplotlib.pylab as plt 
import scipy.special as spe
import doctest
import logging
from   bitstring  import BitString
import datetime as dat
from pylayers.util.project import *
import shutil
import sys
import zipfile
#
# getlong 
# getshort
# getdir
# shp
# dimcmp
# tstincl
# findpos
# ininter
#
###################################
#
# Wave Related functions
#
##############
def delay(p1,p2):
    """ delay in ns between 2 point 

    Parameters
    ----------

    p1  ndarray (1x2)
        point 1 coordinates 
    p2  ndarray (1x2)
        point 2 coordinates

    Examples
    --------

    >>> p1 = np.array([0,0])
    >>> p2 = np.array([0,0.3])
    >>> tau = delay(p1,p2)
    >>> assert tau==1.,"Warning : speed of light has changed"

    """
    v   = p1-p2
    d2  = np.dot(v,v)
    d   = np.sqrt(d2)
    tau = d/0.3
    return(tau)
#####################
def lt2dic(lt):
    """ convert list of tuple to dictionary 

    Parameters
    ----------

    lt : list of tuple
        the first element of tuple is th ekey dictionary

    Examples
    --------

       >>> lt = [ ('1','1 2 3'),('2','1.5 2'),('3','4.78 89.0 2')]
       >>> d = lt2dic(lt)
    """
    dic = {}
    for tup in lt:
        dic[tup[0]]=tup[1]
    return(dic)    

def lt2idic(lt):
    """ convert list of tuple to dictionary 

    Parameters
    ----------

    lt : list 

    Examples 
    --------

    >>> from pylayers.util.pyutil import *
    >>> lt = [ ('1','1 2 3'),('2','1.5 2 3'),('3','4.78 89.0 2')]
    >>> d = lt2idic(lt)

    """
    dic = {}
    for tup in lt:
        val = tup[1].split(' ')
        dic[int(tup[0])]=np.array([float(val[0]),float(val[1]),float(val[2])])
    return(dic)    

def getlong(shortname,directory):
    """  get a long name

    This function allows to construct the long file name relatively
    to a current project directory which is stored in the environment
    variable $BASENAME

    Parameters
    ----------
    shortname : string
        short name of the file

    dir       : directory in $BASENAME or $PYLAYERS
    source    : string 

    """
    try:
        basename=os.environ['BASENAME']
    except:
        logging.critical("BASENAME environment variable should be defined")
        #basename=os.environ['HOME']+"/Pyproject"

    longname = basename+'/'+directory+'/'+shortname
    return(longname)

def getshort(longname):
    shortname=os.path.split(longname)[1]
    return(shortname)

def getdir(longname):
    rac=os.path.split(longname)[0]
    dir=os.path.split(rac)[1]
    return(dir)
        
def shp(arr):
    """ return dimension of an array
    Parameters
    ----------
    arr : ndarray

    Returns 
    -------
    shp : tuple 
    
    Examples
    --------

    >>> import pylayers.util.pyutil as pyu 
    >>> import numpy as np 
    >>> from scipy import *
    >>> a = np.arange(10)
    >>> pyu.shp(a)
    (1, 10)
    >>> b = randn(2,2)
    >>> pyu.shp(b)
    (2, 2)

    """
    ndim = arr.ndim
    if ndim>1:
        shp = np.shape(arr)
    else:
        shp = (1,len(arr))
    return(shp)             

def dimcmp(ar1,ar2):
    """ compare shape of arrays  
    
    Parameters
    ----------
    ar1 : ndarray
    ar2 : ndarray

    Returns
    -------
    return code : int 
        0 arrays are not compatible
        1 arrayis have same dimension 
        2 second argument has greater dimension 
        3 first argument has greater dimension
    """
    sh1 = shp(ar1)
    sh2 = shp(ar2)
    if (sh1[0]==sh2[0]):
        return(1)
    if ((sh1[0]!=1)&(sh2[0]!=1)):
        return(0)
    if (sh2[0]>sh1[0]):
        return(2)
    else:
        return(3)

def tstincl(ar1,ar2):
    """ testincl(ar1,ar2)

    Parameters
    ----------
    ar1 : ndarray
    ar2 : ndarray

    Returns
    -------
    0 : if ar1 and ar2 have no points in common
    1 : if  ar2 includes ar1
    2 : else

    Summary
    -------
        test wheteher the ar1 interval is included in interval ar2  

    """
    if ((ar1[0]>=ar2[0])&(ar1[-1]<=ar2[-1])):
        return(1)
    if ((ar1[0]>ar2[-1]) or (ar2[0]>ar1[-1])):
        return(0)
    else:
        return(2)       

def findpos(ar,val):
    """
    findpos(ar,val)

    return the i position in array ar, such that ar[i] == val 
    if the value is not find, return 'value not found'
    """
    dim=len(ar)
    i=0
    find=0
    while (i<dim):
        if (ar[i]==val):
            pos=i
            find=1
            i=dim
        else:
            i=i+1
    if (find==1):
        return(pos)
    else:
        return('value not found')

def ininter(ar,val1,val2):
    """
    ininter(ar,val1,val2):

    This function return the set of samples from array ar 
    which are included in the interval [val1 val2] 

    Usage Case :

    """
    criterium= (ar>=val1)&(ar<=val2)
    return(ar[criterium])

def cshift(l, offset):
    """ ndarray circular shift     
        Parameters
         ----------
      l       : ndarray 
      offset  : int
      
      The offset value can be either positive or negative and the applied
      offset value is applied modulo the length of l
          

      >>> a = np.array([1,2,3])
      >>> b = cshift(a,1)
      >>> c = cshift(a,-1)
      >>> d = cshift(a,4)
    """
    offset %= len(l)
    return np.concatenate((l[-offset:], l[:-offset]))

def LegFunc(nn,ntrunc,theta,phi):
    """ Compute Legendre functions Ylm(theta,phi)

    Parameters
    ----------
    nn        :  integer 
    ntrunc    :  integer 
    theta     :  np.array(1xNtheta)
    theta     :  np.array(1xNtheta)
    phi       :  np.array(1xNtheta)
    """
    m=array(zeros(nn),dtype=integer)
    l=array(zeros(nn),dtype=integer)
    val=r_[0:ntrunc+1:1]
    k=0
    pas=ntrunc+1
    start=0
    stop=0
    while (stop<nn):
        stop=start+pas
        m[start:stop]=val[k]
        l[start:stop]=r_[k:ntrunc+1:1]
        k=k+1
        start=stop
        pas=pas-1
    Ylm=[]
    for i in range(nn):
        ylm = sph_harm(m[i],l[i],phi,theta)
        Ylm.append(ylm)
    Ylm=array(Ylm)
    return(Ylm)

def ExpFunc (x,y):
    """ Exponential fitting
    Parameters
    ----------
    x : np.array
    y : np.array
    
    Returns
    -------
    a : estimation of \\alpha
    b : estimation of \\beta

    Notes
    -----
    find a possible equation of an exponential function of the form :

        .. math::    y = \\alpha e^{- \\beta x}    

    Examples
    --------

    >>> a = 3
    >>> b = 2
    >>> x = sp.rand(100)
    >>> n = 0.3*sp.randn(100)
    >>> y  = a*np.exp(-b*x) +  abs(n)
    >>> alpha,beta = ExpFunc(x,y)

    """
    z = np.log(y)

    (a,b) = sp.polyfit(x,z,1)

    z2 = sp.polyval([a,b],x)

    alpha = np.exp(b)
    beta = -a

    return(alpha,beta)

def InvFunc (x,z):
    """
    (alfa,beta) = InvFunc (x,y)

    find a possible equation of an inverse proportional function of the form :
         y = alfa/x + beta    
    """
    y = 1./x
    (a,b)=polyfit(y,z,1)

    return(a,b)

def PowFunc (x,y):
    """
    (alfa,beta) = PowFunc (x,y)

    find a possible equation of a power function of the form :
         y = alfa/ x**beta    
    """
    t = 1./x
    z=log(y)
    u=log(t)
    (a,b)=polyfit(u,z,1)

    beta=a
    alpha=exp(b)

    return(alpha,beta)

def randcol(Nc):
    """ get random color 

    Parameters 
    -----------
    Nc :  int 
        Number of color 

    Returns
    -------
    col : list
       A list of colors.

    Example
    -------
    >>> from pylayers.util.pyutil import *
    >>> import matplotlib.pyplot as plt
    >>> col = randcol(100)
    """
    col=[]
    lin=np.linspace(255,16777215,Nc)

    for i in range(Nc):

        hexa=hex(lin[i])
        if hexa[-1] == 'L':
            lh=len(hexa[2:-1])
            hexa='#' +'0'*(6-lh) + hexa[2:-1]
        elif len(hexa)<8:
             hexa='#' +'0'*(6-len(hexa)) +hexa[2:]
        col.append(hexa[0:7])  
    return(col)


def coldict(): 
    """ Color dictionary 
    html color 

    Notes 
    -----

    'Link on html color<http://html-color-codes.blogspot.com/>'_

    """
    cold={}
    cold['black']=  '#000000'
    cold['grey']=  '#BEBEBE'
    cold['DimGrey']=  '#696969'
    cold['LightGray']=  '#D3D3D3'
    cold['LightSlateGrey']=  '#778899'
    cold['SlateGray']=  '#708090'
    cold['SlateGray1']=  '#C6E2FF'
    cold['SlateGray2']=  '#B9D3EE'
    cold['SlateGray3']=  '#9FB6CD'
    cold['SlateGray4']=  '#6C7B8B'
    cold['SlateGrey']=  '#708090'
    cold['grey0']=  '#000000'
    cold['grey1']=  '#030303'
    cold['grey2']=  '#050505'
    cold['grey3']=  '#080808'
    cold['grey4']=  '#0A0A0A'
    cold['grey5']=  '#0D0D0D'
    cold['grey6']=  '#0F0F0F'
    cold['grey7']=  '#121212'
    cold['grey8']=  '#141414'
    cold['grey9']=  '#171717'
    cold['grey10']=  '#1A1A1A'
    cold['grey11']=  '#1C1C1C'
    cold['grey12']=  '#1F1F1F'
    cold['grey13']=  '#212121'
    cold['grey14']=  '#242424'
    cold['grey15']=  '#262626'
    cold['grey16']=  '#292929'
    cold['grey17']=  '#2B2B2B'
    cold['grey18']=  '#2E2E2E'
    cold['grey19']=  '#303030'
    cold['grey20']=  '#333333'
    cold['grey21']=  '#363636'
    cold['grey22']=  '#383838'
    cold['grey23']=  '#3B3B3B'
    cold['grey24']=  '#3D3D3D'
    cold['grey25']=  '#404040'
    cold['grey26']=  '#424242'
    cold['grey27']=  '#454545'
    cold['grey28']=  '#474747'
    cold['grey29']=  '#4A4A4A'
    cold['grey30']=  '#4D4D4D'
    cold['grey31']=  '#4F4F4F'
    cold['grey32']=  '#525252'
    cold['grey33']=  '#545454'
    cold['grey34']=  '#575757'
    cold['grey35']=  '#595959'
    cold['grey36']=  '#5C5C5C'
    cold['grey37']=  '#5E5E5E'
    cold['grey38']=  '#616161'
    cold['grey39']=  '#636363'
    cold['grey40']=  '#666666'
    cold['grey41']=  '#696969'
    cold['grey42']=  '#6B6B6B'
    cold['grey43']=  '#6E6E6E'
    cold['grey44']=  '#707070'
    cold['grey45']=  '#737373'
    cold['grey46']=  '#757575'
    cold['grey47']=  '#787878'
    cold['grey48']=  '#7A7A7A'
    cold['grey49']=  '#7D7D7D'
    cold['grey50']=  '#7F7F7F'
    cold['grey51']=  '#828282'
    cold['grey52']=  '#858585'
    cold['grey53']=  '#878787'
    cold['grey54']=  '#8A8A8A'
    cold['grey55']=  '#8C8C8C'
    cold['grey56']=  '#8F8F8F'
    cold['grey57']=  '#919191'
    cold['grey58']=  '#949494'
    cold['grey59']=  '#969696'
    cold['grey60']=  '#999999'
    cold['grey61']=  '#9C9C9C'
    cold['grey62']=  '#9E9E9E'
    cold['grey63']=  '#A1A1A1'
    cold['grey64']=  '#A3A3A3'
    cold['grey65']=  '#A6A6A6'
    cold['grey66']=  '#A8A8A8'
    cold['grey67']=  '#ABABAB'
    cold['grey68']=  '#ADADAD'
    cold['grey69']=  '#B0B0B0'
    cold['grey70']=  '#B3B3B3'
    cold['grey71']=  '#B5B5B5'
    cold['grey72']=  '#B8B8B8'
    cold['grey73']=  '#BABABA'
    cold['grey74']=  '#BDBDBD'
    cold['grey75']=  '#BFBFBF'
    cold['grey76']=  '#C2C2C2'
    cold['grey77']=  '#C4C4C4'
    cold['grey78']=  '#C7C7C7'
    cold['grey79']=  '#C9C9C9'
    cold['grey80']=  '#CCCCCC'
    cold['grey81']=  '#CFCFCF'
    cold['grey82']=  '#D1D1D1'
    cold['grey83']=  '#D4D4D4'
    cold['grey84']=  '#D6D6D6'
    cold['grey85']=  '#D9D9D9'
    cold['grey86']=  '#DBDBDB'
    cold['grey87']=  '#DEDEDE'
    cold['grey88']=  '#E0E0E0'
    cold['grey89']=  '#E3E3E3'
    cold['grey90']=  '#E5E5E5'
    cold['grey91']=  '#E8E8E8'
    cold['grey92']=  '#EBEBEB'
    cold['grey93']=  '#EDEDED'
    cold['grey94']=  '#F0F0F0'
    cold['grey95']=  '#F2F2F2'
    cold['grey96']=  '#F5F5F5'
    cold['grey97']=  '#F7F7F7'
    cold['grey98']=  '#FAFAFA'
    cold['grey99']=  '#FCFCFC'
    cold['grey100']=  '#FFFFFF'
    cold['AliceBlue']=  '#F0F8FF'
    cold['BlueViolet']=  '#8A2BE2'
    cold['CadetBlue']=  '#5F9EA0'
    cold['CadetBlue1']=  '#98F5FF'
    cold['CadetBlue2']=  '#8EE5EE'
    cold['CadetBlue3']=  '#7AC5CD'
    cold['CadetBlue4']=  '#53868B'
    cold['CornflowerBlue']=  '#6495ED'
    cold['DarkSlateBlue']=  '#483D8B'
    cold['DarkTurquoise']=  '#00CED1'
    cold['DeepSkyBlue']=  '#00BFFF'
    cold['DeepSkyBlue1']=  '#00BFFF'
    cold['DeepSkyBlue2']=  '#00B2EE'
    cold['DeepSkyBlue3']=  '#009ACD'
    cold['DeepSkyBlue4']=  '#00688B'
    cold['DodgerBlue']=  '#1E90FF'
    cold['DodgerBlue1']=  '#1E90FF'
    cold['DodgerBlue2']=  '#1C86EE'
    cold['DodgerBlue3']=  '#1874CD'
    cold['DodgerBlue4']=  '#104E8B'
    cold['LightBlue']=  '#ADD8E6'
    cold['LightBlue1']=  '#BFEFFF'
    cold['LightBlue2']=  '#B2DFEE'
    cold['LightBlue3']=  '#9AC0CD'
    cold['LightBlue4']=  '#68838B'
    cold['LightCyan']=  '#E0FFFF'
    cold['LightCyan1']=  '#E0FFFF'
    cold['LightCyan2']=  '#D1EEEE'
    cold['LightCyan3']=  '#B4CDCD'
    cold['LightCyan4']=  '#7A8B8B'
    cold['LightSkyBlue']=  '#87CEFA'
    cold['LightSkyBlue1']=  '#B0E2FF'
    cold['LightSkyBlue2']=  '#A4D3EE'
    cold['LightSkyBlue3']=  '#8DB6CD'
    cold['LightSkyBlue4']=  '#607B8B'
    cold['LightSlateBlue']=  '#8470FF'
    cold['LightSteelBlue']=  '#B0C4DE'
    cold['LightSteelBlue1']=  '#CAE1FF'
    cold['LightSteelBlue2']=  '#BCD2EE'
    cold['LightSteelBlue3']=  '#A2B5CD'
    cold['LightSteelBlue4']=  '#6E7B8B'
    cold['MediumAquararine']=  '#66CDAA'
    cold['MediumBlue']=  '#0000CD'
    cold['MediumSlateBlue']=  '#7B68EE'
    cold['MediumTurquoise']=  '#48D1CC'
    cold['MidnightBlue']=  '#191970'
    cold['NavyBlue']=  '#000080'
    cold['PaleTurquoise']=  '#AFEEEE'
    cold['PaleTurquoise1']=  '#BBFFFF'
    cold['PaleTurquoise2']=  '#AEEEEE'
    cold['PaleTurquoise3']=  '#96CDCD'
    cold['PaleTurquoise4']=  '#668B8B'
    cold['PowderBlue']=  '#B0E0E6'
    cold['RoyalBlue']=  '#4169E1'
    cold['RoyalBlue1']=  '#4876FF'
    cold['RoyalBlue2']=  '#436EEE'
    cold['RoyalBlue3']=  '#3A5FCD'
    cold['RoyalBlue4']=  '#27408B'
    cold['RoyalBlue5']=  '#002266'
    cold['SkyBlue']=  '#87CEEB'
    cold['SkyBlue1']=  '#87CEFF'
    cold['SkyBlue2']=  '#7EC0EE'
    cold['SkyBlue3']=  '#6CA6CD'
    cold['SkyBlue4']=  '#4A708B'
    cold['SlateBlue']=  '#6A5ACD'
    cold['SlateBlue1']=  '#836FFF'
    cold['SlateBlue2']=  '#7A67EE'
    cold['SlateBlue3']=  '#6959CD'
    cold['SlateBlue4']=  '#473C8B'
    cold['SteelBlue']=  '#4682B4'
    cold['SteelBlue1']=  '#63B8FF'
    cold['SteelBlue2']=  '#5CACEE'
    cold['SteelBlue3']=  '#4F94CD'
    cold['SteelBlue4']=  '#36648B'
    cold['aquamarine']=  '#7FFFD4'
    cold['aquamarine1']=  '#7FFFD4'
    cold['aquamarine2']=  '#76EEC6'
    cold['aquamarine3']=  '#66CDAA'
    cold['aquamarine4']=  '#458B74'
    cold['azure']=  '#F0FFFF'
    cold['azure1']=  '#F0FFFF'
    cold['azure2']=  '#E0EEEE'
    cold['azure3']=  '#C1CDCD'
    cold['azure4']=  '#838B8B'
    cold['blue']=  '#0000FF'
    cold['blue1']=  '#0000FF'
    cold['blue2']=  '#0000EE'
    cold['blue3']=  '#0000CD'
    cold['blue4']=  '#00008B'
    cold['cyan']=  '#00FFFF'
    cold['cyan1']=  '#00FFFF'
    cold['cyan2']=  '#00EEEE'
    cold['cyan3']=  '#00CDCD'
    cold['cyan4']=  '#008B8B'
    cold['navy']=  '#000080'
    cold['turquoise']=  '#40E0D0'
    cold['turquoise1']=  '#00F5FF'
    cold['turquoise2']=  '#00E5EE'
    cold['turquoise3']=  '#00C5CD'
    cold['turquoise4']=  '#00868B'
    cold['DarkSlateGray']=  '#2F4F4F'
    cold['DarkSlateGray1']=  '#97FFFF'
    cold['DarkSlateGray2']=  '#8DEEEE'
    cold['DarkSlateGray3']=  '#79CDCD'
    cold['DarkSlateGray4']=  '#528B8B'
    cold['RosyBrown']=  '#BC8F8F'
    cold['RosyBrown1']=  '#FFC1C1'
    cold['RosyBrown2']=  '#EEB4B4'
    cold['RosyBrown3']=  '#CD9B9B'
    cold['RosyBrown4']=  '#8B6969'
    cold['SaddleBrown']=  '#8B4513'
    cold['SandyBrown']=  '#F4A460'
    cold['beige']=  '#F5F5DC'
    cold['brown']=  '#A52A2A'
    cold['brown1']=  '#FF4040'
    cold['brown2']=  '#EE3B3B'
    cold['brown3']=  '#CD3333'
    cold['brown4']=  '#8B2323'
    cold['burlywood']=  '#DEB887'
    cold['burlywood1']=  '#FFD39B'
    cold['burlywood2']=  '#EEC591'
    cold['burlywood3']=  '#CDAA7D'
    cold['burlywood4']=  '#8B7355'
    cold['chocolate']=  '#D2691E'
    cold['chocolate1']=  '#FF7F24'
    cold['chocolate2']=  '#EE7621'
    cold['chocolate3']=  '#CD661D'
    cold['chocolate4']=  '#8B4513'
    cold['peru']=  '#CD853F'
    cold['tan']=  '#D2B48C'
    cold['tan1']=  '#FFA54F'
    cold['tan2']=  '#EE9A49'
    cold['tan3']=  '#CD853F'
    cold['tan4']=  '#8B5A2B'
    cold['DarkGreen']=  '#006400'
    cold['DarkKhaki']=  '#BDB76B'
    cold['DarkOliveGreen']=  '#556B2F'
    cold['DarkOliveGreen1']=  '#CAFF70'
    cold['DarkOliveGreen2']=  '#BCEE68'
    cold['DarkOliveGreen3']=  '#A2CD5A'
    cold['DarkOliveGreen4']=  '#6E8B3D'
    cold['DarkSeaGreen']=  '#8FBC8F'
    cold['DarkSeaGreen1']=  '#C1FFC1'
    cold['DarkSeaGreen2']=  '#B4EEB4'
    cold['DarkSeaGreen3']=  '#9BCD9B'
    cold['DarkSeaGreen4']=  '#698B69'
    cold['ForestGreen']=  '#228B22'
    cold['GreenYellow']=  '#ADFF2F'
    cold['LawnGreen']=  '#7CFC00'
    cold['LightSeaGreen']=  '#20B2AA'
    cold['LimeGreen']=  '#32CD32'
    cold['MediumSeaGreen']=  '#3CB371'
    cold['MediumSpringGreen']=  '#00FA9A'
    cold['MintCream']=  '#F5FFFA'
    cold['OliveDrab']=  '#6B8E23'
    cold['OliveDrab1']=  '#C0FF3E'
    cold['OliveDrab2']=  '#B3EE3A'
    cold['OliveDrab3']=  '#9ACD32'
    cold['OliveDrab4']=  '#698B22'
    cold['PaleGreen']=  '#98FB98'
    cold['PaleGreen1']=  '#9AFF9A'
    cold['PaleGreen2']=  '#90EE90'
    cold['PaleGreen3']=  '#7CCD7C'
    cold['PaleGreen4']=  '#548B54'
    cold['SeaGreen']=  '#2E8B57'
    cold['SeaGreen1']=  '#54FF9F'
    cold['SeaGreen2']=  '#4EEE94'
    cold['SeaGreen3']=  '#43CD80'
    cold['SeaGreen4']=  '#2E8B57'
    cold['SpringGreen']=  '#00FF7F'
    cold['SpringGreen1']=  '#00FF7F'
    cold['SpringGreen2']=  '#00EE76'
    cold['SpringGreen3']=  '#00CD66'
    cold['SpringGreen4']=  '#008B45'
    cold['YellowGreen']=  '#9ACD32'
    cold['chartreuse']=  '#7FFF00'
    cold['chartreuse1']=  '#7FFF00'
    cold['chartreuse2']=  '#76EE00'
    cold['chartreuse3']=  '#66CD00'
    cold['chartreuse4']=  '#458B00'
    cold['green']=  '#00FF00'
    cold['green1']=  '#00FF00'
    cold['green2']=  '#00EE00'
    cold['green3']=  '#00CD00'
    cold['green4']=  '#008B00'
    cold['khaki']=  '#F0E68C'
    cold['khaki1']=  '#FFF68F'
    cold['khaki2']=  '#EEE685'
    cold['khaki3']=  '#CDC673'
    cold['khaki4']=  '#8B864E'
    cold['DarkOrange']=  '#FF8C00'
    cold['DarkOrange1']=  '#FF7F00'
    cold['DarkOrange2']=  '#EE7600'
    cold['DarkOrange3']=  '#CD6600'
    cold['DarkOrange4']=  '#8B4500'
    cold['DarkSalmon']=  '#E9967A'
    cold['LightCoral']=  '#F08080'
    cold['LightSalmon']=  '#FFA07A'
    cold['LightSalmon1']=  '#FFA07A'
    cold['LightSalmon2']=  '#EE9572'
    cold['LightSalmon3']=  '#CD8162'
    cold['LightSalmon4']=  '#8B5742'
    cold['PeachPuff']=  '#FFDAB9'
    cold['PeachPuff1']=  '#FFDAB9'
    cold['PeachPuff2']=  '#EECBAD'
    cold['PeachPuff3']=  '#CDAF95'
    cold['PeachPuff4']=  '#8B7765'
    cold['bisque']=  '#FFE4C4'
    cold['bisque1']=  '#FFE4C4'
    cold['bisque2']=  '#EED5B7'
    cold['bisque3']=  '#CDB79E'
    cold['bisque4']=  '#8B7D6B'
    cold['coral']=  '#FF7F50'
    cold['coral1']=  '#FF7256'
    cold['coral2']=  '#EE6A50'
    cold['coral3']=  '#CD5B45'
    cold['coral4']=  '#8B3E2F'
    cold['honeydew']=  '#F0FFF0'
    cold['honeydew1']=  '#F0FFF0'
    cold['honeydew2']=  '#E0EEE0'
    cold['honeydew3']=  '#C1CDC1'
    cold['honeydew4']=  '#838B83'
    cold['orange']=  '#FFA500'
    cold['orange1']=  '#FFA500'
    cold['orange2']=  '#EE9A00'
    cold['orange3']=  '#CD8500'
    cold['orange4']=  '#8B5A00'
    cold['salmon']=  '#FA8072'
    cold['salmon1']=  '#FF8C69'
    cold['salmon2']=  '#EE8262'
    cold['salmon3']=  '#CD7054'
    cold['salmon4']=  '#8B4C39'
    cold['sienna']=  '#A0522D'
    cold['sienna1']=  '#FF8247'
    cold['sienna2']=  '#EE7942'
    cold['sienna3']=  '#CD6839'
    cold['sienna4']=  '#8B4726'
    cold['DeepPink']=  '#FF1493'
    cold['DeepPink1']=  '#FF1493'
    cold['DeepPink2']=  '#EE1289'
    cold['DeepPink3']=  '#CD1076'
    cold['DeepPink4']=  '#8B0A50'
    cold['HotPink']=  '#FF69B4'
    cold['HotPink1']=  '#FF6EB4'
    cold['HotPink2']=  '#EE6AA7'
    cold['HotPink3']=  '#CD6090'
    cold['HotPink4']=  '#8B3A62'
    cold['IndianRed']=  '#CD5C5C'
    cold['IndianRed1']=  '#FF6A6A'
    cold['IndianRed2']=  '#EE6363'
    cold['IndianRed3']=  '#CD5555'
    cold['IndianRed4']=  '#8B3A3A'
    cold['LightPink']=  '#FFB6C1'
    cold['LightPink1']=  '#FFAEB9'
    cold['LightPink2']=  '#EEA2AD'
    cold['LightPink3']=  '#CD8C95'
    cold['LightPink4']=  '#8B5F65'
    cold['MediumVioletRed']=  '#C71585'
    cold['MistyRose']=  '#FFE4E1'
    cold['MistyRose1']=  '#FFE4E1'
    cold['MistyRose2']=  '#EED5D2'
    cold['MistyRose3']=  '#CDB7B5'
    cold['MistyRose4']=  '#8B7D7B'
    cold['OrangeRed']=  '#FF4500'
    cold['OrangeRed1']=  '#FF4500'
    cold['OrangeRed2']=  '#EE4000'
    cold['OrangeRed3']=  '#CD3700'
    cold['OrangeRed4']=  '#8B2500'
    cold['PaleVioletRed']=  '#DB7093'
    cold['PaleVioletRed1']=  '#FF82AB'
    cold['PaleVioletRed2']=  '#EE799F'
    cold['PaleVioletRed3']=  '#CD6889'
    cold['PaleVioletRed4']=  '#8B475D'
    cold['VioletRed']=  '#D02090'
    cold['VioletRed1']=  '#FF3E96'
    cold['VioletRed2']=  '#EE3A8C'
    cold['VioletRed3']=  '#CD3278'
    cold['VioletRed4']=  '#8B2252'
    cold['firebrick']=  '#B22222'
    cold['firebrick1']=  '#FF3030'
    cold['firebrick2']=  '#EE2C2C'
    cold['firebrick3']=  '#CD2626'
    cold['firebrick4']=  '#8B1A1A'
    cold['pink']=  '#FFC0CB'
    cold['pink1']=  '#FFB5C5'
    cold['pink2']=  '#EEA9B8'
    cold['pink3']=  '#CD919E'
    cold['pink4']=  '#8B636C'
    cold['red']=  '#FF0000'
    cold['red1']=  '#FF0000'
    cold['red2']=  '#EE0000'
    cold['red3']=  '#CD0000'
    cold['red4']=  '#8B0000'
    cold['tomato']=  '#FF6347'
    cold['tomato1']=  '#FF6347'
    cold['tomato2']=  '#EE5C42'
    cold['tomato3']=  '#CD4F39'
    cold['tomato4']=  '#8B3626'
    cold['DarkOrchid']=  '#9932CC'
    cold['DarkOrchid1']=  '#BF3EFF'
    cold['DarkOrchid2']=  '#B23AEE'
    cold['DarkOrchid3']=  '#9A32CD'
    cold['DarkOrchid4']=  '#68228B'
    cold['DarkViolet']=  '#9400D3'
    cold['LavenderBlush']=  '#FFF0F5'
    cold['LavenderBlush1']=  '#FFF0F5'
    cold['LavenderBlush2']=  '#EEE0E5'
    cold['LavenderBlush3']=  '#CDC1C5'
    cold['LavenderBlush4']=  '#8B8386'
    cold['MediumOrchid']=  '#BA55D3'
    cold['MediumOrchid1']=  '#E066FF'
    cold['MediumOrchid2']=  '#D15FEE'
    cold['MediumOrchid3']=  '#B452CD'
    cold['MediumOrchid4']=  '#7A378B'
    cold['MediumPurple']=  '#9370DB'
    cold['MediumPurple1']=  '#AB82FF'
    cold['MediumPurple2']=  '#9F79EE'
    cold['MediumPurple3']=  '#8968CD'
    cold['MediumPurple4']=  '#5D478B'
    cold['lavender']=  '#E6E6FA'
    cold['magenta']=  '#FF00FF'
    cold['magenta1']=  '#FF00FF'
    cold['magenta2']=  '#EE00EE'
    cold['magenta3']=  '#CD00CD'
    cold['magenta4']=  '#8B008B'
    cold['maroon']=  '#B03060'
    cold['maroon1']=  '#FF34B3'
    cold['maroon2']=  '#EE30A7'
    cold['maroon3']=  '#CD2990'
    cold['maroon4']=  '#8B1C62'
    cold['orchid']=  '#DA70D6'
    cold['orchid1']=  '#FF83FA'
    cold['orchid2']=  '#EE7AE9'
    cold['orchid3']=  '#CD69C9'
    cold['orchid4']=  '#8B4789'
    cold['plum']=  '#DDA0DD'
    cold['plum1']=  '#FFBBFF'
    cold['plum2']=  '#EEAEEE'
    cold['plum3']=  '#CD96CD'
    cold['plum4']=  '#8B668B'
    cold['purple']=  '#A020F0'
    cold['purple1']=  '#9B30FF'
    cold['purple2']=  '#912CEE'
    cold['purple3']=  '#7D26CD'
    cold['purple4']=  '#551A8B'
    cold['thistle']=  '#D8BFD8'
    cold['thistle1']=  '#FFE1FF'
    cold['thistle2']=  '#EED2EE'
    cold['thistle3']=  '#CDB5CD'
    cold['thistle4']=  '#8B7B8B'
    cold['violet']=  '#EE82EE'
    cold['AntiqueWhite']=  '#FAEBD7'
    cold['AntiqueWhite1']=  '#FFEFDB'
    cold['AntiqueWhite2']=  '#EEDFCC'
    cold['AntiqueWhite3']=  '#CDC0B0'
    cold['AntiqueWhite4']=  '#8B8378'
    cold['FloralWhite']=  '#FFFAF0'
    cold['GhostWhite']=  '#F8F8FF'
    cold['NavajoWhite']=  '#FFDEAD'
    cold['NavajoWhite1']=  '#FFDEAD'
    cold['NavajoWhite2']=  '#EECFA1'
    cold['NavajoWhite3']=  '#CDB38B'
    cold['NavajoWhite4']=  '#8B795E'
    cold['OldLace']=  '#FDF5E6'
    cold['WhiteSmoke']=  '#F5F5F5'
    cold['gainsboro']=  '#DCDCDC'
    cold['ivory']=  '#FFFFF0'
    cold['ivory1']=  '#FFFFF0'
    cold['ivory2']=  '#EEEEE0'
    cold['ivory3']=  '#CDCDC1'
    cold['ivory4']=  '#8B8B83'
    cold['linen']=  '#FAF0E6'
    cold['seashell']=  '#FFF5EE'
    cold['seashell1']=  '#FFF5EE'
    cold['seashell2']=  '#EEE5DE'
    cold['seashell3']=  '#CDC5BF'
    cold['seashell4']=  '#8B8682'
    cold['snow']=  '#FFFAFA'
    cold['snow1']=  '#FFFAFA'
    cold['snow2']=  '#EEE9E9'
    cold['snow3']=  '#CDC9C9'
    cold['snow4']=  '#8B8989'
    cold['wheat']=  '#F5DEB3'
    cold['wheat1']=  '#FFE7BA'
    cold['wheat2']=  '#EED8AE'
    cold['wheat3']=  '#CDBA96'
    cold['wheat4']=  '#8B7E66'
    cold['white']=  '#FFFFFF'
    cold['BlanchedAlmond']=  '#FFEBCD'
    cold['DarkGoldenrod']=  '#B8860B'
    cold['DarkGoldenrod1']=  '#FFB90F'
    cold['DarkGoldenrod2']=  '#EEAD0E'
    cold['DarkGoldenrod3']=  '#CD950C'
    cold['DarkGoldenrod4']=  '#8B6508'
    cold['LemonChiffon']=  '#FFFACD'
    cold['LemonChiffon1']=  '#FFFACD'
    cold['LemonChiffon2']=  '#EEE9BF'
    cold['LemonChiffon3']=  '#CDC9A5'
    cold['LemonChiffon4']=  '#8B8970'
    cold['LightGoldenrod']=  '#EEDD82'
    cold['LightGoldenrod1']=  '#FFEC8B'
    cold['LightGoldenrod2']=  '#EEDC82'
    cold['LightGoldenrod3']=  '#CDBE70'
    cold['LightGoldenrod4']=  '#8B814C'
    cold['LightGoldenrodYellow']=  '#FAFAD2'
    cold['LightYellow']=  '#FFFFE0'
    cold['LightYellow1']=  '#FFFFE0'
    cold['LightYellow2']=  '#EEEED1'
    cold['LightYellow3']=  '#CDCDB4'
    cold['LightYellow4']=  '#8B8B7A'
    cold['PaleGoldenrod']=  '#EEE8AA'
    cold['PapayaWhip']=  '#FFEFD5'
    cold['cornsilk']=  '#FFF8DC'
    cold['cornsilk1']=  '#FFF8DC'
    cold['cornsilk2']=  '#EEE8CD'
    cold['cornsilk3']=  '#CDC8B1'
    cold['cornsilk4']=  '#8B8878'
    cold['gold']=  '#FFD700'
    cold['gold1']=  '#FFD700'
    cold['gold2']=  '#EEC900'
    cold['gold3']=  '#CDAD00'
    cold['gold4']=  '#8B7500'
    cold['goldenrod']=  '#DAA520'
    cold['goldenrod1']=  '#FFC125'
    cold['goldenrod2']=  '#EEB422'
    cold['goldenrod3']=  '#CD9B1D'
    cold['goldenrod4']=  '#8B6914'
    cold['moccasin']=  '#FFE4B5'
    cold['yellow']=  '#FFFF00'
    cold['yellow1']=  '#FFFF00'
    cold['yellow2']=  '#EEEE00'
    cold['yellow3']=  '#CDCD00'
    cold['yellow4']=  '#8B8B00'
    cold['copper']=  '#B87333'
    cold['gold']=  '#CD7F32'
    cold['silver']=  '#E6E8FA'
#    cold['red']=array([1,0,0])
#    cold['blue']=array([0,0,1])
#    cold['green']=array([0,1,0])
#    cold['black']=array([0,0,0])
#    cold['white']=array([0,0,0])
#    cold['maroon']=array([0.5,0,0])
#    cold['fuchsia']=array([1,0,1])
#    cold['purple']=array([0.5,0,0.5])
#    cold['lightblue']=array([0.67,0.84,0.9])
#    cold['cyan']=array([0,1,1])
#    cold['silver']=array([0.752,0.752,0.752])

    return cold

def createtrxfile(_filename,freq,phi,theta,Fpr,Fpi,Ftr,Fti):
    """
    Create antenna trx file    
    Usage:createtrxfile(filename,freq,phi,theta,Fpr,Fpi,Ftr,Fti)
    """
    filename=getlong(_filename,"ant")
    fo=open(filename,'w')

    for i in range(size(ravel(freq))):
        fo.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\n"%(ravel(freq)[i],ravel(phi)[i],ravel(theta)[i],ravel(Fpr)[i],ravel(Fpi)[i],ravel(Ftr)[i],ravel(Fti)[i]))
    fo.close()

def rgb(valex,out='int'):
    """
         convert a hexadecimal color into a (r,g,b) array
         >>> import pylayers.util.pyutil as pyu
         >>> coldic = pyu.coldict()
         >>> val = rgb(coldic['gold'],'float')
    """
    r = int(valex[1:3],16)
    g = int(valex[3:5],16)
    b = int(valex[5:7],16)

    col = np.array([r,g,b])
    if out == 'float':
        col = col/255.

    return(col)


def nbint(a):
    """
        calculate the number of intervals in a sequence of integer
    """
    b = a[1:]-a[0:-1]
    u = np.nonzero(b!=1)[0]
    return len(u)+1


def encodmtlb(lin):
    """
    encodmtlbi(lin) : encode python list of string in Matlab format

    Exemple : 

    lin = ['aaa','bbbbbbb','ccc','dd']
    F   = {}
    F['lin']=encodmtl(lin)
    io.savemat('encodmtlb_ex.mat',F)

    Principle: 

        The List is read column by column and written line by line in a same NxM matrix.
        If char does not exist it is replaced by space. 
    """
    
    # Determiner la taille (NxM) de la matrice
    N = len(lin)

    # Determiner la chaine de longueur maximale
    M  = 0
    lout = []
    str  = ''
    for  i in range(N):
        m = len(lin[i])
        if (m>M):
            M=m
        
    for  j in range(M):
        for  i in range(N):
            m = len(lin[i])
            k = j*N+i
            if (j>=m):
                c = ' '
            else:
                 c = lin[i][j]    
                
            str = str + c
            if mod(k+1,M)==0:
                lout.append(str)
                str=''
                     
    return(lout)

def sqrte(z):
    """ Evanescent SQRT for waves problems
     Sophocles J. Orfanidis - 1999-2008 - www.ece.rutgers.edu/~orfanidi/ewa
     Parameters
     ----------
     z : 
         array of complex numbers
     Returns
     -------
     y : 
         square root of z
       
         Notes: for z = a-j*b, y is defined as follows:
              [ sqrt(a-j*b),  if b~=0  
          y = [ sqrt(a),      if b==0 and a>=0              
              [ -j*sqrt(|a|), if b==0 and a<0   (i.e., the negative of what the ordinary SQRT gives)
         this definition is necessary to produce exponentially-decaying evanescent waves 
             (under the convention exp(j*omega*t) for harmonic time dependence)
        it is equivalent to the operation y = conj(sqrt(conj(a-j*b))), 
         but it fixes a bug in the ordinary SQRT in MATLAB arising whenever the real part is negative 
         and the imaginary part is an array with some zero elements. For example, compare the outputs:
         conj(sqrt(conj(-1 - array([0,1])*1j))) =      0 + 1.0000i,    sqrte(-1 - [0; 1]*j) =      0 - 1.0000i
                                          0.4551 - 1.0987i                            0.4551 - 1.0987i
        but

           conj(sqrt(conj(-1 + 0*1j))) = 0 - 1.000i,               sqrte(-1 + 0*j) = 0 - 1.000i
    """
    sh = np.shape(z)
    rz = np.ravel(z)
    y  = np.sqrt(rz)
    u  = np.nonzero((np.imag(rz)==0) & (np.real(rz)<0))[0]
    y[u] = -1j * np.sqrt(abs(rz[u]))
    y  = y.reshape(sh)
    return(y)

def untie(a,b):
    """
    
    Parameters
    ----------
    a
    b
    Returns
    -------
    boolean 
    a
    r 

    """
    la    = len(a)
    lb    = len(b)
    u     = np.intersect1d(a,b)
    lu    = len(u)
    #print lu
    #print min(la,lb)/2
    if lu >= min(la,lb)/2:
        # segment de a non commun avec b
        aa    = a[~np.in1d(a,u)]
        # segment de b non commun avec a
        bb    = b[~np.in1d(b,u)]
        r     = np.hstack((aa,bb))
        if la<lb:
            return(True,a,r)
        else:
            return(True,b,r)
    else:
        return(False,-1,-1)

def corrcy(a,b):
    """ cyclic matching correlation

    Parameters
    ----------
    a : array
    b : array

    Returns
    -------
    tk : 

    Example
    -------

    >>> a  = [1,2,3,4]
    >>> b  =  [1,2,3,4]
    >>> tk = corrcy(a,b)
    >>> assert tk[0]==4,'Problem in corrcy'

    """
    na = len(a)
    nb = len(b)
    tk = np.array([])
    if na>nb:
        for k in range(na):
            cak  = np.hstack((a[k::],a[0:k]))
            #print cak[0,nb]
            #print b
            diff = cak[0:nb]- b
            u    = np.nonzero(diff==0)[0]
            l    = len(u)
            tk   = np.hstack((tk,l))
    else:
        for k in range(nb):
            cbk  = np.hstack((b[k::],b[0:k]))
            #print cbk[0:na]
            #print a
            diff = cbk[0:na]- a
            u    = np.nonzero(diff==0)[0]
            l    = len(u)
            tk   = np.hstack((tk,l))

    return(tk)

def foo(var1, var2, long_var_name='hi') :
    r"""A one-line summary that does not use variable names or the
    function name.

    Several sentences providing an extended description. Refer to
    variables using back-ticks, e.g. `var`.

    Parameters
    ----------

    var1 : array_like
        Array_like means all those objects -- lists, nested lists, etc. --
        that can be converted to an array.  We can also refer to
        variables like `var1`.
    var2 : int
        The type above can either refer to an actual Python type
        (e.g. ``int``), or describe the type of the variable in more
        detail, e.g. ``(N,) ndarray`` or ``array_like``.
    Long_variable_name : {'hi', 'ho'}, optional
        Choices in brackets, default first when optional.

    Returns
    -------

    describe : type
        Explanation
    output : type
        Explanation
    tuple : type
        Explanation
    items : type
        even more explaining

    Other Parameters
    ----------------

    only_seldom_used_keywords : type
        Explanation
    common_parameters_listed_above : type
        Explanation

    Raises
    ------

    BadException
        Because you shouldn't have done that.

    See Also
    --------
    otherfunc : relationship (optional)
    newfunc : Relationship (optional), which could be fairly long, in which
              case the line wraps here.
    thirdfunc, fourthfunc, fifthfunc

    Notes
    -----

    Notes about the implementation algorithm (if needed).

    This can have multiple paragraphs.

    You may include some math:

    .. math:: X(e^{j\omega } ) = x(n)e^{ - j\omega n}

    And even use a greek symbol like :math:`\\omega` inline.

    .. plot:: 

        import matplotlib.pyplot as plt
        import numpy as np
        x = np.random.randn(1000)
        plt.hist( x, 20)
        plt.grid()
        plt.title(r'Normal: $\mu=%.2f, \sigma=%.2f$'%(x.mean(), x.std()))
        plt.show()

    References
    ----------

    Cite the relevant literature, e.g. [1]_.  You may also cite these
    references in the notes section above.

    .. [1] O. McNoleg, "The integration of GIS, remote sensing,
       expert systems and adaptive co-kriging for environmental habitat
       modelling of the Highland Haggis using object-oriented, fuzzy-logic
       and neural-network techniques," Computers & Geosciences, vol. 22,
       pp. 585-588, 1996.

    Examples
    --------

    These are written in doctest format, and should illustrate how to
    use the function.

    >>> a=[1,2,3]

    """

    pass
    

def cdf(x,color='b',label=" ",lw=1,xlabel="x",ylabel="CDF",logx=False):
    """ plot the cumulative density function of x

    Parameters
    ----------

    x :  np.array  (N)
    color : string 
        color symbol 
    label : string
        label
    lw: float  
        linewidth
    xlabel : string 
        xlabel 
    ylabel : string 
        ylabel 

    Examples
    --------

    .. plot::
        :include-source:

        >>> from matplotlib.pyplot import * 
        >>> import pylayers.util.pyutil as pyu
        >>> from scipy import *
        >>> import matplotlib.pylab as plt
        >>> x = randn(100)
        >>> pyu.cdf(x)
        >>> plt.show()

    """
    x  = np.sort(x)
    n  = len(x)
    x2 = np.repeat(x, 2)
    y2 = np.hstack([0.0, np.repeat(np.arange(1,n) / float(n), 2), 1.0])
    if logx:
        plt.semilogx(x2,y2,color=color,label=label,linewidth=lw)
    else:
        plt.plot(x2,y2,color=color,label=label,linewidth=lw)
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

def bitreverse(N=256,nbit=9):
    """ 
    Parameters
    ----------
    N : ideally a power of 2

    Returns
    -------
    t : list of the N integers in time reverse order

    Notes 
    -----
    This function is used for example in buildGv. 
    One error has been fixed  by forbidding the value 0 
    The value 0 is not return

    """
    t = []
    for k in np.arange(N-1)+1:
        b = BitString(uint=k,length=nbit) 
        b.reverse()
        b.ror(1)
        t.append(b.uint)
    return(np.array(t))


def timestamp(now):
    
    dt = dat.datetime.now()
    dn = str(dat.timedelta(seconds=float(now))).split(':')
    return (dt.strftime('%Y-%m-%d ')+dn[0]+':'+dn[1] +':'+dn[2][:2] +dn[2][2:5])

def writemeca(ID,time,p,v,a):
    """
    write mecanic information into text file:
        output/TruePosition.txt
        output/UWBSensorMeasurements.txt
    """


    ### TruePosition
    if not os.path.isfile(basename+'/' + pstruc['DIRNETSAVE'] +'/TruePosition.txt'):
        entete = 'TruePositionID,NodeID, Timestamp, X,Y,Z,ReferencePointID\n'
        file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/TruePosition.txt','w')
        file.write(entete)
        data = '1,'+str(ID) +','+ str(timestamp(time)) +',' + str(p[0])+',' +str(p[1])+','+',\n'
        file.write(data)
        file.close()
    else:
        file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/TruePosition.txt','r')
        lst=file.readlines()
        file.close()
        data = str(eval(lst[-1].split(',')[0])+1) +','+str(ID) +','+ str(timestamp(time)) +',' + str(p[0])+ ',' +str(p[1])+','+',\n'
        file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/TruePosition.txt','a')
        file.write(data)
        file.close()

    ### UWBSensorMeasurements
    if not os.path.isfile(basename+'/' + pstruc['DIRNETSAVE'] +'/UWBSensorMeasurements.txt'):
        entete = 'UWBSensorMeasurementsID,NodeID, Timestamp, UWB_MagX,UWB_MagY,UWB_MagZ,UWB_AccX,UWB_AccY,UWB_AccZ,UWB_GyroX,UWB_GyroY,UWB_GyroZ\n'
        file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/UWBSensorMeasurements.txt','w')
        file.write(entete)
        data = '1,'+str(ID) +','+ str(timestamp(time)) +',' + str(v[0])+',' +str(v[1])+',,'+str(a[0])+','+str(a[1])+',,,,\n'
        file.write(data)
        file.close()
    else:
        file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/UWBSensorMeasurements.txt','r')
        lst=file.readlines()
        file.close()
        data = str(eval(lst[-1].split(',')[0])+1)+',' +str(ID) +','+ str(timestamp(time)) +',' + str(v[0])+',' +str(v[1])+',,'+str(a[0])+','+str(a[1])+',,,,\n'
        file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/UWBSensorMeasurements.txt','a')
        file.write(data)
        file.close()



def writenet(net,t):
    """
    write network information into text file:
        netsave/ZIGLinkMeasurements.txt
        netsave/UWBLinkMeasurements.txt
    """
    for e in net.edges_iter(data=True):
        ### ZIGLinkMeasurements
        if not os.path.isfile(basename+'/' + pstruc['DIRNETSAVE'] +'/ZIGLinkMeasurements.txt'):
            entete = 'ZIGLinkMeasurementsID,NodeID, ZIG_PeerID, ZIG_RSSI, Timestamp\n'
            file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/ZIGLinkMeasurements.txt','w')
            file.write(entete)
            data = '1,'+ e[0] +','+ e[1] +',' + str(e[2]['Pr'][0]) +',' +timestamp(t.now()) +',\n'
            file.write(data)
            file.close()
        else:
            file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/ZIGLinkMeasurements.txt','r')
            lst=file.readlines()
            file.close()
            data = str(eval(lst[-1].split(',')[0])+1)+','+ e[0] +','+ e[1] +',' + str(e[2]['Pr'][0]) +',' +timestamp(t.now()) +',\n'
            file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/ZIGLinkMeasurements.txt','a')
            file.write(data)
            file.close()

        ### UWBLinkMeasurements
        if not os.path.isfile(basename+'/' + pstruc['DIRNETSAVE'] +'/UWBLinkMeasurements.txt'):
            entete = 'UWBLinkMeasurementsID, NodeID, Timestamp, UWB_PeerID, UWB_Dist, UWB_BER, UWB_FER, UWB_CIR\n'
            file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/UWBLinkMeasurements.txt','w')
            file.write(entete)
            data = '1,'+ e[0] +','+ timestamp(t.now()) +',' +e[1] +','+ str(e[2]['d']) +',,,,\n'
            file.write(data)
            file.close()
        else:
            file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/UWBLinkMeasurements.txt','r')
            lst=file.readlines()
            file.close()
            data = str(eval(lst[-1].split(',')[0])+1)+','+ e[0] +','+ timestamp(t.now()) +',' +e[1] +','+ str(e[2]['d']) +',,,,\n'
            file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/UWBLinkMeasurements.txt','a')
            file.write(data)
            file.close()


#        self.insertitem1("ACOLinkMeasurements",('NodeID',
#                                            'ACO_PeerID',
#                                            'ACO_RSSI',
#                                            'Timestamp'),
#                                            (eval(e[0]),
#                                            eval(e[1]),
#                                            e[2]['Pr'][0],
#                                            pyu.timestamp(t)))
#        self.insertitem1("CEALinkMeasurements",('NodeID',
#                                            'Timestamp',
#                                            'CEA_PeerID',
#                                            'CEA_Dist'),
#                                            (eval(e[0]),
#                                            pyu.timestamp(t),
#                                            eval(e[1]),
#                                            e[2]['d']))


def writenode(agent):
    '''
    write Nodes.txt
    '''
    if not os.path.isfile(basename+'/' + pstruc['DIRNETSAVE'] +'/Nodes.txt'):
        entete = 'NodeID, NodeName, NodeOwner, NodeDescription, NodeOwnerID, Mobile OrAnchor, TrolleyID\n'
        file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/Nodes.txt','w')
        file.write(entete)
        file.close()

    data = str(eval(agent.ID)) +','+ agent.name + ',,,,' + str(agent.MoA) +',\n'
    file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/Nodes.txt','a')
    file.write(data)
    file.close()

def writeDetails(t,description='simulation', location ='Rennes'):
    '''
    write MeasurementsDetails.txt
    '''
    if not os.path.isfile(basename+'/' + pstruc['DIRNETSAVE'] +'/MeasurementsDetails.txt'):
        entete = 'MeasurementsDetailsID, MeasurementsDate, MeasurementsDescription, MeasurementsLocation\n'
        file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/MeasurementsDetails.txt','w')
        file.write(entete)
        file.close()

    data = '1' +','+ timestamp(t.now()) + ', ' +description + location +',\n'
    file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/MeasurementsDetails.txt','a')
    file.write(data)
    file.close()

def zipd(path, zipfilename):
    """
    add a folder to a zipfile
    Parameters
    ----------

        filepath : string
        zipfilename : string
    """
    zip_file = zipfile.ZipFile(zipfilename, 'a')
    for (dirpath, dirnames, filenames) in os.walk(path):
        for dirc in dirnames:
            zip_file.write(os.path.join(dirpath, dirc),
            os.path.join(os.path.basename(path),
            os.path.join(dirpath, dirc)[len(path):]))
        for fil in filenames:
            zip_file.write(os.path.join(dirpath, fil),
            os.path.join(os.path.basename(path),
            os.path.join(dirpath, fil)[len(path):]))
    zip_file.close()

def unzipd(path, zipfilename):
    """
    unzip a zipfile to a folder
    Parameters
    ----------

        filepath : string
        zipfilename : string
    """
    zip_file = zipfile.ZipFile(zipfilename)
    if not os.path.isdir(path):
        os.makedirs(path)    

    for each in zip_file.namelist():
        print each
        if not each.endswith('/'): 
            root, name = os.path.split(each)
            directory = os.path.normpath(os.path.join(path, root))
            if not os.path.isdir(directory):
                os.makedirs(directory)
            file(os.path.join(directory, name),
                'wb').write(zip_file.read(each))

def unzipf(path, filepath, zipfilename):
    """
    unzip a file from zipfile to a folder
    Parameters
    ----------

        filepath : string
        zipfilename : string
    """
    zip_file = zipfile.ZipFile(zipfilename)
    if not os.path.isdir(path):
        os.makedirs(path)    

    for each in zip_file.namelist():
        if each == filepath and not each.endswith('/'): 
            root, name = os.path.split(each)
            directory = os.path.normpath(os.path.join(path, root))
            if not os.path.isdir(directory):
                os.makedirs(directory)
            file(os.path.join(directory, name),
                'wb').write(zip_file.read(each))

def rotate_line(A,B,theta):
    """
    rotation of a line [AB] of an angle theta with A fixed
    Parameters
    ----------

        A: ndarray
        B: ndarray
        theta: float

    Returns
    -------

        Br: ndarry
    """
    if np.shape(B)!=(2,1):
        B.reshape((2,1))
    R = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    Br=np.dot(R,B)
    return Br

def extract_block_diag(A,M,k=0):
    """Extracts blocks of size M from the kth diagonal
    of square matrix A, whose size must be a multiple of M."""

    # Check that the matrix can be block divided
    if A.shape[0] != A.shape[1] or A.shape[0] % M != 0:
        raise StandardError('Matrix must be square and a multiple of block size')

    # Assign indices for offset from main diagonal
    if abs(k) > M - 1:
        raise StandardError('kth diagonal does not exist in matrix')
    elif k > 0:
        ro = 0
        co = abs(k)*M 
    elif k < 0:
        ro = abs(k)*M
        co = 0
    else:
        ro = 0
        co = 0
    blocks = np.array([A[i+ro:i+ro+M,i+co:i+co+M] 
                       for i in range(0,len(A)-abs(k)*M,M)])
    return blocks

def fill_block_diag(A, blocks,M,k=0):
    """fill A with blocks of size M from the kth diagonal
    """

    # Check that the matrix can be block divided
    if A.shape[0] != A.shape[1] or A.shape[0] % M != 0:
        raise StandardError('Matrix must be square and a multiple of block size')

    # Assign indices for offset from main diagonal
    if abs(k) > M - 1:
        raise StandardError('kth diagonal does not exist in matrix')
    elif k > 0:
        ro = 0
        co = abs(k)*M 
    elif k < 0:
        ro = abs(k)*M
        co = 0
    else:
        ro = 0
        co = 0
    for i in range(0,len(A)-abs(k)*M,M):
        A[i+ro:i+ro+M,i+co:i+co+M]=blocks[int(i/M),:,:] 
    return A

if __name__ == "__main__":
    doctest.testmod()
