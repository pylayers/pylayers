# -*- coding:Utf-8 -*-

# from pylayers.antprop.loss import *
# NPT=10000
# x=np.array([0,0,8])
# x=x.reshape(3,1)
# y = np.ones((3,NPT))
# y[0,:]=0
# y[1,:]=np.arange(NPT)
# y[2,:]=2
# g0=1
# g1=1
# fGHz=2.4
# r=two_rays_simpleloss(x,y,g0,g1,mode='PL',fGHz)
# pp = PL(fGHz,x,y,2)
# plt.semilogx(10*np.log10(r),label='two-ray model')
# plt.semilogx(-pp[0,:],label='one slope model')
# plt.axis([10,NPT,-150,-50])
# plt.show()


def doe(r,h0,h1,k=4/3.) :
    """

    """


    r0 = 6371e3 # earth radius
    re = k*r0 # telecom earth radius

    # r = distance curviligne entre TXetRX / geodesic

    p = 2/(np.sqrt(3))*np.sqrt(re*(h0+h1)+(r**2/4.)) # eq 8.45
    eps = np.arcsin(2*re*r*(h1-h0)/p**3) # eq 8.46



    # distance of reflection on curved earth
    r1 = r/2 - p*np.sin(eps/3) # eq 8.44

    r2 = r -r1

    phi1 = r1/re # 8.47
    phi2 = r2/re # 8.48

    R1 = np.sqrt(h0**2+4*re*(re+h0)*(np.sin(phi1/2))**2) # 8.51
    R2 = np.sqrt(h1**2+4*re*(re+h1)*(np.sin(phi2/2))**2) #8.52

    Rd = np.sqrt((h1-h0)**2+4*(re+h1)*(re+h0)*np.sin((phi1+phi2)/2.)**2) # 8.53

    # tangente angle on earth
    psy = np.arcsin((h1/R1)-R1/(2*re)) # eq 8.55
    deltaR = 4*R1*R2*np.sin(psy)**2/(R1+R2+Rd)

    dloss = Rd
    dref = R1+R2

    return dloss,dref






from pylayers.antprop.loss import *
from geopy.geocoders import Nominatim
import pylayers.gis.gisutil as gu


# geolocator = Nominatim()
# info,latlon=geolocator.geocode('ietr rennes')
# lat0=latlon[0]
# lon0=latlon[1]

# lat0 = 48.1168958
# lon0 = -1.64065078258313

# h0 = 20
# h1 = 20
# fGHz=2.4
# l = []
# r=[]
# d=[]
# for k in np.arange(0,100 ,1e-2):
#     r.append(two_ray_curvedearth(lat0,lon0,lat0,lon0+k,h0,h1))
#     d.append(gu.distance_on_earth(lat0,lon0,lat0,lon0+k))
# plt.semilogx(d,r,label='two-ray model')





# def lossref_compute(P,h0,h1,k=4/3.) :
#     """
#     compute loss and reflection rays on curved earth

#     Parameters
#     ----------

#     P : float | list 

#         if len(P) == 1 => P is a distance
#         if len(P) == 4 => P is a list of [lon0,lat0,lon1,lat1]

#         where :
#         lat0 : float | string
#             latitude first point (decimal | deg min sec Direction)
#         lat1 : float | string
#             latitude second point (decimal | deg min sec Direction)
#         lon0 : float | string
#             longitude first point (decimal | deg min sec Direction)
#         lon1 : float | string
#             longitude second point (decimal | deg min sec Direction)
#     h0 : float:
#         height of 1st point 
#     h1 : float:
#         height of 2nd point 
#     k : electromagnetic earth factor


#     Returns
#     -------
    
#     dloss : 
#         length of direct path (meter)
#     dref :
#         length of reflective path (meter)

#     References
#     ----------
#     B. R. Mahafza, Radar systems analysis and design using MATLAB, Third edition. Boca Raton ; London: CRC/Taylor & Francis, chapter 8, 2013.

#     """


#     if isinstance(P,float) or isinstance(P,int) :
#         # P is a distance
#         r=P
#         mode = 'dist'
#     elif isinstance(P,np.ndarray) or isinstance(P,list):
#         if len(P) == 1:
#             # P is a distance
#             r=P
#             mode = 'dist'
#         elif len(P) == 4:
#             # P is a lonlat
#             lat0=P[0]
#             lon0=P[1]
#             lat1=P[2]
#             lon1=P[2]
#             mode = 'lonlat'
#         else :
#             raise AttributeError('P must be a list [lat0,lon0,lat1,lon0] or a distance')
#     else :
#         raise AttributeError('Invalid P format ( list | ndarray )')

#     r0 = 6371e3 # earth radius
#     re = k*r0 # telecom earth radius


#     if mode == 'lonlat':
#         # r = distance curviligne entre TXetRX / geodesic
#         r = gu.distance_on_earth(lat0, lon0, lat1, lon1)
#     else :
#         r=P


#     p = 2/(np.sqrt(3))*np.sqrt(re*(h0+h1)+(r**2/4.)) # eq 8.45
#     eps = np.arcsin(2*re*r*(h1-h0)/p**3) # eq 8.46



#     # distance of reflection on curved earth
#     r1 = r/2 - p*np.sin(eps/3) # eq 8.44

#     r2 = r -r1

#     phi1 = r1/re # 8.47
#     phi2 = r2/re # 8.48

#     R1 = np.sqrt(h0**2+4*re*(re+h0)*(np.sin(phi1/2))**2) # 8.51
#     R2 = np.sqrt(h1**2+4*re*(re+h1)*(np.sin(phi2/2))**2) #8.52

#     Rd = np.sqrt((h1-h0)**2+4*(re+h1)*(re+h0)*np.sin((phi1+phi2)/2.)**2) # 8.53

#     # tangente angle on earth
#     psy = np.arcsin((h1/R1)-R1/(2*re)) # eq 8.55
#     deltaR = 4*R1*R2*np.sin(psy)**2/(R1+R2+Rd)

#     dloss = Rd
#     dref = R1+R2

#     return dloss,dref

# def two_ray_curvedearth(P,h0,h1,fGHz=2.4,**kwargs):
#     """


#     Parameters
#     ----------

#     P : float | list 

#         if len(P) == 1 => P is a distance
#         if len(P) == 4 => P is a list of [lon0,lat0,lon1,lat1]

#         where :
#         lat0 : float | string
#             latitude first point (decimal | deg min sec Direction)
#         lat1 : float | string
#             latitude second point (decimal | deg min sec Direction)
#         lon0 : float | string
#             longitude first point (decimal | deg min sec Direction)
#         lon1 : float | string
#             longitude second point (decimal | deg min sec Direction)
#     h0 : float:
#         height of 1st point 
#     h1 : float:
#         height of 2nd point 
#     fGHz : float
#         frequency (GHz)


#     k : float
#         electromagnetic earth factor
#     Gt : float
#         Transmitter Antenna Gain (dB)
#     Gr : float
#         Receiver Antenna Gain (dB)
#     gamma : 
#         Reflexion coeff(-1)

#     mode : PL | E (default : PL)
#         return Energy (E) or Path loss/power loss (PL)
#     dB : boolean (True)
#         return result in dB


#     Example
#     -------

#     .. plot::
#         :include-source:
        
#         >>> from pylayers.antprop.loss import *
#         >>> import pylayers.gis.gisutil as gu
#         >>> lat0 = 48.1168958
#         >>> lon0 = -1.64065078258313
#         >>> h0 = 20
#         >>> h1 = 20
#         >>> fGHz=2.4
#         >>> r=[]
#         >>> d=[]
#         >>> for k in np.arange(0,100 ,1e-2):
#         >>>     r.append(two_ray_curvedearth(lat0,lon0,lat0,lon0+k,h0,h1))
#         >>>     d.append(gu.distance_on_earth(lat0,lon0,lat0,lon0+k))
#         >>> plt.semilogx(d,r,label='two-ray model')

#     """


#     defaults = { 'Gt':0.,
#                  'Gr':0.,
#                  'k':4/3.,
#                  'gamma': -1.,
#                  'mode':'PL',
#                  'dB':True
#                }

#     for k in defaults:
#         if k not in kwargs:
#             kwargs[k]=defaults[k]

#     Gt=kwargs.pop('Gt')
#     Gr=kwargs.pop('Gr')

#     Gt = 10**((1.*Gt)/10.)
#     Gr = 10**((1.*Gr)/10.)
#     k=kwargs.pop('k')
#     gamma=kwargs.pop('gamma')
    
#     dloss,dref = lossref_compute(P,h0,h1,k)


#     lossterm = np.exp(2.j*np.pi*dloss*fGHz/0.3) / (1.*dloss)
#     refterm = np.exp(2.j*np.pi*dref*fGHz/0.3) / (1.*dref)


    

#     E = np.sqrt(Gt*Gr)*0.3/(4*np.pi*fGHz)* (lossterm + gamma*refterm)

#     if kwargs['mode'] == 'E':
#         return E
#     if kwargs['mode'] == 'PL':
#         if kwargs['dB'] :
#             return 20*np.log10(E)
#         else:
#             return E**2



# fGHz=2.4
# p0=np.array(([0,0,20]))
# p1=np.array(([0,1,20]))
# p0=p0.reshape(3,1)
# p1=p1.reshape(3,1)
# OR = []
# TRF = []
# TRC = []


# for d in np.arange(1,10000,0.5):
#     p1[1,:]=d
#     OR.append(-PL(fGHz,p0[:2,:],p1[:2,:],2)[0])
#     TRF.append(two_rays_flatearth(p0[:,0],p1[:,0],0.,0.,fGHz))
#     TRC.append(two_ray_curvedearth(d,p0[2,:],p1[2,:],fGHz))
#     # l.append(-( pl0+ 20*np.log10(k)))
# plt.semilogx(TRF,label='two-ray model flat earth')
# plt.semilogx(TRC,label='two-ray model curved earth')
# plt.semilogx(OR,label='one-ray model')
# plt.legend()
# plt.show()





from pylayers.antprop.loss import *
import matplotlib.pyplot as plt
fGHz=3.5
p0=np.array(([0,1,52]))

p1=np.array(([0,0,12]))
p0=p0.reshape(3,1)
p1=p1.reshape(3,1)
OR = [] # One Ray model
TRF = [] # Two Ray model on flat earth
TRC = [] # Two Ray model on curved earth
Gt=12
Gr=12
dr= np.arange(1,1000,1.0)
for d in dr:
    p1[1,:]=d
    OR.append(-PL(fGHz,p0[:2,:],p1[:2,:],2)[0]+Gt+Gr)
    TRF.append(two_rays_flatearth(p0[:,0],p1[:,0],Gt=Gt,Gr=Gr,fGHz=fGHz))
    TRC.append(two_ray_curvedearth(d,p0[2,:],p1[2,:],Gt=Gt,Gr=Gr,fGHz=fGHz))
plt.semilogx(TRF,label='two-ray model flat earth')
plt.semilogx(dr,20*np.log10(1./dr)+20,label='1/d')
plt.semilogx(dr,20*np.log10(1./dr**2)+20,label='1/d**2')

plt.semilogx(dr,TRC,label='two-ray model curved earth')
plt.semilogx(dr,OR,label='one-ray model')
plt.legend()
plt.show()