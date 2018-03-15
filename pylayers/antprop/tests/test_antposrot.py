from pylayers.antprop.antenna import *
from pylayers.util.geomutil import *
import mayavi.mlab as mlab

deg_to_rad = np.pi/180.
R = 0.660
H = 1.5
theta = np.array([30,60,90])*deg_to_rad
phi = np.array([45,135,225,315])*deg_to_rad
x = R*np.sin(theta[:,None])*np.cos(phi[None,:])
y = R*np.sin(theta[:,None])*np.sin(phi[None,:])
z = H + R*np.cos(theta[:,None])*np.ones(len(phi))[None,:] 
p = np.vstack((x.ravel(),y.ravel(),z.ravel()))
th = (theta[:,None]*np.ones(len(phi))[None,:]).ravel()
ph = (np.ones(len(theta))[:,None]*phi[None,:]).ravel()
po = np.array([0,0,H])
tA =  []
z = np.array([0,0,1])
tx = np.arange(-0.3,0.3,0.005)
ty = np.arange(-0.3,0.3,0.005)
try:
    del pg 
except:
    pass
for x in tx:
    for y in ty:
        pt = np.array([x,y,H])[None,:]
        try:
            pg = np.vstack((pg,pt))
        except:
            pg = pt

try:
    del E
except:
    pass
#pg=np.array([[0,0,1.5]])
for k in range(12):
    #T  = MEulerAngle(ph[k],th[k],0)
    v = po-p[:,k]
    rn = v/np.sqrt(np.sum(v*v,axis=0))
    phi = np.cross(z,rn)
    phin = phi/np.sqrt(np.sum(phi*phi,axis=0))
    thetan = np.cross(phin,rn)
    T = np.hstack((phin[:,None],rn[:,None],thetan[:,None]))
    n = np.random.rand(1)
    #if n>0.5: 
    #    A = AntPosRot('CEA_OTA_p2.vsh3',p=p[:,k],T = T)
    #else:
    A = AntPosRot('CEA_OTA.vsh3',p=p[:,k],T = T)
    tA.append(A)
#pg=np.array([[0,0,1.5],[0.1,0.1,1.5]])
    try:
        E = E + A.field(pg)
    except:
        E = A.field(pg)

U = E.reshape(len(tx),len(ty),3,9)
mU = np.abs(U)
vmax = np.max(mU)
for k in range(9):
    plt.figure(figsize=(20,10))
    plt.subplot(231)
    plt.imshow(np.abs(U[:,:,0,k]),extent=[-0.3,0.3,-0.3,0.3],vmax=vmax,cmap='jet')
    plt.colorbar()
    plt.title('|Ex| : f = '+str(A.fGHz[k])+' GHz',fontsize=18)
    plt.subplot(232)
    plt.imshow(np.abs(U[:,:,1,k]),extent=[-0.3,0.3,-0.3,0.3],vmax=vmax,cmap='jet')
    plt.colorbar()
    plt.title('|Ey| : f = '+str(A.fGHz[k])+' GHz',fontsize=18)
    plt.subplot(233)
    plt.imshow(np.abs(U[:,:,2,k]),extent=[-0.3,0.3,-0.3,0.3],vmax=vmax,cmap='jet')
    plt.colorbar()
    plt.title('|Ez| : f = '+str(A.fGHz[k])+ ' GHz',fontsize=18)
    plt.subplot(234)
    plt.imshow(np.angle(U[:,:,0,k]),extent=[-0.3,0.3,-0.3,0.3],vmax=np.pi,vmin=-np.pi,cmap='jet')
    plt.colorbar()
    plt.title('arg(Ex) : f = '+str(A.fGHz[k])+' GHz',fontsize=18)
    plt.subplot(235)
    plt.imshow(np.angle(U[:,:,1,k]),extent=[-0.3,0.3,-0.3,0.3],vmax=np.pi,vmin=-np.pi,cmap='jet')
    plt.colorbar()
    plt.title('arg(Ey) : f = '+str(A.fGHz[k])+' GHz',fontsize=18)
    plt.subplot(236)
    plt.imshow(np.angle(U[:,:,2,k]),extent=[-0.3,0.3,-0.3,0.3],vmax=np.pi,vmin=-np.pi,cmap='jet')
    plt.colorbar()
    plt.title('arg(Ez) : f = '+str(A.fGHz[k])+' GHz',fontsize=18)
    plt.savefig('figOTA_pz_'+str(k)+'.png')
#vmax= = np.max(abs(E))
#plt.subplot(321)
#plt.imshow(np.abs(E[:,:,Ea
#A._show3(scale=0.2)
#A._show3()
#B = AntPosRot('cst',p = np.array([-1,-2,3]),T = np.eye(3))
#B._show3()
