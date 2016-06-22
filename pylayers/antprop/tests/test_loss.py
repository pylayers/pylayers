from pylayers.antprop.loss import *
import pdb

fGHz = np.array([24,25,26])
h = 0.75
w = 0.3
tx = np.array([0,0,10])    # over axis z
rx = np.array([0,0,-10])   # over axis z
pg = np.array([0,0,0])
uw = np.array([1,0,0])     # over axis x
uh = np.array([0,1,0])     # over axis y

L = LossMetisShadowing(fGHz,tx,rx,pg,uw,uh,w,h)

fGHz = np.array([24,25,26])
h = np.array([0.75])
w = np.array([0.3])
tx = np.array([0,0,10]).reshape(3,1)    # over axis z
rx = np.array([0,0,-1]).reshape(3,1)   # over axis z
pg = np.array([0,0,0]).reshape(3,1)
uw = np.array([1,0,0]).reshape(3,1)     # over axis x.reshape(3,1)
uh = np.array([0,1,0]).reshape(3,1)     # over axis y
#Nseg = 10
#Nscreen = 5 
#tx = 10*np.random.rand(3,Nseg)    # over axis z
#rx = -10*np.random.rand(3,Nseg)
#pg = 5*np.random.rand(3,Nscreen)
#uw = np.random.rand(3,Nscreen)     # over axis x
#uh = np.random.rand(3,Nscreen)     # over axis y
#muw = np.sqrt(np.sum(uw*uw,axis=0))
#muh = np.sqrt(np.sum(uh*uh,axis=0))
#uwn = uw/muw
#uhn = uh/muh
#pwh = np.sum(uwn*uhn,axis=0)
#uhnn = uhn-pwh*uwn
#muhnn = np.sqrt(np.sum(uhnn*uhnn,axis=0))
#uhnnn = uhnn/muhnn
#w = 2*np.random.rand(Nscreen)
#h = 2*w
#L2 = LossMetisShadowing2(fGHz,tx,rx,pg,uwn,uhnnn,w,h)
L2 = LossMetisShadowing2(fGHz,tx,rx,pg,uw,uh,w,h)
