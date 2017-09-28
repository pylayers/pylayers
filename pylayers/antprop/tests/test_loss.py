from pylayers.antprop.loss import *
import matplotlib.cm as cm
import pdb

<<<<<<< HEAD
fGHz = np.array([24,25,26])
# fGHz = fGHz=np.arange(18,28,0.1)
=======
fGHz = np.array([3,25,60])
>>>>>>> b27014e9a940b367ff58465846bdb402cd4a9d8e
h = 0.75
w = 0.3
tx = np.array([0,0,10])    # over axis z
rx = np.array([0,0,-10])   # over axis z

pg = np.array([0,0,0])
uw = np.array([1,0,0])     # over axis x
uh = np.array([0,1,0])     # over axis y

L = LossMetisShadowing(fGHz,tx,rx,pg,uw,uh,w,h)

#######################################################

# fGHz = np.array([24,25,26])
h = np.array([0.75])
w = np.array([0.3])
tx = np.array([0,0,10]).reshape(3,1)     # over axis z
rx = np.array([0,0,-10]).reshape(3,1)    # over axis z
# rx = np.array([0,0,-1]).reshape(3,1)   # near the screen

pg = np.array([0,0,0]).reshape(3,1)
uw = np.array([1,0,0]).reshape(3,1)     # over axis x.reshape(3,1)
uh = np.array([0,1,0]).reshape(3,1)     # over axis y
L2 = LossMetisShadowing2(fGHz,tx,rx,pg,uw,uh,w,h)

#######################################################


# Nseg = 10
# Nscreen = 5 
# tx = 10*np.random.rand(3,Nseg)     # over axis z
# rx = -10*np.random.rand(3,Nseg)
# pg = 5*np.random.rand(3,Nscreen)
# uw = np.random.rand(3,Nscreen)     # over axis x
# uh = np.random.rand(3,Nscreen)     # over axis y

# # normalization of the unitary vector

# muw = np.sqrt(np.sum(uw*uw,axis=0))
# muh = np.sqrt(np.sum(uh*uh,axis=0))
# uwn = uw/muw
# uhn = uh/muh
# pwh = np.sum(uwn*uhn,axis=0)
# uhnn = uhn-pwh*uwn
# muhnn = np.sqrt(np.sum(uhnn*uhnn,axis=0))
# uhnnn = uhnn/muhnn
# w = 2*np.random.rand(Nscreen)
# h = 2*w
# L2 = LossMetisShadowing2(fGHz,tx,rx,pg,uwn,uhnnn,w,h)

# print "L   :",L
# print "L2  :",L2
# plt.imshow(L2[:,:,1],cmap=cm.jet,lower='origin')
# plt.colorbar()
# plt.axis('auto')


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# # ax.scatter(tx,rx,Nseg,marker='o', cmap=cm)
# use_colours = {"L": "red", "A": "green", "B": "blue"}
# ax.scatter(tx,rx,Nscreen,c=[use_colours[x[0]] for x in d],marker='o', cmap='hot')

# ax.set_xlabel('Tx')
# ax.set_ylabel('Rx')
# ax.set_zlabel('Nseg')
# plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=8, bbox_to_anchor=(0, 0))