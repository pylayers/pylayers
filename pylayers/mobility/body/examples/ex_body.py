from pylayers.mobility.body.body import *


# select a number of frame 
nframes = 126

# create a body 
# each agent should have a body which may be 
# extracted from the same motion capture file 
# for different agent
# A BodyCylinder object is a member of a mobile agent
#

#
# create a bodyCylinder object
#
bc = BodyCylinder()

#
# load a .c3dmotion capture file
# this update the g.pos 
#
bc.loadC3D(filename='07_01.c3d',nframes=nframes)
#
# plot the graph representation of the object
# 
plt.figure()
plt.title("Body graph 15 Nodes - 10 Edges ")
nx.draw(bc.g,bc.g.pos)

plt.axis('scaled')

############# verification de la base
#B.LoadMotion()
#c3d_frames = B.d
#B.CylinderModel()
#fr_id  = 125
#B.cylinder_basis0()
#B.cylinder_basis_k (frameId = fr_id)
#P = np.array([1,1,1])
#for cy_id in range(10):
#A = np.zeros(shape = (3))
#u0 = B.basis0[cy_id,0:3]
#v0 = B.basis0[cy_id,3:6]
#w0 = B.basis0[cy_id,6:]
#fig  = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#SF = 1
#plt.plot((A[0],SF*u0[0]+A[0]),(A[1],SF*u0[1]+A[1]),(A[2],SF*u0[2]+A[2]), 'y')
#plt.plot((A[0],SF*v0[0]+A[0]),(A[1],SF*v0[1]+A[1]),(A[2],SF*v0[2]+A[2]),'r')
#plt.plot((A[0],SF*w0[0]+A[0]),(A[1],SF*w0[1]+A[1]),(A[2],SF*w0[2]+A[2]),'g')
#uk = B.basisk[cy_id,0:3]
#vk = B.basisk[cy_id,3:6]
#wk = B.basisk[cy_id,6:]
#plt.plot((A[0],SF*uk[0]+A[0]),(A[1],SF*uk[1]+A[1]),(A[2],SF*uk[2]+A[2]), 'y--')
#plt.plot((A[0],SF*vk[0]+A[0]),(A[1],SF*vk[1]+A[1]),(A[2],SF*vk[2]+A[2]),'r--')
#plt.plot((A[0],SF*wk[0]+A[0]),(A[1],SF*wk[1]+A[1]),(A[2],SF*wk[2]+A[2]),'g--')
#plt.autoscale()
#plt.title(str(cy_id))
#plt.show()
############   Antenna on cylinders


# Load a .c3dmotion capture file
bc.loadC3D(filename='07_01.c3d',nframes=nframes)

#
# extract a numpy array (3 x np x nf )
#
#   3  : x y  z
#   np : number of points
#   nf : number of frames
#

c3dframe = bc.d
bc.geomfile(0)
bc.movie()
pg = np.sum(bc.d,axis=1)

#
# Convert c3d file into a 10 (4+4+1+1) cylinder model
#
#  4 arms 4 legs 1 trunk 1 head
#
#  create member c array
#  c : array(shape  =  (nc,8)), Id (1), A coordinate (3) , B coordinate (3) , radius (1)
#

#bc.antennas()

#
# create a reference basis on each cylinder
#
# WARNING : arbitray orientation for the first frame
#

#bc.cylinder_basis0()
#
##
## Below an example of positioning antenna around cylnders
##
##
#
#fr_ind = 0
#cyl = B.c[:, :, fr_ind]
## scale factor
#SF = 1
#
## 10 antennas on 10 cylinders
#
#cy_id = range(0, 10)
#l = np.ones((10))
#alpha = np.pi / 2 * np.ones((10))
#
#
##
## Display Frame 0 
##
#B.cyl_antenna(cy_id, l, alpha, 0)
#
#fig = plt.figure()
#
#ax = fig.add_subplot(111, aspect='equal', projection='3d')
#
#for i in range(B.c.shape[0]):
#    pltu.cylinder(fig, cyl[i, 1:4], cyl[i, 4:7], cyl[i, 7])
#    A3 = cyl[i, 1:4]
#    u3 = B.basis0[i, 0:3]
#    v3 = B.basis0[i, 3:6]
#    w3 = B.basis0[i, 6:]
#
#    plt.plot((A3[0], SF * u3[0] + A3[0]), (
#        A3[1], SF * u3[1] + A3[1]), (A3[2], SF * u3[2] + A3[2]), 'y--')
#    plt.plot((A3[0], SF * v3[0] + A3[0]), (
#        A3[1], SF * v3[1] + A3[1]), (A3[2], SF * v3[2] + A3[2]), 'r--')
#    plt.plot((A3[0], SF * w3[0] + A3[0]), (
#        A3[1], SF * w3[1] + A3[1]), (A3[2], SF * w3[2] + A3[2]), 'g--')
#    plt.plot((A3[0], B.ant[i, 0] + A3[0]), (
#        A3[1], B.ant[i, 1] + A3[1]), (A3[2], B.ant[i, 2] + A3[2]), 'c*--')
#    #~ plt.axis('tight')
#    #~ fig  = plt.figure()
#    #~ plt.plot((A3[0],SF*u3[0]+A3[0]),(A3[1],SF*u3[1]+A3[1]), 'y--')
#    #~ plt.plot((A3[0],SF*v3[0]+A3[0]),(A3[1],SF*v3[1]+A3[1]),'r--')
#    #~ plt.plot((A3[0],SF*w3[0]+A3[0]),(A3[1],SF*w3[1]+A3[1]),'g--')
#        #~
#    #~ plt.plot((A3[0],B.ant[0]+ A3[0]) ,(A3[1],B.ant[1]+ A3[1]),'c*--')
#    #~ plt.axis('scaled')
#
##
## Display Frame 60 
##
#fr_ind2 = 60
#cyl1 = B.c[:, :, fr_ind2]
#B.cylinder_basis_k(frameId=fr_ind2)
#B.cyl_antenna(cy_id, l, alpha, fr_ind2)
#SF = 1
#
#fig = plt.figure()
#for i in range(B.c.shape[0]):
#    pltu.cylinder(fig, cyl1[i, 1:4], cyl1[i, 4:7], cyl1[i, 7])
#    A3 = cyl1[i, 1:4]
#    u3 = B.basisk[i, 0:3]
#    v3 = B.basisk[i, 3:6]
#    w3 = B.basisk[i, 6:]
#
#    plt.plot((A3[0], SF * u3[0] + A3[0]), (
#        A3[1], SF * u3[1] + A3[1]), (A3[2], SF * u3[2] + A3[2]), 'y--')
#    plt.plot((A3[0], SF * v3[0] + A3[0]), (
#        A3[1], SF * v3[1] + A3[1]), (A3[2], SF * v3[2] + A3[2]), 'r--')
#    plt.plot((A3[0], SF * w3[0] + A3[0]), (
#        A3[1], SF * w3[1] + A3[1]), (A3[2], SF * w3[2] + A3[2]), 'g--')
#    plt.plot((A3[0], B.ant[i, 0] + A3[0]), (
#        A3[1], B.ant[i, 1] + A3[1]), (A3[2], B.ant[i, 2] + A3[2]), 'c*--')
#    #~ plt.axis('tight')
#
#plt.show()
#
#####################
#B.LoadMotion()
#f = B.d
#cycle_distance = dist(f[:, 0, 0], f[:, 0, -1])
#print 'distance parcouru pendant un cycle de marche = ', cycle_distance
##~ for fr_id in range(0,124):
#    #~ inter_frames_dist = dist(f[:,0,fr_id],f[:,0,fr_id  +1])
#    #~ print ' frame ID = ' , fr_id ,  'distance parcouru pendant un step = ', inter_frames_dist
##~
##~
#
######## test translation and rotation
#fig = plt.figure()
#
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(f[0, :, 0], f[1, :, 0], f[2, :, 0], c='g', marker='o')
#ax.scatter(f[0, :, -1], f[1, :, -1], f[2, :, -1], c='g', marker='*')
#
## translation
##~ p = np.array([0,0,0])
##~ f1 = translate(f,p)
##~ ax.scatter(f1[0,:,0],f1[1,:,0],f1[2,:,0] )
##~ ax.scatter(f1[0,:,-1],f1[1,:,-1],f1[2,:,-1] )
## rotantion
##~ fig = plt.figure()
##~
##~ ax = fig.add_subplot(111, projection='3d')
#f2 = rotation(f, np.pi / 2)
#ax.scatter(f2[0, :, 0], f2[1, :, 0], f2[2, :, 0], c='r', marker='o')
#ax.scatter(f2[0, :, -1], f2[1, :, -1], f2[2, :, -1], c='r', marker='*')
#plt.show()
