from numpy import *
from scipy import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pylab as plt

def Tangent_Point(t,center, radius):
	
	"""
	t : is the transmitter coordinates 
	center : coordinates of the circle center
	radius :  radius of the circle
	
	""" 
	t = t - center 
	alpha = t[0]
	beta  = t[1]
	d = sqrt(alpha **2+ beta**2)
	delta = 4* alpha**2*radius**2*(beta**2+alpha**2-radius**2)
	print 'delta  = ', delta
	message  = ''
	T1 = array([0,0])
	T2 = array([0,0])
	if delta < 0:
		message = ' the transmitter is in the cercle'
	elif delta == 0 :
		message  = 'the transmitter is on the circle'
	else:
		y1 = (beta*radius**2+alpha*radius*sqrt(d**2-radius**2))/d**2
		y2 = (beta*radius**2-alpha*radius*sqrt(d**2-radius**2))/d**2
		
		if alpha ==0:
			x1 = sqrt(radius**2 - y1**2)
			x2 = -sqrt(radius**2 - y1**2)
		else:
			x1 = (radius**2- beta*y1)/alpha
			x2 = (radius**2- beta*y2)/alpha	

	
	T1 = array([x1,y1]) + center
	T2 = array([x2,y2]) + center
	print 'T1 = ', T1
	print 'T2 = ', T2
	
	return message,T1,T2
	
def Reflection_Point (t,r,center, radius):
	t = t -center
	r = r-center
	alpha = t[0] - r[0]
	beta  = t[1] - r[1]
	message = ''
	R1 = array([0,0])
	R2 = array([0,0])
	if alpha ==0:
		if beta ==0:
			message = 'the transmitter and the receiver are the same '
		else:
			x1 = radius
			x2 = -radius 
			y1 = 0
			y2 = 0 
		
	
	else:
		y1 = sqrt(1/(alpha**2 + beta**2))*alpha*radius 
		y2 = -sqrt(1/(alpha**2 + beta**2))*alpha*radius
		x1 = -beta*y1/alpha
		x2 = -beta*y2/alpha
		
	R1 = array([x1,y1])
	R2 = array([x2,y2])
	print 'R1 = ', R1
	print 'R2 = ', R1 
	return message , R1,R2
	

#
# Lancer de rayons sur cylindre
#
R  = 1
pa = 10*rand(3)
pb = 10*rand(3)
tx = 10*rand(3)
rx = 10*rand(3)
v  = pa-pb
vn = v/sqrt(dot(v,v))
t  = rand(3)
tn = t/sqrt(dot(t,t))
w  = tn-dot(vn,tn)*vn
wn = w/sqrt(dot(w,w))
un = cross(vn,wn)

xs = []
ys = []
zs = []
for alpha in arange(0,2*pi,0.2):
    for t in arange(0,1,0.05):
        tp = pa+t*(pb-pa)+R*cos(alpha)*wn+R*sin(alpha)*un
        xs.append(tp[0])
        ys.append(tp[1])
        zs.append(tp[2])

fig = plt.figure()
ax  = fig.add_subplot(111,projection='3d')
ax.scatter(xs,ys,zs,s=10,marker='o',c='b')
ax.scatter(array([tx[0]]),array([tx[1]]),array([tx[2]]),s=100,marker='^',c='r')
ax.scatter(array([rx[0]]),array([rx[1]]),array([rx[2]]),s=100,marker='^',c='g')
#
# Projection dans le plan du cylindre 
# 
#~ 
#~ ptx = (tx-pa) - dot(tx,vn)*vn
#~ prx = (rx-pa) - dot(rx,vn)*vn

dtu = dot((tx-pa),un)
dtw = dot((tx-pa),wn)
dru = dot((rx-pa),un)
drw = dot((rx-pa),wn)

#~ fig = plt.figure()
#~ plt.plot(R*cos(arange(0,2*pi,0.1)),R*sin(arange(0,2*pi,0.1)))
#~ plt.scatter(dtu,dtw,marker='o',c='r')
#~ plt.scatter(dru,drw,marker='o',c='g')
#~ plt.axis('equal')
ptx= array([dtu,dtw])
prx= array([dru,drw])
msgt,Tt1,Tt2 = Tangent_Point (ptx,array([0,0]), R)
msgr,Tr1,Tr2 = Tangent_Point (prx,array([0,0]), R)
#~ plt.plot([ptx[0],Tt1[0]],[ptx[1],Tt1[1]],'r')
#~ plt.plot([ptx[0],Tt2[0]],[ptx[1],Tt2[1]],'r')

Ct1 = array([Tt1[0],Tt1[1], dot(tx,vn)]) 
Ct2 = array([Tt2[0],Tt2[1], dot(tx,vn)]) 

Cr1 = array([Tr1[0],Tr1[1], dot(rx,vn)]) 
Cr2 = array([Tr2[0],Tr2[1], dot(rx,vn)])

mat  = array([un,wn,vn])
Ct1 = dot(Ct1,mat)
Ct2 = dot(Ct2,mat)
Cr1 = dot(Cr1,mat)
Cr1 = dot(Cr1,mat)

ax.scatter(array([Ct1[0]]),array([Ct1[1]]),array([Ct1[2]]),s=100,marker='o',c='r')
ax.scatter(array([Ct2[0]]),array([Ct2[1]]),array([Ct2[2]]),s=100,marker='o',c='r')

ax.scatter(array([Cr1[0]]),array([Cr1[1]]),array([Cr1[2]]),s=100,marker='o',c='g')
ax.scatter(array([Cr2[0]]),array([Cr2[1]]),array([Cr2[2]]),s=100,marker='o',c='g')

#~ plt.plot([prx[0],Tr1[0]],[prx[1],Tr1[1]],'g')
#~ plt.plot([prx[0],Tr2[0]],[prx[1],Tr2[1]],'g')

#~ plt.plot([tx[0],R1[0]],[tx[1],R1[1]],'-r')
#~ plt.plot([rx[0],R1[0]],[rx[1],R1[1]],'-r')
#~ plt.plot([tx[0],R2[0]],[tx[1],R2[1]],'-g')
#~ plt.plot([rx[0],R2[0]],[rx[1],R2[1]],'-g')


plt.show()
