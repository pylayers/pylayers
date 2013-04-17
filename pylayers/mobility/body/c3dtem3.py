import c3d as c3d
import numpy as np
import scipy.stats as sp
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D

import pdb





file_name = '07_01.c3d' 
h = c3d.Header(open(file_name))
r = c3d.Reader(open(file_name))
frm=[]
var = r.read_frames()
i=0
sizes = [1200,1200,1200]
for pnt, ang in var:
	print i
	frm.append((pnt, ang))
	x = pnt[:,0]
	y = pnt[:,1]
	z = pnt[:,2]
	if not x.any() and not y.any() and not z.any():
		print x
		print y
		print z
	'''print x
	print y
	print z'''
	fig = plt.figure()
	ax  = Axes3D(fig)
	ax.plot3D(x,y,z, 'bo',markersize = 5)
	plt.xlabel('X')
	plt.ylabel('Y')
	#~ ax.set_xlim3d([0,250])
	#~ ax.set_ylim3d([0,250])
	#~ ax.set_zlim3d([0,180])
	#~ plt.plot(x,z,'o')
	#~ plt.axis([-2400,2400,-200,1800])
	#~ plt.savefig('Figures/'+str(i)+".png", format='png')
	plt.show()

			

