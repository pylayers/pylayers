# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of RGPA.

#Foobar is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#Foobar is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
#authors :
#Nicolas AMIOT		: nicolas.amiot@univ-rennes1.fr
#Bernard UGUEN		: bernard.uguen@univ-rennes1.fr
#Mohamed LAARAIEDH	: mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
import numpy as np
import scipy as sp 
import time
import sys
#sys.path.append('./Geometric_Loc/Util') 
from pylayers.location.geometric.util.BoxN import *
#from BoxN import *



class Constraint(object):
	""" Constraint class

	

	attributes
	----------

	Constraint.C_Id : Contraint Identity automatically increased for new instanciation

	type	: contraint type
	time	: time stamp
	lbox	: LBoxN intialisation
	runable : boolean. information on constraint center ( TO DO)
	p	: Constraint center
	validity: validity duration ( not use yet)
	self.Id	: Constraint.C_Id
	ndim	: Constraint dimension


	parmsh : dictionary
		keys :	['display']=True     # launch geomview interactively 
			['mode']='Current'   # display current box or full constraint
			['scene']=True      # display whole scene
			['point']=True       # display constraint point(s)
			['boxes']=True       # display constraint box
			['estimated']=True  # display estimated point
			['estimated_LS']=False # display estimated point with LS method
			['quadric']=True   # display sphere or hyperbola 
			['grid']=True       # display grid 
			['grav']=True 


	methods
	-------
		
	info()     : Display information about constraint
	show3()    : display constraint on Geomview. Parameters tahnks to self.parmsh


	TODO
	----
	ini file for geom view !

	"""
	C_Id = 0   # Global counter of existing constraints 

	def __init__(self,type,p=0):
		
		self.type       = type
		self.time       = time.time()
		self.lbox       = LBoxN([])
		self.evaluated  = False
		if (isinstance(p,np.ndarray)):
			self.runable = True
			self.p = p 
		else:	
			self.runable = False
		self.validity   = 20                  # not used 
		self.Id         = Constraint.C_Id  # Id of constraint is set to counter + 1
		self.ndim       = len(p.T)
		#
		# Parameters for show3 with geomview 
		#
 		self.parmsh={}
		self.parmsh['display']=True     # launch geomview interactively 
		self.parmsh['mode']='Current'   # display current box or full constraint
		self.parmsh['scene']=True      # display whole scene
		self.parmsh['point']=True       # display constraint point(s)
		self.parmsh['boxes']=True       # display constraint box
		self.parmsh['estimated']=True  # display estimated point
		self.parmsh['estimated_LS']=False # display estimated point with LS method
		self.parmsh['quadric']=True   # display sphere or hyperbola 
		self.parmsh['grid']=True       # display grid 
		self.parmsh['grav']=True       # display box gravity center 




	def info(self):
		""" display info on constraint
		"""
		print "Type         : ",self.type
		print "--------------------------"
		print "Time         : ",self.time
		print "validity (s) : ",self.validity
		
		if ((self.runable)):
			print "Origin : ",self.p
		
		if self.evaluated:
			Npts = np.shape(self.g.p)[0]
			print "Nb valid points in volume    : ",Npts," voxel"
			print "Taille kO: ",Npts*12/(2**10)  ," kO"

		if self.type=="TOA":
			self.estvol()
			print "Volume Estimatif",self.estvlm
		 	print "Toa (ns)", self.value	
		 	print "std (ns)", self.std
		 	print "vcw     ", self.vcw
		 	print "Range(m)", self.range	
		 	print "sstd (m)", self.sstd
		print "-------------------"
		self.lbox.info()
		print "-------------------"



	def show3(self):
		""" display constraint on Geomview
		
		The filename is boxes{Id}.list 

		Id is the Id of the current constraint



		"""
		fname    = 'boxes'+str(self.Id)
		filename = "./geom/"+fname+".list"
		fd = open(filename,"w")
		fd.write("LIST\n")

		#
		# Display scene
		#
		if self.parmsh['scene']:
			fd.write("{<scene.list}\n")

		#
		# Display boxes
		#

		if self.parmsh['mode']=='Full':
			H_Id=self.history[-1].Id+1
			#
			# Display all boxes
			#
			for k in range(len(self.history)):		
				cons = self.history[k]   # constraint k 
				lb   = cons.lbox
				lb.parmsh['display']=False
				filename2=lb.show3(Id=cons.Id)
				fd.write("{<"+filename2+"}\n")
		elif self.parmsh['mode']=='Current':
			#
			# Display current box
			#
			color = ['m','g','c','y','m','b','r','m','g','c','y','orange','skyblue']	
			#color = ['skyblue','skyblue','orange']
				
			lb = self.lbox
			lb.parmsh['display']=False
#			if self.type=='TOA':
#				filename2 = lb.show3(Id=[self.Id],col='skyblue')
#			elif self.type=='RSS':
#				filename2 = lb.show3(Id=[self.Id],col='orange')
#			elif self.type=='TDOA':
#				filename2 = lb.show3(Id=[self.Id],col='c')	
#			else :
#				filename2 = lb.show3(Id=[self.Id],col='y')	
			filename2 = lb.show3(Id=[self.Id],col=color[self.Id])				
			fd.write("{<"+filename2+"}\n")
			#
			# Display Spherical Constraint
			#
			if self.parmsh['quadric']:
				if self.type=='TOA':
					c1 = str(self.range+self.vcw*self.sstd)
					c2 = str(max(0,self.range-self.vcw*self.sstd))
					try:
						c3 = str(self.p[0])+" "+str(self.p[1])+" "+str(self.p[2])
					except : 
						c3 = str(self.p[0])+" "+str(self.p[1])+" "+str(0)
					fd.write("{appearance {-edge  patchdice	10 10 material {alpha 0.2}} {SPHERE "+ c1 +" "+ c3+ " }}\n")
					fd.write("{appearance {-edge  patchdice	10 10 material {alpha 0.2}} {SPHERE "+ c2 +" "+ c3+ " }}\n")
		fd.close()		
		#
		# Display points
		#
		if self.parmsh['point']:
			if self.type!='Fusion':
				g.cloud(self.p,name=fname,color='g',dice=6,access='append')

		if self.evaluated:
			if self.parmsh['estimated']:
				g.cloud(self.pe,name=fname,color='b',dice=6,access='append')
			
			if self.parmsh['estimated_LS']:
				g.cloud(self.p_LS,name=fname,color='r',dice=6,access='append')

			if self.parmsh['grid']:
				if self.evaluated:
					g.cloud(self.g.p,name=fname,color='k',dice=2,R=0.1,access='append')

		if self.parmsh['display']:
			chaine = "geomview  -nopanel  -b 1 1 1 " + filename + " 2>/dev/null &"
			os.system(chaine)
		else:
			return(fname)
	
