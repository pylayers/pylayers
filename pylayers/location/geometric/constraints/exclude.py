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
from pylayers.location.geometric.util.boxn import *
from pylayers.location.geometric.constraints.constraint import *


class Exclude(Constraint):		
	"""
	"""
	def __init__(self,nodes):
		"""
		Id : Cell Id
		p  : Cell position 
		v  : Cell direction 
		alpha : sector angle
		Rmax  : Cell Radius (m) 
		"""


		self.std = 1.0
		self.vcw = 3.0
		
		box         = BoxN(np.vstack((np.min(nodes,axis=0),np.max(nodes,axis=0))),ndim=np.shape(nodes)[1]) 
		self.p      = box.ctr
		

		Constraint.__init__(self,'Exclude',self.p)
		self.lbox   = LBoxN([box],ndim=np.shape(self.p)[0])


	def annulus_bound(self):
		"""
		annulus_bound():
		Compute the minimum and maximum distance of the enclosing annulus of the constraint for a given self.vcw
		"""
		pass

	def rescale(self):
		"""
		"""
		pass

	def valid(self,b):
		"""
		valid(b) : check if box b is valid for the given constraint 
		
		A box is valid if it not not valid

		A box is not valid if all distances are greater than rangemax
			           or all distances are less than rangemin				
		"""

		p0 = self.lbox.box[0]
		eps = 0.00000000001
		i=b.intersect(p0)

		if np.sum(i.bd == b.bd) > 5:
			return(True)
		elif i.vol < eps :#(bmax<cmin)|(bmin>cmax): 
			return('out')			
		else :
			return(False)



	def valid_v(self,v):
		"""
		
		"""

		DDbound = []

		DDbound.append(np.sum(v>=self.lbox.bd[0,:],axis=1)>2)
		DDbound.append(np.sum(v<=self.lbox.bd[1,:],axis=1)>2)

		return DDbound


	def limit_valid_v(self,v):
		"""
		
		"""
		return v[np.nonzero( (np.sum(v>=self.lbox.bd[0,:],axis=1)>2) & (np.sum(v<=self.lbox.bd[1,:],axis=1)>2) )]



	def rescale(self,vcw):
		return True


