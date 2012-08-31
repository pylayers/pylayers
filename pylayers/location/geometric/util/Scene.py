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
import pdb
import Geomview as g




class Scene(object):
	"""
	Scene creator for Geomview display


	:Attributes:
		an	: np.array
			anchors nodes
		bn	: np.array
			blind nodes

	:Methods:
		generate :
			create scene

	"""
	def __init__(self,an,bn):
		self.an  = an
		self.bn  = bn		

	def generate(self):
		g.cloud(self.an,name='scene',color='r',dice=6,R=1,access='new')
		try :
			ll=np.shape(self.bn)[1]
			for i in range (ll):
				g.cloud(self.bn[i,:],display=False,name='scene',color='g',dice=6,access='append')
		except :
			g.cloud(self.bn,display=False,name='scene',color='g',dice=6,access='append')

