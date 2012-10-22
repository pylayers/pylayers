# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of Network.

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
#Nicolas AMIOT        : nicolas.amiot@univ-rennes1.fr
#Bernard UGUEN        : bernard.uguen@univ-rennes1.fr
#####################################################################
import os
import pkgutil

import numpy as np
import scipy as sp
import ConfigParser
import itertools

import pylayers.util.pyutil as pyu

import pylayers.util.project
from   pylayers.gis.layout import Layout
import pylayers.antprop.slab
from pylayers.antprop.multiwall import *

from pylayers.util.project import *
from   pylayers.network.model import Model


import pdb


class EMSolver(object):
    """ Invoque an electromagnetic solver
    """

    def __init__(self,L=Layout()):

        self.config     = ConfigParser.ConfigParser()
        self.fileini='EMSolver.ini'
        self.config.read(pyu.getlong(self.fileini,pstruc['DIRSIMUL']))
        self.ems_opt = dict(self.config.items('EMS_config'))
        self.toa_opt = dict(self.config.items('TOA'))
        self.rss_opt = dict(self.config.items('RSS'))


        self.EMS_method = self.ems_opt['method']
        self.sigmaTOA      = float(self.toa_opt['sigmatoa']) # meters !!!!!!

        self.model={}
        self.method     = self.rss_opt['method'] # mean, median , mode

        self.L=L

    
    def save_model(self,RAT,model):
        """save a RAT model
            
        """
        fileini = pyu.getlong(self.fileini, pstruc['DIRSIMUL'])
        fd = open(fileini, "a")
        nconfig     = ConfigParser.ConfigParser()
        nconfig.add_section(RAT+'_PLM')
        nconfig.set(RAT+'_PLM','sigrss', str(model.sigrss))
        nconfig.set(RAT+'_PLM','f', str(model.f))
        nconfig.set(RAT+'_PLM','rssnp', str(model.rssnp))
        nconfig.set(RAT+'_PLM','d0', str(model.d0))
        nconfig.write(fd)
        fd.close()



    def load_model(self,RAT):
        ratopt=dict(self.config.items(RAT+'_PLM'))
        self.model[RAT]=Model(f=eval(ratopt['f']),rssnp=eval(ratopt['rssnp']),d0=eval(ratopt['d0']),sigrss=eval(ratopt['sigrss']),method=self.method)

    def solve(self,p,e,LDP,RAT):
        """compute and return a LDP value thanks to a given method

        Attributes
        ----------
            n1p : np.array
                node 1 postion
            n2p : np.array
                node 2 postion
            LDP : string
                Type of LDP ( TOA, Pr, .... any other are to be add in teh todo list)
        
        Returns
        -------
            value : float
                A LDP value :     * A time in ns for LDP ='TOA'
                        * A received power in dBm for LDP ='Pr'
            std : float
                A LDP value standard deviation:     * A time in ns for LDP ='TOA'
                                    * A received power in dBm for LDP ='Pr'

        """
        try:
            model= self.model[RAT]
        except:
            try:
                self.load_model(RAT)
                model=self.model[RAT]
            except:
                self.model[RAT]=Model()
                self.save_model(RAT,self.model[RAT])
                model=self.model[RAT]


        if self.EMS_method == 'direct':
            dd={} # distance dictionnary

            if len(e) > 0:
                lp=np.array([np.array((p[e[i][0]],p[e[i][1]])) for i in range(len(e))])
                d=np.sqrt(np.sum((lp[:,0]-lp[:,1])**2,axis=1))
                if LDP == 'TOA':
                    std = self.sigmaTOA*sp.randn(len(d))
                    return ([[max(0.0,(d[i]+std[i])*0.3),self.sigmaTOA*0.3] for i in range(len(d))],d)

                elif LDP == 'Pr':
                    std = self.model.sigrss*sp.randn(len(d))
                    r=model.getPL(d,model.sigrss)
                    return ([[r[i]-model.PL0,model.sigrss] for i in range(len(d))],d)
                
            



                else :
                    raise NameError('invalid LDP name')
            else :
                return ([[0.],[0.]])


        elif self.EMS_method == 'multiwall':

            dd={} # distance dictionnary

            if len(e) > 0:
                lp=np.array([np.array((p[e[i][0]],p[e[i][1]])) for i in range(len(e))])
                d=np.sqrt(np.sum((lp[:,0]-lp[:,1])**2,axis=1))
#                if LDP == 'TOA':
#                    std = self.sigmaTOA*sp.randn(len(d))
#                    return ([[max(0.0,(d[i]+std[i])*0.3),self.sigmaTOA*0.3] for i in range(len(d))],d)

                if LDP == 'Pr':

                    pa = np.vstack(p.values())
                    pn = p.keys()
                    lpa = len(pa)
                    Lwo = []
                    for i in range(lpa-1):
                        Lwo.extend(Loss0_v2(self.L,pa[i+1:lpa],model.f,pa[i])[0])
                    return ([[-Lwo[i]-model.PL0,model.sigrss] for i in range(len(Lwo))],d)
            
                elif LDP == 'TOA': #### NOT CORRECT !
                    std = self.sigmaTOA*sp.randn(len(d))
                    return ([[max(0.0,(d[i]+std[i])*0.3),self.sigmaTOA*0.3] for i in range(len(d))],d)


            else :
                return ([[0.],[0.]])


        elif self.method == 'Pyray':
            print 'Okay, I think we\'ve got something to append in the TODO list'


        else :
            raise NameError('invalid method name')






