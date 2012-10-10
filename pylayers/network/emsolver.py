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

import pylayers.util.pyutil as pyu

import pylayers.util.project
from   pylayers.gis.layout import Layout
import pylayers.antprop.slab
from pylayers.antprop.multiwall import *

from   pylayers.network.model import Model


import pdb


class EMSolver(object):
    """ Invoque an electromagnetic solver
    """

    def __init__(self,L=Layout()):

        self.config     = ConfigParser.ConfigParser()
        self.config.read(pyu.getlong('EMSolver.ini','ini'))
#        self.config.read(pkgutil.get_loader('pylayers').filename +'/ini/EMSolver.ini')
        self.ems_opt = dict(self.config.items('EMS_config'))
        self.toa_opt = dict(self.config.items('TOA'))
        self.plm_opt = dict(self.config.items('PL_MODEL'))
        self.EMS_method = self.ems_opt['method']
        self.sigmaTOA      = float(self.toa_opt['sigmatoa']) # meters !!!!!!
        self.sigmaRSS    = float(self.plm_opt['sigmarss'])# dBm !!!!!!                          
        self.f             = float(self.plm_opt['f'])
        self.RSSnp        = float(self.plm_opt['rssnp'])
        self.d0               = float(self.plm_opt['d0'])
        self.PL_method     = self.plm_opt['method'] # mean, median , mode
        self.L=L


    def solve(self,p,e,LDP):
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

    
#        if self.method == 'direct':

#            d =np.linalg.norm(n1p-n2p)
#            
#            if LDP == 'TOA':
#                sp.random.seed(0)
#                std = self.sigmaTOA*sp.randn()
#                return ((d + std)*0.3,self.sigmaTOA)

#            elif LDP == 'Pr':
#                sp.random.seed(0)
#                std = self.sigmaRSS*sp.randn()
#                M = Model.Model(method='mean')
#                r,rssstd=M.imeasure(d,self.sigmaRSS)
#                return (r,rssstd)
#                
#            else :
#                raise NameError('invlaid LDP name')
        #sp.random.seed(0)
        if self.EMS_method == 'direct':
            dd={} # distance dictionnary

            if len(e) > 0:
                lp=np.array([np.array((p[e[i][0]],p[e[i][1]])) for i in range(len(e))])
                d=np.sqrt(np.sum((lp[:,0]-lp[:,1])**2,axis=1))
            #    dd=dict(zip(e,d))
                if LDP == 'TOA':
                    std = self.sigmaTOA*sp.randn(len(d))
                    return ([[max(0.0,(d[i]+std[i])*0.3),self.sigmaTOA*0.3] for i in range(len(d))],d)

                elif LDP == 'Pr':
                    std = self.sigmaRSS*sp.randn(len(d))
                    M = Model(method=self.PL_method,f=self.f,RSSnp=self.RSSnp,d0=self.d0)
                    r=M.getPL(d,self.sigmaRSS)
                    return ([[r[i],self.sigmaRSS] for i in range(len(d))],d)
                
            



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
                    pdb.set_trace()
                    pa=np.vstack(p.values())
                    for px in pa:
                        Lwo,Lwp=Loss0_v2(self.L,pa,self.f,px)
#                    std = self.sigmaRSS*sp.randn(len(d))
#                    M = Model(method=self.PL_method,f=self.f,RSSnp=self.RSSnp,d0=self.d0)
#                    r=M.getPL(d,self.sigmaRSS)
#                    return ([[r[i],self.sigmaRSS] for i in range(len(d))],d)
            
                elif LDP == 'TOA': #### NOT CORRECT !
                    std = self.sigmaTOA*sp.randn(len(d))
                    return ([[max(0.0,(d[i]+std[i])*0.3),self.sigmaTOA*0.3] for i in range(len(d))],d)


            else :
                return ([[0.],[0.]])


        elif self.method == 'Pyray':
            print 'Okay, I think we\'ve got something to append in the TODO list'


        else :
            raise NameError('invalid method name')






