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
import pylayers.antprop.slab

from pylayers.gis.layout import Layout
import pylayers.antprop.loss as lo
from pylayers.util.project import *
from pylayers.network.model import PLSmodel


import pdb


class EMSolver(object):
    """ Invoque an electromagnetic solver

    Attributes
    ----------

    config
    fileini
    ems_opt
    toa_opt
    EMS_method
       
    sigmaTOA

    Notes
    -----

    During dynamic simulation, on a regular basis nodes needs radio information
    about their neighbors. This radio information can be obtained with
    different technique from statistical model to site specific simulation
    using ray tracing.

    EMS_method = { 'multiwall','raytracing'}



    """

    def __init__(self,L=Layout()):
        self.config  = ConfigParser.ConfigParser()
        self.fileini ='EMSolver.ini'
        self.config.read(pyu.getlong(self.fileini,pstruc['DIRSIMUL']))
        self.ems_opt = dict(self.config.items('EMS_config'))
        self.toa_opt = dict(self.config.items('TOA'))

        self.EMS_method = self.ems_opt['method']
        self.sigmaTOA = float(self.toa_opt['sigmatoa']) # meters !!!!!!

        self.model = {}

        self.L = L


    def save_model(self,RAT,model):
        """save a RAT model

        Parameters
        ----------

        RAT : string
            RAT name
        model : dictionnary
            parameters of the PL model
           
           
        """

        fileini = pyu.getlong(self.fileini, pstruc['DIRSIMUL'])
        fd = open(fileini, "a")
        nconfig     = ConfigParser.ConfigParser()
        nconfig.add_section(RAT+'_PLM')
        nconfig.set(RAT+'_PLM','sigrss', str(model.sigrss))
        nconfig.set(RAT+'_PLM','f', str(model.f))
        nconfig.set(RAT+'_PLM','rssnp', str(model.rssnp))
        nconfig.set(RAT+'_PLM','d0', str(model.d0))
        nconfig.set(RAT+'_PLM','method', str(model.method))
        nconfig.write(fd)
        fd.close()



    def load_model(self,RAT):
        """ load a path loss shadowing model for a RAT

        Parameters
        ----------

        RAT  : string
            RAT name


        """

        ratopt = dict(self.config.items(RAT+'_PLM'))
       
        # Path Loss Shadowing model
        self.model[RAT] = PLSmodel(f = eval(ratopt['f']),
                                   rssnp = eval(ratopt['rssnp']),
                                   d0 = eval(ratopt['d0']),
                                   sigrss = eval(ratopt['sigrss']),
                                   method = ratopt['method'])


    def solve(self,p,e,LDP,RAT,epwr,sens):
        """ computes and returns a LDP value a given method

        Parameters
        ----------
       
        p : np.array
        e : np.array
        LDP : string
            Type of LDP ( TOA, Pr, .... any other are to be add in teh todo list)
        epwr : list
           nodes emmited power
        sens : list
            nodes sensitivity
               
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
                model = self.model[RAT]
            except:
                self.model[RAT] = PLSmodel()
                self.save_model(RAT,self.model[RAT])
                model = self.model[RAT]


#        if self.EMS_method == 'direct':
#            dd={} # distance dictionnary

#            if len(e) > 0:
#                lp=np.array([np.array((p[e[i][0]],p[e[i][1]])) for i in range(len(e))])
#                d=np.sqrt(np.sum((lp[:,0]-lp[:,1])**2,axis=1))
#                if LDP == 'TOA':
#                    std = self.sigmaTOA*sp.randn(len(d))
#                    return ([[max(0.0,(d[i]+std[i])*0.3),self.sigmaTOA*0.3] for i in range(len(d))],d)

#                elif LDP == 'Pr':
#                    std = self.model.sigrss*sp.randn(len(d))
#                    r=model.getPL(d,model.sigrss)
#                    return ([[- r[i]-model.PL0,model.sigrss] for i in range(len(d))],d)
#               
#           



#                else :
#                    raise NameError('invalid LDP name')
#            else :
#                return ([[0.],[0.]])


        if self.EMS_method == 'multiwall':

            dd={} # distance dictionnary
            if len(e) > 0:

                lp = np.array([np.array((p[e[i][0]],p[e[i][1]])) for i in range(len(e))])
                # MW is 2D only now.
                # This explain the following licenses
                lp = lp[:,:,:2]
                dim = lp.shape[2]
                # euclidian distance
                d = np.sqrt(np.sum((lp[:,0]-lp[:,1])**2,axis=1))
                slp = np.shape(lp)[1]

                # evaluation of all LDPs
                if LDP=='all':
                    pa = np.vstack(p.values())
                    lpa = len(pa)
                    Pr = []
                    TOA = []
                    lsens = np.array(())
                    loss = np.array(())
                    frees = np.array(())
                    lepwr = np.array(())

                    for i in range(lpa-1):
                        # excess time of flight + losses computation
                        #
                        # Losst returns 4 parameters
                        #   Lo Lp Edo Edp
                        #
                        #
                        MW = lo.Losst(self.L,model.f,pa[i+1:lpa,:dim].T,pa[i,:dim])
                        # MW = lo.Loss0_v2(self.L,pa[i+1:lpa],model.f,pa[i])
                        # loss free space

                        frees=np.hstack((frees,lo.PL(np.array([model.f]),pa[i+1:lpa,:dim].T,pa[i,:dim].reshape(2,1),model.rssnp)[0] ))
                        # Pr.extend(lepwr - MW[0] - frees)
                        # WARNING : only one polarization is taken into
                        # account here
                        # save losses computation

                        loss = np.hstack((loss,MW[0][0]))
                        # save excess tof computation
                        TOA  = np.hstack((TOA,MW[2][0]))

                    # emmited power for the first nodes of computed edges
                    lepwr1 = [epwr[i[0]][RAT] for i in e]
                    lepwr2 = [epwr[i[1]][RAT] for i in e]
                    Pr = lepwr1 - loss - frees

                    # concatenate reverse link
                    Pr = np.hstack((Pr, lepwr2 - loss - frees))
                    P = np.outer(Pr,[1,1])
                    P[:,1] = model.sigrss
                    lsens = [sens[i[0]][RAT] for i in e] + [sens[i[1]][RAT] for i in e]

                    # visibility or not
                    v = P[:,0] > lsens

                    # same toa for link and reverse link
                    T = np.hstack((TOA+d/0.3,TOA+d/0.3))
                    T=np.outer(T,[1,1])
                    T[:,1]=self.sigmaTOA*0.3
                    d=np.hstack((d,d))
                    return (P,T,d,v)




#                elif LDP == 'Pr':
#                    pa = np.vstack(p.values())
#                    pn = p.keys()
#                    lpa = len(pa)
#                    Lwo = []
#                    frees=[]
#                    lepwr=[]
#                    for i in range(lpa-1):
#                        lo.append(Loss0_v2(self.L,pa[i+1:lpa],model.f,pa[i]))
#                        Lwo.extend(Loss0_v2(self.L,pa[i+1:lpa],model.f,pa[i])[0])
#                        frees.extend(PL(pa[i+1:lpa],model.f,pa[i],model.rssnp))
#                        lepwr.extend(epwr[i+1:lpa])
#                    return ([[lepwr[i] - Lwo[i]-frees[i],model.sigrss] for i in range(len(Lwo))],d)
#           
#                elif LDP == 'TOA': #### NOT CORRECT !
#                    std = self.sigmaTOA*sp.randn(len(d))
#                    return ([[max(0.0,(d[i]+std[i])*0.3),self.sigmaTOA*0.3] for i in range(len(d))],d)


            else :
                return (np.array((0.,0.)),np.array((0.,0.)),np.array((0.,0.)),np.array((0.,0.)))


        elif self.method == 'raytracing':
            print 'Okay, I think we\'ve got something to append in the TODO list'


        else :
            raise NameError('invalid method name')






