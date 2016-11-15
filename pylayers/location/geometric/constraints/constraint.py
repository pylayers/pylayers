# -*- coding:Utf-8 -*-
"""

Constraint Class
================

.. autosummary::
    :toctree: generated/

    Constraint.__init__
    Constraint.__repr__
    Constraint.updc
    Constraint.info
    Constraint.show3

"""
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
#Nicolas AMIOT          : nicolas.amiot@univ-rennes1.fr
#Bernard UGUEN          : bernard.uguen@univ-rennes1.fr
#Mohamed LAARAIEDH      : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
from pylayers.util.project import *
import numpy as np
import scipy as sp
import time
import sys
from pylayers.location.geometric.util.boxn import *
from pylayers.location.geometric.util import geomview as g


class Constraint(object):
    """ Constraint class



    Attributes
    ----------

    Constraint.C_Id : Contraint Identity automatically increased for new instanciation

    type    : str
        contraint type (TOA,TDOA,RSS)

    p       : nd.array
        Constraint center

    self.Id : str
        Constraint.C_Id

    ndim    : int
        Constraint dimension (2D/3D)

    origin : dict
        origin of the contraint, (used for simulnet simulation 
        'id' : id of the node generating the constraint
        'link' : edge of network graph Gr
        'rat' : RAT on which constraint has been obtained
        'ldp' : type of observable


    Parameters
    ----------

    time    : float
        time stamp
    lbox    : LBoxN
        LBoxN intialisation
        self.usable = True
        self.visible = True
        self.obsolete = False

    runable : boolean
        A constriat is runable if it has a center
    usable: boolean
        a constraint is usable is it has a std AND a value AND is visible
    visible : boolean
        is the constraint visible (used for simulnet simulation 


    parmsh : dictionary
            keys :  ['display']=True     # launch geomview interactively
                    ['mode']='Current'   # display current box or full constraint
                    ['scene']=True      # display whole scene
                    ['point']=True       # display constraint point(s)
                    ['boxes']=True       # display constraint box
                    ['estimated']=True  # display estimated point
                    ['estimated_LS']=False # display estimated point with LS method
                    ['quadric']=True   # display sphere or hyperbola
                    ['grid']=True       # display grid
                    ['grav']=True


    Methods
    -------

    info()     : Display information about constraint
    show3()    : display constraint on Geomview. Parameters tahnks to self.parmsh

    See Also
    --------

    pylayers.simul.simulnet
    pylayers.location.localization


    """
    def __init__(self, type, id='0', p=np.array(()), origin={}):

        self.type = type
        self.time = time.time()
        self.lbox = LBoxN([])
        self.evaluated = False
        if p.any():
            self.runable = True
        else:
            self.runable = False
        self.usable = True
        self.visible = True
        self.obsolete = False
        self.p = p
        self.validity = 20                  # not used
        self.id = id  # Id of constraint is set to counter + 1
        self.ndim = len(p.T)
        self.origin = origin
        #
        # Parameters for show3 with geomview
        #
        self.parmsh = {}
        self.parmsh['display'] = True     # launch geomview interactively
        self.parmsh['mode'] = 'Current'   # display current box or full constraint
        self.parmsh['scene'] = False      # display whole scene
        self.parmsh['point'] = True       # display constraint point(s)
        self.parmsh['boxes'] = True       # display constraint box
        self.parmsh['estimated'] = True  # display estimated point
        self.parmsh['estimated_LS'] = False  # display estimated point with LS method
        self.parmsh['quadric'] = True   # display sphere or hyperbola   
        self.parmsh['grid'] = True       # display grid
        self.parmsh['grav'] = True       # display box gravity center
        if self.origin == {}:
            if self.type == 'RSS':
                typ = 'Pr'
            else :
                typ =self.type
            self.origin = {'id':self.id, 'ldp':typ, 'link':'','rat':''}

    def __repr__(self):
        np.set_printoptions(precision=3)
        s = '{0:4} | {1:6} |{2:4} | {3:4} | {4:15}| {5:9}| {6:5}| {7:7}| {8:6}|'.format('node','peer','type', 'rat', 'p', 'value', 'std', 'runable' , 'usable' )

        node = self.origin['id']
        peer = self.origin['link']
        rat  = self.origin['rat']
        if self.type != 'TDOA':
            s = s + '\n' + '{0:4} | {1:6} |{2:4} | {3:4} | {4:15}| {5:9}| {6:5}| {7:7}| {8:6}|'.format(node,peer,self.type,rat, self.p, self.value, self.std, self.runable, self.usable)
        else:
            s = s + '\n' + '{0:4} | {1:6} |{2:4} | {3:4} | {4:15}| {5:9}| {6:5}| {7:7}| {8:6}|'.format(node,peer,self.type,rat, self.p[0], self.value, self.std, self.runable, self.usable)
            s = s + '\n' + '                            '+str(self.p[1])

        return s
        # np.set_printoptions(precision=3)
        # s = '{0:4} | {1:15}| {2:9}| {3:5}| {4:7}| {5:6}| {6:8}| {7:9}'.format('type', 'p', 'value', 'std', 'runable' , 'usable' , 'obsolete' , 'evaluated')
        # if self.type != 'TDOA':
        #     s = s + '\n' + '{0:4} | {1:15}| {2:9}| {3:5}| {4:7}| {5:6}| {6:8}| {7:9}'.format(self.type, self.p, self.value, self.std, self.runable, self.usable , self.obsolete , self.evaluated)
        # else:
        #     s = s + '\n' + '{0:4} | {1:15}| {2:9}| {3:5}| {4:7}| {5:6}| {6:8}| {7:9}'.format(self.type, self.p[0], self.value, self.std, self.runable, self.usable , self.obsolete , self.evaluated)
        #     s = s + '\n' + '       '+str(self.p[1])

        return s

    def updc(self,name='p',value=np.array(())):
        """ update values of a constraint

        Example :
        ---------

        >>> from pylayers.location.geometric.constraints.toa import *
        >>> T=TOA()
        >>> T.usable
        True
        >>> T.p
        array([], dtype=float64)


        >>> T.updc('p',np.array((10,3)))
        >>> T.p
        array([10,  3])
        >>> T.updc('value',[np.nan])
        >>> T.usable
        False

        """
        if not isinstance(value,np.ndarray):
            value = np.array([value])

        if np.sum(np.isnan(value))<1 and value.size>0:
            setattr(self,name,value)

            # once p is set, constraint is runable
            if name =='p':
                self.runable = True

            # once value is set and constraint has a std value
            # constraint is usable
            if name =='value':
                if np.sum(np.isnan(self.std))<1 and self.runable:
                    self.usable = True

            # once std is set and constraint has a value
            # constraint is usable
            if name =='std':
                if np.sum(np.isnan(self.value))<1 and self.runable:
                    self.usable = True


        else :
            if name =='p':
                self.runable = False
            elif name =='value' or name =='std':
                self.usable = False


    def info_old(self):
        """ display info on constraint
        """
        print "Type         : ", self.type
        print "--------------------------"
        print "Time         : ", self.time
        print "validity (s) : ", self.validity

        if ((self.runable)):
            print "Origin : ", self.p

        if self.evaluated:
            Npts = np.shape(self.g.p)[0]
            print "Nb valid points in volume    : ", Npts, " voxel"
            print "Taille kO: ", Npts * 12 / (2 ** 10), " kO"

        if self.type == "TOA":
            self.estvol()
            print "Estimated Volume", self.estvlm
            print "Toa (ns)", self.value
            print "std (ns)", self.std
            print "vcw     ", self.vcw
            print "Range(m)", self.range
            print "sstd (m)", self.sstd
        print "-------------------"
        self.lbox.info()
        print "-------------------"

    # def info(self):


    #     print '{0:3} , {1:10}, {2:10}, {3:7}, {4:1}, {5:1}, {6:1}, {7:1}'.format('type', 'p', 'value', 'std', 'runable' , 'usable' , 'obsolete' , 'evaluated')
    #     np.set_printoptions(precision=3)
    #     print '{0:3} , {1:10}, {2:10}, {3:7}, {4:1}, {5:1}, {6:1}, {7:1}'.format(self.type, self.p, self.value, self.std, self.runable, self.usable , self.obsolete , self.evaluated)
    def info(self):
        """

        """
        print '{0:4} , {1:15}, {2:5}, {3:5}, {4:7}, {5:6}, {6:8}, {7:9}'.format('type', 'p', 'value', 'std', 'runable' , 'usable' , 'obsolete' , 'evaluated')
        np.set_printoptions(precision=3)
        print '{0:4} , {1:15}, {2:5}, {3:5}, {4:7}, {5:6}, {6:8}, {7:9}'.format(self.type, self.p, self.value, self.std, self.runable, self.usable , self.obsolete , self.evaluated)


    def _show3(self):


        if self.runable:
            

            #
            # Display scene
            #
            # if self.parmsh['scene']:

            if self.parmsh['mode'] == 'Current':
                #
                # Display current box
                #
                color = ['m', 'g', 'c', 'y', 'm', 'b', 'r',
                         'm', 'g', 'c', 'y', 'orange', 'skyblue']
                #color = ['skyblue','skyblue','orange']
                if self.parmsh['boxes']:
                    lb = self.lbox
                    lb.parmsh['display'] = False
                    filename2 = lb._show3(Id=[self.id], col='m')  # )color[self.Id])
                #
                # Display Spherical Constraint
                #
                if self.parmsh['quadric']:
                    if self.type == 'TOA':
                        # Create inner and outer sphere
                        ro = self.range + self.vcw * self.sstd
                        ri = max(0, self.range - self.vcw * self.sstd)
                        pi = np.pi
                        cos = np.cos
                        sin = np.sin
                        phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
                        xi = ri*sin(phi)*cos(theta)
                        yi = ri*sin(phi)*sin(theta)
                        zi = ri*cos(phi)
                        xo = ro*sin(phi)*cos(theta)
                        yo = ro*sin(phi)*sin(theta)
                        zo = ro*cos(phi)

                        si = mlab.mesh(xi+self.p[0] , yi+self.p[1] , zi+self.p[2],opacity=0.1,color=(0,1,0.5))
                        so = mlab.mesh(xo+self.p[0] , yo+self.p[1] , zo+self.p[2],opacity=0.1,color=(0,1,1))

            #
            # Display points
            #
            if self.parmsh['point']:
                if self.type != 'Fusion':
                    mlab.points3d(self.p[0],self.p[1],self.p[2],color=(1,1,1),opacity=1)
                    # mlab.text3d(self.p[0],self.p[1],self.p[2],'p',color=(0,0,0))
                
            if self.evaluated:
                if self.parmsh['estimated']:
                    mlab.points3d(self.pe[0],self.pe[1],self.pe[2],color=(0,0,1))
                    mlab.text3d(self.pe[0],self.pe[1],self.pe[2],'pe',color=(0,0,0))
                if self.parmsh['estimated_LS']:
                    mlab.points3d(self.p_LS[0],self.p_LS[1],self.p_LS[2],color=(1,0,1))
                    mlab.text3d(self.p_LS[0],self.p_LS[1],self.p_LS[2],'p_LS',color=(0,0,0))
        else:
            print 'constraint is not runnable. It can not be displayed'



    def show3(self):
        """ display constraint on Geomview

        The filename is boxes{Id}.list

        Id is the Id of the current constraint



        """
        if self.runable:
            fname = 'boxes' + str(self.id)
            filename = basename + "/geom/" + fname + ".list"
            fd = open(filename, "w")
            fd.write("LIST\n")

            #
            # Display scene
            #
            if self.parmsh['scene']:
                fd.write("{<scene.list}\n")

            #
            # Display boxes
            #

    #               if self.parmsh['mode']=='Full':
    #                       H_Id=self.history[-1].Id+1
    #                       #
    #                       # Display all boxes
    #                       #
    #                       for k in range(len(self.history)):
    #                               cons = self.history[k]   # constraint k
    #                               lb   = cons.lbox
    #                               lb.parmsh['display']=False
    #                               filename2=lb.show3(Id=cons.Id)
    #                               fd.write("{<"+filename2+"}\n")
            elif self.parmsh['mode'] == 'Current':
                #
                # Display current box
                #
                color = ['m', 'g', 'c', 'y', 'm', 'b', 'r',
                         'm', 'g', 'c', 'y', 'orange', 'skyblue']
                #color = ['skyblue','skyblue','orange']
                if self.parmsh['boxes']:
                    lb = self.lbox
                    lb.parmsh['display'] = False
                    filename2 = lb.show3(Id=[self.id], col='m')  # )color[self.Id])
                    fd.write("{<" + filename2 + "}\n")
                #
                # Display Spherical Constraint
                #
                if self.parmsh['quadric']:
                    if self.type == 'TOA':
                        c1 = str(self.range + self.vcw * self.sstd)
                        c2 = str(max(0, self.range - self.vcw * self.sstd))
                        try:
                            c3 = str(self.p[0]) + " " + str(self.p[1]) + " " + \
                                str(self.p[2])
                        except:
                            c3 = str(self.p[0]) + " " + str(self.p[
                                1]) + " " + str(0)
                        fd.write("{appearance {-edge  patchdice 10 10 material {alpha 0.2}} {SPHERE " + c1 + " " + c3 + " }}\n")
                        fd.write("{appearance {-edge  patchdice 10 10 material {alpha 0.2}} {SPHERE " + c2 + " " + c3 + " }}\n")
            fd.close()
            #
            # Display points
            #
            if self.parmsh['point']:
                if self.type != 'Fusion':
                    g.cloud(self.p, name=fname, color='g',
                            dice=6, access='append')

            if self.evaluated:
                if self.parmsh['estimated']:
                    g.cloud(self.pe, name=fname, color='b',
                            dice=6, access='append')

                if self.parmsh['estimated_LS']:
                    g.cloud(self.p_LS, name=fname, color='r',
                            dice=6, access='append')

                if self.parmsh['grid']:
                    if self.evaluated:
                        g.cloud(self.g.p, name=fname, color='k', dice=2,
                                R=0.1, access='append')

            if self.parmsh['display']:
                print filename
                chaine = "geomview  -nopanel  -b 1 1 1 " + \
                    filename + " 2>/dev/null &"
                os.system(chaine)
            else:
                return(fname)
        else:
            print 'constraint is not runnable. It can not be displayed'
