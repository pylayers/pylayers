#!/usr/bin/python
# -*- coding: utf-8 -*-
#
"""
.. currentmodule:: pylayers.signal.device


This module describes the radio devices to be used for electromagentic simulations

Device Class
============

.. autosummary::
    :toctree: generated/

    Device.__init__

"""

import doctest
import numpy as np
import matplotlib.pylab as plt
import json


from pylayers.util.project import *
import pylayers.util.pyutil as pyu
# Handle Antenna
from pylayers.antprop.antenna import Antenna
# Handle Standards
from pylayers.signal.standard import Wstandard
try:
    from tvtk.api import tvtk
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi import mlab
except:
    print 'Layout:Mayavi is not installed'
import pdb


class Device(object):

    """ Device Class

    """

    def __init__(self, devname='Telephone1', ID =0, owner='',typ='ap'):
        """ init of Device class

        Attributes
        ----------

        owner : string
            owner of the device
        name : string
            name of the device
        ID : int 
            device id
        ant : dict
            dictionnary of used antennas with:
            key : 
                antenna name
            values : 
                p : relative position of antenna on the device
                T : Rotaion matrice of antenna (acs)
        dim : ndarray
            dimension of the device
        typ : string ('ap'|....)
            type of the device
        wstd : dict
            dictionnary of wireless standards
            key : 
                wireless standard name
            values : dict of std information (from Wstandard class)

        See Also
        --------

        pylayers.signal.standard


        """
        self.owner = owner
        self.name = devname
        self.ID = ID
        self.typ = typ
        self.load(devname)

    def __repr__(self):

        s = self.name + '\n'
        s = s + ''*len(self.name) + '\n\n'
        s = s + 'Dimensions' + '\n'
        s = s + '==========' + '\n'
        s = s + str(self.dim[0]*1000) + ' (mm) x ' +\
                str(self.dim[1]*1000) + ' (mm) x ' +\
                str(self.dim[2]*1000) + ' (mm)\n\n'

        s = s + 'Wireless Standards' + '\n'
        s = s + '==================' + '\n'
        for k in self.wstd:
            if self.wstd[k]['power']:
                power = 'on'
            else :
                power ='off'
            ant = str([a for a in self.wstd[k]['ant']])
            s = s + str(k) + '\n'
            s = s + '-'*len(k) + '\n'
            s = s + '{0:5} | {1:7} |{2:10} | {3:10} | {4:10} '.format('power', 'channel', 'modulation', 'code rate', 'antenna(s)') + '\n'
            s = s +'{0:5} | {1:7} |{2:10} | {3:10} | {4:10} '.format(power, self.wstd[k]['chan'], self.wstd[k]['mod'], self.wstd[k]['cr'], ant) + '\n\n'



        # s = s + '{0:7} | {1:20} |{2:20} '.format('on/off', 'standard', 'antenna(s)') + '\n'
        # for k in self.std:
        #     if self.std[k]['on']:
        # for k in self.wstd:
        #     if self.wstd[k]['on']:
        #         power = 'on'
        #     else :
        #         power ='off'
        #     ant = str([a for a in self.wstd[k]['ant']])
        #     s = s +'{0:7} | {1:20} |{2:20} '.format(power, str(k), ant ) + '\n'
        s = s + '\n\nAntennas' + '\n'
        s = s + '========' + '\n'
        for k in self.ant:
            s = s + str(k) + '\n'
            s = s + '-'*len(k) + '\n'
            s = s + 'On device Relative position: \n' + str(self.ant[k]['p']) + '\n'
            s = s + 'On device Rotation Matrice: \n' + str(self.ant[k]['T']) + '\n\n'

        return s

    def load(self, devname):
        """ load a device

        Parameters
        ----------

        devname : string
            Name of the device to be loaded

        """

        fp = open(pyu.getlong('devices.json', pstruc['DIRSIMUL']))
        dev = json.load(fp)
        fp.close()

        fp = open(pyu.getlong('wstd.json',pstruc['DIRSIMUL']))
        wstds = json.load(fp)
        fp.close()


        dim = dev[devname]['dimensions']
        ant = dev[devname]['antennas']
        wstd = dev[devname]['standards']
        # meter conversion
        self.dim = np.array((dim['height'], dim['width'], dim['depth'])) / 1000
        self.ant = {}

        for k in ant.keys():
            self.ant[k] = {}
            self.ant[k]['p'] = np.array(eval(ant[k]['p']))
            self.ant[k]['T'] = np.array(eval(ant[k]['T']))
        self.wstd = {}
        for k in wstd.keys():
            self.wstd[k] = {}
            W = Wstandard(k)
            self.wstd[k]['power']=True
            self.wstd[k]['ant'] = wstd[k]['antenna']
            self.wstd[k]['chan'] = W.chan.keys()[0]
            self.wstd[k]['epwr'] = W.power(self.wstd[k]['chan'])
            self.wstd[k]['sens'] = [0]*len(self.wstd[k]['ant'])
            self.wstd[k]['cr'] = W['crate'][0]
            self.wstd[k]['mod'] = W['modulation'][0]
            self.wstd[k]['cca'] = ''


