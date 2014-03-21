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

    def __init__(self, devname='Telephone1'):
        """ init
        """
        self.name = devname
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
        for k in self.std:
            if self.std[k]['on']:
                power = 'on'
            else :
                power ='off'
            ant = str([a for a in self.std[k]['ant']])
            s = s + str(k) + '\n'
            s = s + '-'*len(k) + '\n'
            s = s + '{0:5} | {1:7} |{2:10} | {3:10} | {4:10} '.format('power', 'channel', 'modulation', 'code rate', 'antenna(s)') + '\n'
            s = s +'{0:5} | {1:7} |{2:10} | {3:10} | {4:10} '.format(power, self.std[k]['chan'], self.std[k]['mod'], self.std[k]['cr'], ant) + '\n\n'



        # s = s + '{0:7} | {1:20} |{2:20} '.format('on/off', 'standard', 'antenna(s)') + '\n'
        
        # for k in self.std:
        #     if self.std[k]['on']:
        #         power = 'on'
        #     else :
        #         power ='off'
        #     ant = str([a for a in self.std[k]['ant']])
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
        stds = json.load(fp)
        fp.close()


        dim = dev[devname]['dimensions']
        ant = dev[devname]['antennas']
        std = dev[devname]['standards']
        # meter conversion
        self.dim = np.array((dim['height'], dim['width'], dim['depth'])) / 1000
        self.ant = {}
        for k in ant.keys():
            self.ant[k] = {}
            self.ant[k]['p'] = np.array(eval(ant[k]['p']))
            self.ant[k]['T'] = np.array(eval(ant[k]['T']))
        self.std = {}
        for k in std.keys():
            self.std[k] = {}
            W = Wstandard(k)
            self.std[k]['on']=True
            self.std[k]['ant'] = std[k]['antenna']
            self.std[k]['chan'] = W.chan.keys()[0]
            self.std[k]['epwr'] = W.power(self.std[k]['chan'])
            self.std[k]['sens'] = [0]*len(self.std[k]['ant'])
            self.std[k]['cr'] = W['crate'][0]
            self.std[k]['mod'] = W['modulation'][0]
            self.std[k]['cca'] = ''


