# -*- coding: utf-8 -*-
#
#
#   Building  Module
#
#
#
r"""

"""
from __future__ import print_function
try:
    from tvtk.api import tvtk
    from mayavi import mlab
except:
    print('Layout:Mayavi is not installed')
import pdb
import sys
import os
import logging
import copy
import glob
import time
import tqdm
import numpy as np
import scipy as sp
import doctest
import matplotlib.pyplot as plt
import pylayers.util.project as pro



class Building(pro.PyLayers,list):

    def __init__(self):

        self.lzfloor = []
        self.lzceil = []
        self.lfilename = []
        self.Nfloor = 0


    def __repr__(self):

        s = 'Building has ' + str(self.Nfloor) + ' floor(s)'
        s = s + '\n'
        s = s + '\n'
        for f in range(self.Nfloor):
            s = s + 'floor ' + str(f).zfill(2) + ' : ' +str(self.lfilename[f])
            s = s + '\n'
            s = s + '--------'
            s = s + '\n'
            s = s + 'height floor :' + str(self.lzfloor[f])
            s = s + '\n'
            s = s + 'height ceil :' + str(self.lzceil[f])
            s = s + '\n'
            s = s + '\n'
        return s

    def __add__(self,lL):
        if not isinstance(lL,list):
            lL = [lL]
        return self.append(lL)

    def pop(self,n):

        Lpop = super(Building, self).pop(n)
        zfpop = self.lzfloor.pop(n)
        zcpop = self.lzceil.pop(n)
        self.lfilename.pop(n)
        self.Nfloor -= 1
        for n in range(n,len(self)):
            self.lzfloor[n] = self.lzfloor[n] - \
                              (zcpop-zfpop) -\
                              sum(Lpop.sl['FLOOR']['lthick']) -\
                              sum(Lpop.sl['CEIL']['lthick'])
            self.lzceil[n] = self.lzceil[n] - \
                              (zcpop-zfpop) -\
                              sum(Lpop.sl['FLOOR']['lthick']) -\
                              sum(Lpop.sl['CEIL']['lthick'])

    def append(self,lL):
        if not isinstance(lL,list):
            lL = [lL]
        for ul, L in enumerate(lL):
            if L.typ != 'indoor':
                raise AttributeError('Layout must be indoor for creating a Building')
            super(Building, self).append(L)
            if ul == 0 and len(self) == 1:
                self.lzfloor=[L.zfloor]
                self.lzceil=[L.zceil]
            else:
                self.lzfloor.append(L.zfloor + \
                                   self.lzceil[-1] + \
                                   sum(self[-1].sl['CEIL']['lthick'])+
                                   sum(L.sl['FLOOR']['lthick'])
                                   )
                self.lzceil.append(L.zceil + \
                                   self.lzceil[-1] + \
                                   sum(self[-1].sl['CEIL']['lthick'])+
                                   sum(L.sl['FLOOR']['lthick'])
                                   )
            self.lfilename.append(L._filename)
            self.Nfloor += 1









