# -*- coding:Utf-8 -*-
from __future__ import print_function

import numpy as np
import doctest
import scipy as sp
from scipy import optimize
import numpy.linalg as nplg
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from pylayers.util.geomutil import dist
from pylayers.location.observables import Observables
import string
import pdb


class Loc(object):
    """
    """
    def __init__(self, O=Observables()):

        self.O = O
