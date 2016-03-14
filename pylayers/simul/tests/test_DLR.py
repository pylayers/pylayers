import mayavi.mlab as mlab
from pylayers.simul.simulem import * #Obsolete
from pylayers.simul.link import *
from pylayers.antprop.rays import *
from pylayers.antprop.channel import *
from pylayers.antprop.signature import *
import pylayers.util.pyutil as pyu
from pylayers.gis.layout import *
from pylayers.util.project import *
import pylayers.signal.bsignal as bs
from datetime import datetime
import time
import pdb
import pickle
import numpy as np
#import matplotlib.pyplot as plt

# Set the frequency range
fGHz = np.arange(2,6,0.01)
# Define the Layout
L = Layout('DLR2.ini')
# Create a link for the Layout
S = DLink(L=L,fGHz=fGHz)
ak,tauk=S.eval(force=['Ct','H'])
S._show3()
