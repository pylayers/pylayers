import matplotlib.pyplot as plt
from pylayers.antprop.antenna import *
A = Antenna()
fig,ax = A.plotG(fGHz=[2,3,4],plan='theta',angdeg=0)
