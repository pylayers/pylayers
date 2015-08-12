from pylayers.antprop.antenna import *
import numpy as np
import matplotlib.pylab as plt
A = Antenna('defant.vsh3')
F = A.eval(grid=True)
