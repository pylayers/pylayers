import numpy as np
import matplotlib.pyplot as plt
from pylayers.signal.bsignal import *
x = np.arange(2,8,0.1)
y = np.ones(len(x))
U = Usignal(x,y)
fig,ax = U.plot()
U.window('hamming')
fig,ax = U.plot()
