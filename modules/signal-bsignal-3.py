import numpy as np
import matplotlib.pyplot as plt
from pylayers.signal.bsignal import *
x = np.arange(2,8,0.1)
y = np.ones(len(x))
U = FUsignal(x,y)
fi = plt.figure()
U.plot()
U.window('hamming')
U.plot()
