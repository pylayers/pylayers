from pylayers.signal.bsignal import *
import matplotlib.pyplot as plt
si = Bsignal()
si.x= np.arange(100)
si.y= np.arange(100)
f,a = si.stem()
