import numpy as np
import matplotlib.pyplot as plt
from pylayers.signal.bsignal import *
x = np.linspace(-10,10,100)
y = np.sin(2*np.pi*12*x)+np.random.normal(0,0.1,len(x))
s = TUsignal(x,y)
fig,ax = s.plot(types=['v'])
txt1 = plt.title('before gating')
plt.show()
s.gating(-3,4)
fig,ax=s.plot(types=['v'])
txt2 = plt.title('after gating')
plt.show()
