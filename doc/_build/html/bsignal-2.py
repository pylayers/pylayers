from pylayers.simul.simulem import *
from matplotlib.pylab import *
S = Simul()
S.load('where2.ini')
st = S.wav.st
sf = S.wav.sf
st.plot()
figure()
sf.plot()
show()