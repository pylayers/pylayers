from pylayers.simul.simulem import *
from matplotlib.pylab import *
S = Simul()
S.load('where2.ini')
vc = S.VC(1,1)
figure(figsize=(16,8))
vc.doadod()
show()