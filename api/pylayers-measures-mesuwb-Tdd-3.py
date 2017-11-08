from pylayers.util.project import *
from pylayers.measures.mesuwb import *
M  = UWBMeasure(2)
T  = M.tdd
s1 = T.show_span()
plt.show()
