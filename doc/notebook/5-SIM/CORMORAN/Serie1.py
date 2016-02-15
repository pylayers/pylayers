from pylayers.measures.cormoran import *
import seaborn as sns
sns.set_style("white")
import seaborn as sns
sns.set_style("white")
# load serie 1
S = CorSer(serie=10,day=12)
AP = S.din.keys()
OB = S.B['Bernard'].dev.keys()
values = S.getlink(AP[0],OB[0])
plt.plot(values)
plt.show()
