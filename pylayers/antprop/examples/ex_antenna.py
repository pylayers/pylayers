from pylayers.antprop.antenna import *
import matplotlib.pyplot as plt

A = Antenna('defant.vsh3')
# Create pattern 
A.Fsynth3()
A.show3()
#A2 = Antenna('mat',)
#A.vshd()
#A.C.show(typ='s1')
#plt.show()
#A.C.s3tos2()

