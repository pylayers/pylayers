from pylayers.antprop.antenna import *
import matplotlib.pyplot as plt

#
# Be careful for trx fil you have to specify the number of column defant.trx
# has 7 column
#

A = Antenna('trx','defant.trx',directory='ant',nf=121,ntheta=37,nphi=72)
A.vshd()
A.C.show(typ='s1')
#plt.show()
#A.C.s3tos2()
