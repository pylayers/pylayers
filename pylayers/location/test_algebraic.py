from __future__ import print_function
import numpy as np
from pylayers.location.algebraic.algebraic import Algloc
from pylayers.location.observables import Observables


an = np.array([[0, 1], [0, 3], [2, 1], [2, 3], [1, 0]])
an = an.T

bn = np.array([1,2])

O = Observables(an=an,bn=bn,mode='toa')

A = Algloc(an_toa=O.an, toa=O.rng + O.noise, toa_std=O.noise_model['std'])

print('Blind node :' + str(bn) + '\n')
print ('TOA')
print ('---\n')
print ('LS')
print(A.ls_locate())
print ('WLS')
print(A.wls_locate())
print ('ML')
print(A.ml_locate())


O = Observables(an=an,bn=bn,mode='rss')

A = Algloc(an_rss=O.an, rss=O.rp , rss_std=O.noise_model['std'],
    rss_np=2.,PL0=40.04,d0=1.)
print ('RSS')
print ('---\n')
print ('LS')
print(A.ls_locate())
print ('WLS')
print(A.wls_locate())
print ('ML')
print(A.ml_locate())