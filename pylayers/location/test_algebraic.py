import numpy as np
from pylayers.location.algebraic.algebraic import Algloc
from pylayers.location.observables import Observables


an = np.array([[0, 1], [0, 3], [2, 1], [2, 3], [1, 0]])
an = an.T

bn = np.array([2,3])

O = Observables(an=an,bn=bn,mode='toa')

A = Algloc(an_toa=O.an, toa=O.rng + O.noise, toa_std=O.noise_model['std'])

print A.ls_locate()