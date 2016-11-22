from __future__ import print_function
import numpy as np
from pylayers.location.algebraic.algebraic import Algloc
from pylayers.location.observables import Observables


an = np.array([[0, 1, 2.], [0, 3, 1], [2, 1, 3],
               [2, 3, -1], [1, 0, 5], [1, 4, 0]])
an = an.T

bn = np.array([1, 1, 2.])

O_toa = Observables(an=an, bn=bn, mode='toa')

A_toa = Algloc(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
               toa_std=O_toa.noise_model['std'], bnGT=bn)

print('Blind node :' + str(bn) + '\n')
print('TOA')
print('---')
print('LS')
print(A_toa.ls_locate())
print('WLS')
print(A_toa.wls_locate())
print('ML')
print(A_toa.ml_locate())
print('CRB')
print(A_toa.crb(A_toa.bnGT))

O_rss = Observables(an=an, bn=bn, mode='rss')

A_rss = Algloc(an_rss=O_rss.an, rss=O_rss.rp, rss_std=O_rss.noise_model['std'],
               rss_np=2., PL0=40.04, d0=1., bnGT=bn)
print('RSS')
print('---')
print('LS')
print(A_rss.ls_locate())
print('WLS')
print(A_rss.wls_locate())
print('ML')
print(A_rss.ml_locate())
print('CRB')
print(A_rss.crb(A_toa.bnGT))

O_tdoa = Observables(an=an, bn=bn, mode='tdoa')


A_tdoa = Algloc(an_tdoa=O_tdoa.an, tdoa=O_tdoa.drng + O_tdoa.noise,
                tdoa_ref=O_tdoa.an_ref, tdoa_std=O_tdoa.noise_model['std'], bnGT=bn)
print('TDOA')
print('---')
print('LS')
print(A_tdoa.ls_locate())
print('WLS')
print(A_tdoa.wls_locate())
print('ML')
print(A_tdoa.ml_locate())
print('CRB')
print(A_tdoa.crb(A_tdoa.bnGT))


A_toa_tdoa = A_toa + A_tdoa
print('TOA + TDOA')
print('---')
print('LS')
print(A_toa_tdoa.ls_locate())
print('WLS')
print(A_toa_tdoa.wls_locate())
print('ML')
print(A_toa_tdoa.ml_locate())
print('CRB')
print(A_toa_tdoa.crb(A_toa_tdoa.bnGT))

A_toa_rss = A_toa + A_rss
print('TOA + RSS')
print('---')
print('LS')
print(A_toa_rss.ls_locate())
print('WLS')
print(A_toa_rss.wls_locate())
print('ML')
print(A_toa_rss.ml_locate())
print('CRB')
print(A_toa_rss.crb(A_toa_rss.bnGT))


A_tdoa_rss = A_tdoa + A_rss

print('TDOA + RSS')
print('---')
print('LS')
print(A_tdoa_rss.ls_locate())
print('WLS')
print(A_tdoa_rss.wls_locate())
print('ML')
print(A_tdoa_rss.ml_locate())
print('CRB')
print(A_tdoa_rss.crb(A_tdoa_rss.bnGT))


A_full = A_tdoa + A_rss + A_toa

print('full')
print('---')
print('LS')
print(A_full.ls_locate())
print('WLS')
print(A_full.wls_locate())
print('ML')
print(A_full.ml_locate())
print('CRB')
print(A_full.crb(A_full.bnGT))


# A_toa_tdoa = Algloc(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
#                     toa_std=O_toa.noise_model[
#                         'std'], an_tdoa=O_tdoa.an, tdoa=O_tdoa.drng + O_tdoa.noise,
#                     tdoa_ref=O_tdoa.an_ref, tdoa_std=O_tdoa.noise_model['std'], bnGT=bn)
# print('TOA + TDOA')
# print('---')
# print('LS')
# print(A_toa_tdoa.ls_locate())
# print('WLS')
# print(A_toa_tdoa.wls_locate())
# print('ML')
# print(A_toa_tdoa.ml_locate())


# A_toa_rss = Algloc(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
#                    toa_std=O_toa.noise_model[
#                        'std'], an_rss=O_rss.an, rss=O_rss.rp, rss_std=O_rss.noise_model['std'],
#                    rss_np=2., PL0=40.04, d0=1., bnGT=bn)

# print('TOA + RSS')
# print('---')
# print('LS')
# print(A_toa_rss.ls_locate())
# print('WLS')
# print(A_toa_rss.wls_locate())
# print('ML')
# print(A_toa_rss.ml_locate())

# A_tdoa_rss = Algloc(an_tdoa=O_tdoa.an, tdoa=O_tdoa.drng + O_tdoa.noise,
#                     tdoa_ref=O_tdoa.an_ref, tdoa_std=O_tdoa.noise_model[
#                         'std'], an_rss=O_rss.an, rss=O_rss.rp, rss_std=O_rss.noise_model['std'],
#                     rss_np=2., PL0=40.04, d0=1., bnGT=bn)

# print('TDOA + RSS')
# print('---')
# print('LS')
# print(A_tdoa_rss.ls_locate())
# print('WLS')
# print(A_tdoa_rss.wls_locate())
# print('ML')
# print(A_tdoa_rss.ml_locate())


# A_full = Algloc(an_tdoa=O_tdoa.an, tdoa=O_tdoa.drng + O_tdoa.noise,
#                 tdoa_ref=O_tdoa.an_ref, tdoa_std=O_tdoa.noise_model[
#                     'std'], an_rss=O_rss.an, rss=O_rss.rp, rss_std=O_rss.noise_model['std'],
#                 rss_np=2., PL0=40.04, d0=1., an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
#                 toa_std=O_toa.noise_model['std'], bnGT=bn)

# print('TDOA + RSS')
# print('---')
# print('LS')
# print(A_full.ls_locate())
# print('WLS')
# print(A_full.wls_locate())
# print('ML')
# print(A_full.ml_locate())
