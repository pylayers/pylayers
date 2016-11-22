from __future__ import print_function
import numpy as np
from pylayers.location.algebraic.algebraic import Algloc
from pylayers.location.observables import Observables

method = ['ls', 'wls', 'ml']

an = np.array([[0, 1, 2.], [0, 3, 1], [2, 1, 3],
               [2, 3, -1], [1, 0, 5], [1, 4, 0]])
an = an.T

bn = np.array([1, 1, 2.])

O_toa = Observables(an=an, bn=bn, mode='toa')

A_toa = Algloc(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
               toa_std=O_toa.noise_model['std'], bnGT=bn)

print('Blind node :' + str(bn) + '\n')

print('---')
print('TOA')
print('---')
for m in method:
    print(m)
    pe = A_toa.locate(method=m)
    print(pe)
print('crb')
print(A_toa.crb(A_toa.bnGT))


O_rss = Observables(an=an, bn=bn, mode='rss')
A_rss = Algloc(an_rss=O_rss.an, rss=O_rss.rp, rss_std=O_rss.noise_model['std'],
               rss_np=2., PL0=40.04, d0=1., bnGT=bn)

print('---')
print('RSS')
print('---')
for m in method:
    print(m)
    pe = A_rss.locate(method=m)
    print(pe)
print('crb')
print(A_rss.crb(A_rss.bnGT))


O_tdoa = Observables(an=an, bn=bn, mode='tdoa')
A_tdoa = Algloc(an_tdoa=O_tdoa.an, tdoa=O_tdoa.drng,
                tdoa_ref=O_tdoa.an_ref, tdoa_std=0.05, bnGT=bn)

print('---')
print('TDOA')
print('---')
for m in method:
    print(m)
    pe = A_tdoa.locate(method=m)
    print(pe)
print('crb')
print(A_tdoa.crb(A_tdoa.bnGT))

A_toa_tdoa = A_toa + A_tdoa

print('---')
print('TOA + TDOA')
print('---')
for m in method:
    print(m)
    pe = A_toa_tdoa.locate(method=m)
    print(pe)
print('crb')
print(A_toa_tdoa.crb(A_toa_tdoa.bnGT))

A_toa_rss = A_toa + A_rss

print('---')
print('TOA + RSS')
print('---')
for m in method:
    print(m)
    pe = A_toa_rss.locate(method=m)
    print(pe)
print('crb')
print(A_toa_rss.crb(A_toa_rss.bnGT))


A_tdoa_rss = A_tdoa + A_rss

print('---')
print('TDOA + RSS')
print('---')
for m in method:
    print(m)
    pe = A_tdoa_rss.locate(method=m)
    print(pe)
print('crb')
print(A_tdoa_rss.crb(A_tdoa_rss.bnGT))


A_full = A_tdoa + A_rss + A_toa

print('---')
print('full')
print('---')
for m in method:
    print(m)
    pe = A_full.locate(method=m)
    print(pe)
print('crb')
print(A_full.crb(A_full.bnGT))
