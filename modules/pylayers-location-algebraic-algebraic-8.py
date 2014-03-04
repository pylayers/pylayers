nRN = 4
dim = 3 # 2 for 2D, 3 for 3D
L = 20.
c = 0.3
BN = L*sp.rand(dim,1)
BN0 = L*sp.rand(dim,1)
RN_TOA = L*sp.rand(dim,nRN)
RN_RSS = L*sp.rand(dim,nRN)
RN_TDOA = L*sp.rand(dim,nRN)

d_TOA = dist(RN_TOA,BN,0) # actual distances
TOF = d_TOA/c # actual TOA
TOA_std = 0.001/c*np.ones(np.shape(TOF))
TOA = TOF + TOA_std

rss_std = 0.001 * np.ones(nRN)
rss_np = 2.645 * np.ones(nRN)
pl0 = 34.7*np.ones(nRN)
d0 = 1.
d_RSS = dist(RN_RSS,BN,0) # actual distances
X = rss_std * np.random.randn(np.shape(pl0)[0])
rss_db = pl0-10*rss_np*np.log10(d_RSS/d0)+X

RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
d = dist(RN_TDOA,BN,0)
dr = dist(RNr_TDOA,BN,0)
TDOF = (d-dr)/c # actual TDOA
TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
TDOA = TDOF + TDOA_std

nodes={}
nodes['BN']= BN
nodes['RN_RSS']= RN_RSS
nodes['RN_TOA']= RN_TOA
nodes['RN_TDOA']= RN_TDOA
nodes['RNr_TDOA']= RNr_TDOA

ldp={}
ldp['RSS'] = rss_db
ldp['RSS_std'] = rss_std
ldp['RSS_np'] = rss_np
ldp['d0'] = d0
ldp['PL0'] = pl0
ldp['TOA'] = TOA
ldp['TOA_std'] = TOA_std
ldp['TDOA'] = TDOA
ldp['TDOA_std'] = TDOA_std

print 'Nodes'
print nodes
print 'LDPs'
print ldp
print 'BN0:initial guess for ML estimator'
print BN0
