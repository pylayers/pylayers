f,a = H.show(cmap='jet')

# We keep only the significant rays

H.cut()
f,a = H.show(cmap='jet')

# We then transform into the time domain

T1 = H.totime()
T1.plot()

H.minphas()
T2 = H.totime()
T2.plot()
