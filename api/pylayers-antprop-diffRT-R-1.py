import numpy as np
th = np.linspace(0,np.pi/2,180)[None,:]
fGHz = 0.3
lamda = 0.3/fGHz
k = np.array([2*np.pi/2])[:,None]
Rs,Rh = R(th,k,9,0,0.01,1,0,0)
