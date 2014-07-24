import numpy as np
import pdb

def eqt(lL,lL0):
    """ encode lon Lat in quad tree integer

    Parameters
    ----------

    lL : nd.array (2xN)
        longitude Latitude
    lL0 : nd.array (,2)
        lower left corner of the 1degree tile

    """
    N = lL.shape[1]
    # offset from the lower left corner
    d = lL-lL0[:,None]
    dui8 = np.floor(d*256).astype('uint8')
    ndui8 = np.unpackbits(dui8).reshape(2,N,8)
    d16 = np.empty((N,16)).astype('int')
    d16[:,1::2]=ndui8[0,:,:]
    d16[:,0::2]=ndui8[1,:,:]
    ud8 = np.packbits(d16,axis=1)
    ud16 = (ud8[:,0]*256+ud8[:,1]).astype('uint16')
    return(ud8,ud16)

def dqt(ud16,lL0):
    """ decode quad tree integer to lon Lat

    Parameters
    ----------

    lL : nd.array (2xN)
        longitude Latitude
    lL0 : nd.array (,2)
        lower left corner of the 1degree tile

    """
    N = len(ud16)
    # offset from the lower left corner
    #d = lL-lL0[:,None]
    #dui8 = np.floor(d*256).astype('uint8')
    uh8  = ud16/256
    ul8  = ud16-uh8*256
    ud8 =  (np.vstack((uh8,ul8)).T).astype('uint8')
    ud16 = np.unpackbits(ud8).reshape(N,16)
    ndu8 = np.empty((2,N,8)).astype('int')
    ndu8[0,:,:]=ud16[:,1::2]
    ndu8[1,:,:]=ud16[:,0::2]
    du8 = np.packbits(ndu8).reshape(2,N)/256.
    lL = lL0[:,None]+du8

    return(lL)

lL0 = np.array([-2,48])
lL = lL0[:,None] + np.random.rand(2,10)

ud8,ud16 = eqt(lL,lL0)
lLb  = dqt(ud16,lL0)
