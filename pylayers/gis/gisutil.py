import numpy as np
import pdb

def ent(lL,lL0):
    """ encode lon Lat in natural integer

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
    lab = map(lambda x:'i'+str(x[0])+'-'+str(x[1]),dui8.T)
    return(lab)

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
    return(ud16 )

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

def ext2qt(extent=np.array([-1.8,-1.7,48.4,48.5]),lL0=np.array([-2,48])):
    """ convert an extent region into a list of qt regions
    """
    lm = extent[0]
    lM = extent[1]
    Lm = extent[2]
    LM = extent[3]
    lL = np.array([[lm,Lm],[lM,LM]]).T
    uf8 = ent(lL,lL0)
    ill = uf8[0].replace('i','').split('-')
    iur = uf8[1].replace('i','').split('-')
    il  = np.arange(eval(ill[0]),eval(iur[0]))
    iL  = np.arange(eval(ill[1]),eval(iur[1]))
    ltile = []
    for l in il:
        for L in iL:
            ltile.append('i'+str(l)+'-'+str(L))

    return(ltile)

def ctrad2qt(extent):
    """ convert center,radius into a list of qt regions
    """
    pass

if __name__=='__main__':
    lL0 = np.array([-2,48])
    lL = lL0[:,None] + np.random.rand(2,5)
    print lL
    ud16 = eqt(lL,lL0)
    un8  = ent(lL,lL0)
    lLb  = dqt(ud16,lL0)
