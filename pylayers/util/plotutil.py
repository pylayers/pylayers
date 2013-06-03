# -*- coding:Utf-8 -*-
import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import doctest
import pdb


def displot(pt, ph,color='black',fig=None,ax =None,linewidth=2):
    """ discontinuous plot

    Parameters
    ----------
    pt:
        tail points array (2 x (2*Nseg))
    ph :
        head points array (2 x (2*Nseg))
    col : string 
        color name

    Returns
    -------
    f,a
        fig and ax

    Examples
    --------

    .. plot::
        :include-source:

        >>> import scipy as sp
        >>> import matplotlib.pyplot as plt
        >>> from pylayers.util.geomutil import *
        >>> N   = 10
        >>> pt  = sp.rand(2,N)
        >>> ph  = sp.rand(2,N)
        >>> f,a = displot(pt,ph)
        >>> txt = plt.title('pylayers.util.geomutil.displot(pt,ph) : plot 10 random segments')

    """
    if fig == None:
        fig = plt.gcf()
        ax  = fig.gca()
    Nseg = np.shape(pt)[1]
    pz = np.empty((2,))
    pn = np.zeros((2,))

    for i in range(Nseg):
        pz = np.vstack((pz, pt[:, i], ph[:, i], pn))

    m1 = np.array([0, 0, 1])
    mask = np.kron(np.ones((2, Nseg)), m1)
    pzz = pz[1:, :].T
    vertices = np.ma.masked_array(pzz, mask)
    ax.plot(vertices[0, :], vertices[1, :], color=color,linewidth=linewidth)
    return fig, ax

def pol3D(fig,rho,theta,phi,sf=False,shade=True,title='pol3D'):
    """ polar 3D  surface plot

    Parameters
    ----------
    rho  : np.array
          t  x p
    theta : np.array
          1 x t
    phi  : np.array
          1 x p

    Example
    -------

    >>> from pylayers.util.plotutil import *
    >>> import numpy as np
    >>> theta = np.linspace(0,np.pi,90)
    >>> phi = np.linspace(0,2*np.pi,180)
    >>> rho = np.ones((len(theta),len(phi)))
    >>> fig=plt.figure()
    >>> pol3D(fig,rho,theta,phi)
    >>> plt.show()

    """
    ax = axes3d.Axes3D(fig)

    # x : r x t x p
    x = rho * np.sin(theta[:,np.newaxis]) *  np.cos(phi[np.newaxis,:])
    y = rho * np.sin(theta[:,np.newaxis]) *  np.sin(phi[np.newaxis,:])
    z = rho * np.cos(theta[:,np.newaxis])

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.plot_surface(x, y, z, rstride=1, cstride=1,
                    cmap=cm.hot_r,linewidth=0,antialiased=False)
    plt.title(title)
    ax.set_xlim3d([-1, 1])
    ax.set_ylim3d([-1, 1])
    ax.set_zlim3d([-1, 1])

    if sf:
        sz = fig.get_size_inches()
        fig.set_size_inches(sz * 1.8)
        figname = figdirV + 'V' + str(n) + str(m)
        fig.savefig(figname + ext1, orientation='portrait')
        fig.savefig(figname + ext2, orientation='portrait')
        fig.savefig(figname + ext3, orientation='portrait')

def cylinder(fig,pa,pb,R):
    """ plot a cylinder
    pa : 3 x nc 
    pb : 3 x nc 
    R  : 1 x Nc

    >>> from pylayers.util.plotutil import *
    >>> import numpy as np
    >>> pa = np.array([0,0,0])
    >>> pb = np.array([0,0,10])
    >>> fig = plt.figure()
    >>> cylinder(fig,pa,pb,3)
    >>> plt.show()

    """
    try:
        nc = np.shape(pa)[1]
    except:
        nc = 1
        pa = pa.reshape(3,1)
        pb = pb.reshape(3,1)
    ax = fig.gca(projection='3d') 
    theta = np.linspace(0, 2 * np.pi, 40)
    # 3 x nc
    v = (pb-pa)
    mv = np.sqrt(np.sum(v*v,axis=0)).reshape(1,nc)
    vn = v/mv
    if nc>1:
        pr = sp.rand(3,nc)
    else:
        pr = sp.rand(3,1)
    pp = pr - np.sum(pr*vn,axis=0)
    mpp = np.sum(pp*pp,axis=0)
    pn = pp/mpp
    qn = np.cross(pn,vn,axisa=0,axisb=0,axisc=0)
    # 3 x nc x ntheta
    p = pa[:,:,np.newaxis] + \
            R*(np.cos(theta[np.newaxis,np.newaxis,:])*pn[:,:,np.newaxis] + \
               np.sin(theta[np.newaxis,np.newaxis,:])*qn[:,:,np.newaxis])
    ax.plot(p[0,0,:], p[1,0,:], p[2,0,:], label='parametric curve',color='b')
    p = pb[:,:,np.newaxis] + \
            R*(np.cos(theta[np.newaxis,np.newaxis,:])*pn[:,:,np.newaxis] + \
               np.sin(theta[np.newaxis,np.newaxis,:])*qn[:,:,np.newaxis])
    ax.plot(p[0,0,:], p[1,0,:], p[2,0,:], label='parametric curve',color='b')
    t = np.arange(0,1,0.05)
    p = pa[:,:,np.newaxis]+t[np.newaxis,np.newaxis,:]*(pb[:,:,np.newaxis]-pa[:,:,np.newaxis])
    ax.plot(p[0,0,:], p[1,0,:], p[2,0,:], label='parametric curve',color='b')

if (__name__ == "__main__"):
    doctest.testmod()
