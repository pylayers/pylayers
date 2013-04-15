# -*- coding:Utf-8 -*-
import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import doctest


def pol3D(fig,rho,theta,phi,sf=False):
    """ polar 3D surface plot

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
    >>> theta = np.linspace(0,np.pi,30)
    >>> phi = np.linspace(0,2*np.pi,60)
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

    ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.hot_r)

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

if (__name__ == "__main__"):
    doctest.testmod()
