# -*- coding:Utf-8 -*-
import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import doctest
import pdb

def mulcplot(x,y,**kwargs):
    """ handling multiple complex variable plots

    Parameters
    ----------

    x : ndarray  (Nc x Nx)
    y : ndarray  (Nv x Ny)

    types : 'm'   : modulus 
            'v'   : value
            'l10' : dB (10 log10)
            'l20' : dB (20 log10)
            'd'   : phase degrees
            'r'   : phase radians
            'du'  : phase degrees unwrap 
            'ru'  : phase radians unwrap 
            'gdn' : group delay (ns) 
            'gdm' : group distance (m) 
            're'  : real part
            'im'  : imaginary part 

    dB : bool
        False
    fig = []    
    ax  = []    
    nlg  : int 
        number of lines 

    Examples
    --------

    >>> from pylayers.util.plotutil import *
    >>> import numpy as np

    Notes
    -----

    If len(y.shape) > 2 the two first axes are used as nlin and ncol this
    takes the priority over the pased values nlin and ncol 


    """
    defaults = {'types':['l20'],
                'titles':[''],
                'labels':[''],
                'xlabels':['time (ns)'],
                'ylabels':['Amplitude (dB)'],
                'ncol':1,
                'nlin':1,
                'fig':[],
                'ax':[],
               }

    # radians to degree coefficient   
    rtd = 180./np.pi

    # smart placement of legend box
    plt.rcParams['legend.loc'] = 'best'
    
    grid = False
    if 'nlin' in kwargs:
        grid=True

    for key, value in defaults.items():
        if key not in kwargs:
            kwargs[key] = value
   
    if 'ylabels' not in kwargs:
        ylabels = []
        for t in kwargs['types']:
            if t=='m':
                ylabels.append('Amplitude'),
            if t=='v':
                ylabels.append('Amplitude'),
            if t=='l10':
                ylabels.append('Amplitude (dB)'),
            if t=='l20':
                ylabels.append('Amplitude (dB)'),
            if t=='d':
                ylabels.append('Phase (deg)'),
            if t=='r':
                ylabels.append('Phase (rad)'),
            if t=='du':
                ylabels.append('Unwrapped Phase (deg)'),
            if t=='ru':
                ylabels.append('Unwrapped Phase (rad)'),
            if t=='re':
                ylabels.append('Real part'),
            if t=='im':
                ylabels.append('Imaginary part'),
            if t=='gdn':
                ylabels.append('Group delay (ns)'),
            if t=='gdm':
                ylabels.append('Group distance (m)'),
    else:
        ylabels = kwargs['ylabels']

    fig = kwargs['fig']
    ax = kwargs['ax']
    nlin = kwargs['nlin']
    ncol = kwargs['ncol']
    types = kwargs['types']
    titles = kwargs['titles']
    labels = kwargs['labels']
    xlabels = kwargs['xlabels']

    ntypes = np.prod(np.array(types).shape)
    ntitles = np.prod(np.array(titles).shape)
    nlabels = np.prod(np.array(labels).shape)
    nxlabels = np.prod(np.array(xlabels).shape)
    nylabels = np.prod(np.array(ylabels).shape)


    # filtering kwargs argument for plot function 
    args ={}
    for k in kwargs:
        if k not in defaults.keys():
            args[k]=kwargs[k]

    shx = x.shape
    shy = y.shape

    if len(shy)>2:
        ydim = shy[0:-1]
        nlin = ydim[0]
        ncol = ydim[1]
    else:
        if not grid:
            if len(shy)>1:
                nlin = shy[0]
                ncol = 1
            else:
                nlin = 1
                ncol = 1
                y = y[np.newaxis,:]

    
    # below is the same axis constraint as for Bsignal object
    assert(shx[0]==shy[-1]), "plotutil : x,y shape incompatibility"

    nfigy = np.prod(np.array(y.shape[0:-1]))

    assert((nfigy==ncol*nlin) | (nfigy==1))
    assert((nlabels==nfigy)|(nlabels==1))
    assert((ntitles==ncol*nlin)|(ntitles==1))
    assert((nxlabels==nfigy)|(nxlabels==1))
    assert((nylabels==nfigy)|(nxlabels==1))

    if ax==[]:    
        fig,ax=plt.subplots(nlin,ncol,sharey=True,sharex=True)
        if (nlin==1)&(ncol==1):
            ax = np.array(ax)[np.newaxis,np.newaxis]
        else:    
            if nlin==1:
                ax = ax[np.newaxis,:]
            if ncol==1:
                ax = ax[:,np.newaxis]
   
    for l in range(nlin):
        for c in range(ncol):
            if len(shy)>2:
                if nlabels>1:
                    lablc = labels[l,c]
                else:
                    lablc = labels[0]
                if types[0]=='v':
                        ax[l,c].plot(x,y[l,c,:],label=lablc,**args)
                if types[0]=='r':
                    ax[l,c].plot(x,np.angle(y[l,c,:]),label=lablc,**args)
                if types[0]=='ru':
                    ax[l,c].plot(x,np.unwrap(np.angle(y[l,c,:])),label=lablc,**args)
                if types[0]=='d':
                    ax[l,c].plot(x,np.angle(y[l,c,:])*rtd,label=lablc,**args)
                if types[0]=='du':
                    ax[l,c].plot(x,np.unwrap(np.angle(y[l,c,:]))*rtd,label=lablc,**args)
                if types[0]=='m':
                    ax[l,c].plot(x,np.abs(y[l,c,:]),label=lablc,**args)
                if types[0]=='l10':
                    ax[l,c].plot(x,10*np.log10(np.abs(y[l,c,:])),label=lablc,**args)
                if types[0]=='l20':
                    ax[l,c].plot(x,20*np.log10(np.abs(y[l,c,:])),label=lablc,**args)
                if types[0]=='re':
                    ax[l,c].plot(x,np.real(y[l,c,:]),label=lablc,**args)
                if types[0]=='im':
                    ax[l,c].plot(x,np.imag(y[l,c,:]),label=lablc,**args)
                if types[0]=='gdn':
                    df  = x[1]-x[0]
                    ax[l,c].plot(x[0:-1],-np.diff(np.unwrap(np.angle(y[l,c,:])))/(2*np.pi*df),label=lablc,**args)
                if types[0]=='gdm':
                    df  = x[1]-x[0]
                    ax[l,c].plot(x[0:-1],-0.3*np.diff(np.unwrap(np.angle(y[l,c,:])))/(2*np.pi*df),label=lablc,**args)
                
                if nxlabels>1:
                    ax[l,c].set_xlabel(xlabels[l,c])
                else:
                    ax[l,c].set_xlabel(xlabels[0])
                if nylabels>1:
                    ax[l,c].set_ylabel(ylabels[l,c])
                else:
                    ax[l,c].set_ylabel(ylabels[0])
                if ntitles>1:    
                    ax[l,c].set_title(titles[l,c])
                else:
                    ax[l,c].set_title(titles[0])
                ax[l,c].legend()
                   
            else:
                k = l*ncol+c
                if types[k%ntypes]=='v':
                    #ax[l,c].plot(x[k%nfigx,:],y[k%nfigy,:],label=labels[k%nlabels],**args)
                    ax[l,c].plot(x,y[k%nfigy,:],label=labels[k%nlabels],**args)
                if types[k%ntypes]=='r':
                    ax[l,c].plot(x,np.angle(y[k%nfigy,:]),label=labels[k%nlabels],**args)
                if types[k%ntypes]=='ru':
                    ax[l,c].plot(x,np.unwrap(np.angle(y[k%nfigy,:])),label=labels[k%nlabels],**args)
                if types[k%ntypes]=='d':
                    ax[l,c].plot(x,np.angle(y[k%nfigy,:])*rtd,label=labels[k%nlabels],**args)
                if types[k%ntypes]=='du':
                    ax[l,c].plot(x,np.unwrap(np.angle(y[k%nfigy,:]))*rtd,label=labels[k%nlabels],**args)
                if types[k%ntypes]=='m':
                    ax[l,c].plot(x,np.abs(y[k%nfigy,:]),label=labels[k%nlabels],**args)
                if types[k%ntypes]=='l10':
                    ax[l,c].plot(x,10*np.log10(np.abs(y[k%nfigy,:])),label=labels[k%nlabels],**args)
                if types[k%ntypes]=='l20':
                    ax[l,c].plot(x,20*np.log10(np.abs(y[k%nfigy,:])),label=labels[k%nlabels],**args)
                if types[k%ntypes]=='re':
                    ax[l,c].plot(x,np.real(y[k%nfigy,:]),label=labels[k%nlabels],**args)
                if types[k%ntypes]=='im':
                    ax[l,c].plot(x,np.imag(y[k%nfigy,:]),label=labels[k%nlabels],**args)
                if types[0]=='gdn':
                    df  = x[1]-x[0]
                    ax[l,c].plot(x[0:-1],-np.diff(np.unwrap(np.angle(y[k%nfigy,:])))/(2*np.pi*df),label=labels[k%nlabels],**args)
                if types[0]=='gdm':
                    df  = x[1]-x[0]
                    ax[l,c].plot(x[0:-1],-0.3*np.diff(np.unwrap(np.angle(y[k%nfigy,:])))/(2*np.pi*df),label=labels[k%nlabels],**args)

                ax[l,c].set_xlabel(xlabels[k%nxlabels])
                ax[l,c].set_ylabel(ylabels[k%nylabels])
                ax[l,c].set_title(titles[k%ntitles])
                ax[l,c].legend()
                #ax[l,c].get_xaxis().set_visible(False)
                #ax[l,c].get_yaxis().set_visible(False)

    #plt.tight_layout()

    return(fig,ax)                  

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
    #ax.plot(p[0,0,:], p[1,0,:], p[2,0,:], label='parametric curve',color='b')
    p = pb[:,:,np.newaxis] + \
            R*(np.cos(theta[np.newaxis,np.newaxis,:])*pn[:,:,np.newaxis] + \
               np.sin(theta[np.newaxis,np.newaxis,:])*qn[:,:,np.newaxis])
    #ax.plot(p[0,0,:], p[1,0,:], p[2,0,:], label='parametric curve',color='b')
    t = np.arange(0,1,0.05)
    p = pa[:,:,np.newaxis]+t[np.newaxis,np.newaxis,:]*(pb[:,:,np.newaxis]-pa[:,:,np.newaxis])
    ax.plot(p[0,0,:], p[1,0,:], p[2,0,:], label='parametric curve',color='b')
    

if (__name__ == "__main__"):
    plt.ion()
    doctest.testmod()
