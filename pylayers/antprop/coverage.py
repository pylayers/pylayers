from pylayers.util.project import *
import pylayers.util.pyutil as pyu
from pylayers.util.utilnet import str2bool
from pylayers.gis.layout import Layout
from pylayers.antprop.multiwall import *
from pylayers.network.model import *


import numpy as np
import matplotlib.pyplot as plt
import ConfigParser
import pdb

class Coverage(object):
    """ Simulation Class

        Methods
        -------
  
        create grid()
            create a uniform grid for evaluating losses
        cover()
            run the coverage 
        showPr()
            display the map of received power
        showLo()
            display the map of loss from orthogonal field
        showLp()
            display the map of loss from parallel field


        Attributes
        ----------
    
        All attributes are read from fileini ino the ini director of projection

        _fileini
            default coverage.ini

        L
            a Layout
        model
            a pylayers.network.model object. 
        xstep
            x step for grid
        ystep
            y step for grid
        tx
            transmitter position
        txpe
            transmitter power emmission level
        show 
            Boolean to automatic display power map

    """


    def __init__(self,_fileini='coverage.ini'):


        self.config = ConfigParser.ConfigParser()
        self.config.read(pyu.getlong(_fileini,pstruc['DIRSIMUL']))
        self.plm = dict(self.config.items('pl_model'))
        self.layoutopt = dict(self.config.items('layout'))
        self.gridopt = dict(self.config.items('grid'))
        self.txopt = dict(self.config.items('tx'))
        self.showopt=dict(self.config.items('show'))

        self.L=Layout(self.layoutopt['filename'])
        self.model=Model(f=eval(self.plm['f']),rssnp=eval(self.plm['rssnp']),d0=eval(self.plm['d0']),sigrss=eval(self.plm['sigrss']))
        self.xstep = eval(self.gridopt['xstep'])
        self.ystep = eval(self.gridopt['ystep'])

        self.tx = np.array((eval(self.txopt['x']),eval(self.txopt['y'])))
        self.txpe = eval(self.txopt['pe'])

        self.show = str2bool(self.showopt['show'])

        try:
            self.L.Gt.nodes()
        except:
            self.L.buildGt()
        self.creategrid()


    def creategrid(self):
        """create a grid
            create a grid for evaluating losses

        """
        mi=np.min(self.L.Gs.pos.values(),axis=0)+0.01
        ma=np.max(self.L.Gs.pos.values(),axis=0)-0.01
        x=np.linspace(mi[0],ma[0],self.xstep)
        y=np.linspace(mi[1],ma[1],self.ystep)
        self.grid=np.array((list(np.broadcast(*np.ix_(x, y)))))




    def cover(self):
        """ cover
        start the cover estimation

        Examples
        --------
    
        .. plot::
        :include-source:

        >>> import matplotlib.pyplot as plt 
        >>> from pylayers.antprop.coverage import * 
        >>> C=Coverage()
        >>> C.cover()
        >>> C.showPr()

        """
        self.Lwo,self.Lwp=Loss0_v2(self.L,self.grid,self.model.f,self.tx)
        self.Pr=self.txpe-self.model.PL0-self.Lwo

    def showPr(self):
        """ show
        show the map of received power
        """

        fig=plt.figure()
        fig,ax=self.L.showGs(fig=fig)
        l=self.grid[0,0]
        r=self.grid[-1,0]
        b=self.grid[0,1]
        t=self.grid[-1,-1]
        cov=ax.imshow(self.Pr.reshape((self.xstep,self.ystep)).T,extent=(l,r,b,t),vmin=-120,origin='lower')
        ax.scatter(self.tx[0],self.tx[1],linewidth=0)
        ax.set_title('Map of received power')
        fig.colorbar(cov)
        if self.show:
            plt.show()

    def showLo(self):
        """ show
        map of Loss from orthogonal field
        """

        fig=plt.figure()
        fig,ax=self.L.showGs(fig=fig)
        l=self.grid[0,0]
        r=self.grid[-1,0]
        b=self.grid[0,1]
        t=self.grid[-1,-1]
        cov=ax.imshow(self.Lwo.reshape((self.xstep,self.ystep)).T,extent=(l,r,b,t),origin='lower')
        ax.scatter(self.tx[0],self.tx[1],linewidth=0)
        ax.set_title('Map of Loss from orthogonal field')
        fig.colorbar(cov)
        if self.show:
            plt.show()


    def showLp(self):
        """ show
        map of Loss from parallel field
        """

        fig=plt.figure()
        fig,ax=self.L.showGs(fig=fig)
        l=self.grid[0,0]
        r=self.grid[-1,0]
        b=self.grid[0,1]
        t=self.grid[-1,-1]
        cov=ax.imshow(self.Lwp.reshape((self.xstep,self.ystep)).T,extent=(l,r,b,t),origin='lower')
        ax.scatter(self.tx[0],self.tx[1],linewidth=0)
        ax.set_title('Map of Loss from parallel field')
        fig.colorbar(cov)
        if self.show:
            plt.show()

#if (__name__ == "__main__"):

#    C=Coverage(tx=np.array((0.3,0.3)))
#    C.cover()
#    fig=plt.figure()
#    fig,ax=C.L.showGs(fig=fig)
##    fig=plt.figure()
##    fig2,ax2=C.L.showGs(fig=fig)

#    l=C.grid[0,0]
#    r=C.grid[-1,0]
#    b=C.grid[0,1]
#    t=C.grid[-1,-1]
#    cc=ax.imshow(C.Pr.reshape((C.xstep,C.ystep)).T,extent=(l,r,b,t),vmin=-120,origin='lower')
#    ax.scatter(C.tx[0],C.tx[1],linewidth=0)
#    fig.colorbar(cc)
##    ax2.imshow(C.Pr.reshape((C.xstep,C.ystep)).T,extent=(l,r,b,t),vmin=-120)
##    ax2.scatter(C.tx[0],C.tx[1],linewidth=0)

#    plt.show()
