from pylayers.util.project import *
import pylayers.util.pyutil as pyu
from pylayers.gis.layout import Layout
from pylayers.antprop.multiwall import *
from pylayers.network.model import *

import numpy as np
import matplotlib.pyplot as plt
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
        L
            a Layout

        method
            the chosen method for evaluating losses
            solely multiwall for now
        model
            a pylayers.network.model object
        xstep
            x step for grid
        ystep
            y step for grid
        tx
            transmitter position
    """

    def __init__(self,L=Layout('TA-Office.str'),model=Model(),method='multiwall',tx=[]):

        self.L=L
        self.method=method
        self.model=model
        self.xstep=10.
        self.ystep=10.
        try:
            self.L.Gt.nodes()
        except:
            self.L.buildGt()
        self.creategrid()
        if tx==[]:
            self.tx=self.grid[0]
        else:
            self.tx=tx




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
        >>> C=Coverage(tx=np.array((3,3)))
        >>> C.cover()
        >>> C.showPr()
        >>> plt.show()

        """
        self.Lwo,self.Lwp=Loss0_v2(self.L,self.grid,self.model.f,self.tx)
        self.Pr=-self.Lwo-self.model.PL0

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
        cov=ax.imshow(self.Pr.reshape((self.xstep,self.ystep)).T,extent=(l,r,b,t))
        ax.scatter(self.tx[0],self.tx[1],linewidth=0)
        ax.set_title('Map of received power')
        fig.colorbar(cov)

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
        cov=ax.imshow(self.Lwo.reshape((self.xstep,self.ystep)).T,extent=(l,r,b,t))
        ax.scatter(self.tx[0],self.tx[1],linewidth=0)
        ax.set_title('Map of Loss from orthogonal field')
        fig.colorbar(cov)

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
        cov=ax.imshow(self.Lwp.reshape((self.xstep,self.ystep)).T,extent=(l,r,b,t))
        ax.scatter(self.tx[0],self.tx[1],linewidth=0)
        ax.set_title('Map of Loss from parallel field')
        fig.colorbar(cov)

if (__name__ == "__main__"):

    C=Coverage(tx=np.array((12,12)))
    C.cover()
    C.showPr()
    plt.show()

