from pylayers.util.project import *
import pylayers.util.pyutil as pyu
from pylayers.util.utilnet import str2bool
from pylayers.gis.layout import Layout
from pylayers.antprop.multiwall import *
from pylayers.network.model import *


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as m
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
        self.rxopt = dict(self.config.items('rx'))
        self.showopt=dict(self.config.items('show'))

        self.L=Layout(self.layoutopt['filename'])
        self.model=Model(f=eval(self.plm['f']),
                         rssnp=eval(self.plm['rssnp']),
                         d0=eval(self.plm['d0']),
                         sigrss=eval(self.plm['sigrss']))
        self.xstep = eval(self.gridopt['xstep'])
        self.ystep = eval(self.gridopt['ystep'])
        # transitter section
        self.tx = np.array((eval(self.txopt['x']),eval(self.txopt['y'])))
        self.ptdbm = eval(self.txopt['ptdbm'])
        self.framelengthbytes = eval(self.txopt['framelengthbytes'])

        # receiver section
        self.rxsens = eval(self.rxopt['sensitivity'])
        kBoltzmann = 1.3806503e-23
        self.bandwidthmhz = eval(self.rxopt['bandwidthmhz'])
        self.temperaturek = eval(self.rxopt['temperaturek'])
        self.noisefactordb = eval(self.rxopt['noisefactordb'])

        Pn = (10**(self.noisefactordb/10.)+1)*kBoltzmann*self.temperaturek*self.bandwidthmhz*1e3
        self.pndbm = 10*np.log10(Pn)+60

        self.show = str2bool(self.showopt['show'])

        try:
            self.L.Gt.nodes()
        except:
            pass
        try:
            self.L.dumpr('t')
        except:
            self.L.buildGt()
            self.L.dumpw('t')
        self.creategrid()


    def creategrid(self):
        """create a grid
            create a grid for various evaluation

        """
        mi=np.min(self.L.Gs.pos.values(),axis=0)+0.01
        ma=np.max(self.L.Gs.pos.values(),axis=0)-0.01
        x=np.linspace(mi[0],ma[0],self.xstep)
        y=np.linspace(mi[1],ma[1],self.ystep)
        self.grid=np.array((list(np.broadcast(*np.ix_(x, y)))))




    def cover(self):
        """ start the coverage calculation

        Examples
        --------
        .. plot::
            :include-source:

            >>> from pylayers.antprop.coverage import *
            >>> C=Coverage()
            >>> C.cover()
            >>> C.showPr()

        """
        self.Lwo,self.Lwp,self.Edo,self.Edp=Loss0_v2(self.L,self.grid,self.model.f,self.tx)
        self.freespace = PL(self.grid,self.model.f,self.tx)
        self.prdbm = self.ptdbm - self.freespace - self.Lwo


    def showPr(self,rxsens=True,nfl=True):
        """ show the map of received power

        Parameters
        ----------

            rxsens : bool
                clip the map with rx sensitivity set in self.rxsens
            fnl : bool
                clip the map with noise floor set in self.pndbm

        """

        fig=plt.figure()
        fig,ax=self.L.showGs(fig=fig)
        l=self.grid[0,0]
        r=self.grid[-1,0]
        b=self.grid[0,1]
        t=self.grid[-1,-1]



#        tCM = plt.cm.get_cmap('jet')
#        tCM._init()
#        alphas = np.abs(np.linspace(.0,1.0, tCM.N))
#        tCM._lut[:-3,-1] = alphas
        title='Map of received power'

        cdict = {
        'red'  :  ((0., 0.5, 0.5), (1., 1., 1.)),
        'green':  ((0., 0.5, 0.5), (1., 1., 1.)),
        'blue' :  ((0., 0.5, 0.5), (1., 1., 1.))
        }
        #generate the colormap with 1024 interpolated values
        my_cmap = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)



        if rxsens :

            ### values between the rx sensitivity and noise floor
            mcPrf = np.ma.masked_where((self.prdbm > self.rxsens) 
                                     & (self.prdbm < self.pndbm),self.prdbm)
            cov1 = ax.imshow(mcPrf.reshape((self.xstep,self.ystep)).T,
                             extent=(l,r,b,t),cmap = my_cmap,
                             vmin=self.rxsens,origin='lower')

            ### values above the sensitivity
            mcPrs = np.ma.masked_where(self.prdbm < self.rxsens,self.prdbm)
            cov = ax.imshow(mcPrs.reshape((self.xstep,self.ystep)).T,
                            extent=(l,r,b,t),
                            cmap = 'jet',
                            vmin=self.rxsens,origin='lower')
            title=title + '\n black : PrdBm < rx sensitivity'

        else :
            cov=ax.imshow(self.prdbm.reshape((self.xstep,self.ystep)).T,
                          extent=(l,r,b,t),
                          cmap = 'jet',
                          vmin=self.PndBm,origin='lower')

        if nfl:
            ### values under the noise floor 
            ### we first clip the value below he noise fllor
            cl = np.nonzero(self.prdbm<=self.pndbm)
            cPr = self.prdbm
            cPr[cl] = self.pndbm
            mcPruf = np.ma.masked_where(cPr > self.pndbm,cPr)
            cov2 = ax.imshow(mcPruf.reshape((self.xstep,self.ystep)).T,
                             extent=(l,r,b,t),cmap = 'binary',
                             vmax=self.pndbm,origin='lower')
            title=title + '\n white : PrdBm < ' + str(self.pndbm) + 'dBm'


        ax.scatter(self.tx[0],self.tx[1],linewidth=0)

        ax.set_title(title)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(cov,cax)
        if self.show:
            plt.show()


    def showTransistionRegion(self):
        """
        Notes
        -----
        See  : Analyzing the Transitional Region in Low Power Wireless Links
                  Marco Zuniga and Bhaskar Krishnamachari

        """
        frameLength = self.framelengthbytes
        PndBm = self.pndbm
        gammaU = 10*np.log10(-1.28*np.log(2*(1-0.9**(1./(8*frameLength)))))
        gammaL = 10*np.log10(-1.28*np.log(2*(1-0.1**(1./(8*frameLength)))))
        PrU = PndBm + gammaU
        PrL = PndBm + gammaL

        fig=plt.figure()
        fig,ax = self.L.showGs(fig=fig)
        l = self.grid[0,0]
        r = self.grid[-1,0]
        b = self.grid[0,1]
        t = self.grid[-1,-1]

        zones = np.zeros(len(self.prdbm))
        uconnected  = np.nonzero(self.prdbm>PrU)[0]
        utransition = np.nonzero((self.prdbm < PrU)&(self.prdbm > PrL))[0]
        udisconnected = np.nonzero(self.prdbm < PrL)[0]

        zones[uconnected] = 1
        zones[utransition] = (self.prdbm[utransition]-PrL)/(PrU-PrL)
        cov = ax.imshow(zones.reshape((self.xstep,self.ystep)).T,
                             extent=(l,r,b,t),cmap = 'BuGn',origin='lower')

        title='PDR region'
        ax.scatter(self.tx[0],self.tx[1],linewidth=0)

        ax.set_title(title)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(cov,cax)
        if self.show:
            plt.show()

    def showEdo(self):
        """ show
        map of Loss from orthogonal field

        """

        fig=plt.figure()
        fig,ax=self.L.showGs(fig=fig)
        l=self.grid[0,0]
        r=self.grid[-1,0]
        b=self.grid[0,1]
        t=self.grid[-1,-1]
        cov=ax.imshow(self.Edo.reshape((self.xstep,self.ystep)).T,extent=(l,r,b,t),origin='lower')
        ax.scatter(self.tx[0],self.tx[1],linewidth=0)
        ax.set_title('Map of excess delay from orthogonal field')
        fig.colorbar(cov)
        if self.show:
            plt.show()

    def showLo(self):
        """ map Losses for orthogonal field
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
        """ map Losses for parallel field
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



if (__name__ == "__main__"):
    C=Coverage()
    C.cover()
    C.showPr()
    C.L.dumpr()
    sigar,sig=C.L.signature(C.grid[2],C.tx)

