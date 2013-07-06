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
    """ Handle Layout Coverage

        Methods
        -------
        create grid()
            create a uniform grid for evaluating losses
        cover()
            run the coverage calculation
        showPower()
            display the map of received power
        showLoss()
            display the map of losses


        Attributes
        ----------
        All attributes are read from fileini ino the ini directory of the
        current project

        _fileini
            default coverage.ini

        L :  a Layout
        model : a pylayers.network.model object.
        nx    : number of point on x
        ny    : number of point on y
        tx    : transmitter position
        txpe  : transmitter power emmission level
        show  : boolean for automatic display power map

    """


    def __init__(self,_fileini='coverage.ini'):


        self.config = ConfigParser.ConfigParser()
        self.config.read(pyu.getlong(_fileini,pstruc['DIRSIMUL']))
        self.plm = dict(self.config.items('pl_model'))
        self.layoutopt = dict(self.config.items('layout'))
        self.gridopt = dict(self.config.items('grid'))
        self.txopt = dict(self.config.items('tx'))
        self.rxopt = dict(self.config.items('rx'))
        self.showopt = dict(self.config.items('show'))

        self.L = Layout(self.layoutopt['filename'])
        self.model = PLSmodel(f=eval(self.plm['fghz']),
                         rssnp=eval(self.plm['rssnp']),
                         d0=eval(self.plm['d0']),
                         sigrss=eval(self.plm['sigrss']))

        self.nx = eval(self.gridopt['nx'])
        self.ny = eval(self.gridopt['ny'])
        self.mode = eval(self.gridopt['full'])
        self.boundary = eval(self.gridopt['boundary'])
        # transitter section
        self.fGHz = eval(self.txopt['fghz'])
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

        self.creategrid(full=self.mode,boundary=self.boundary)

    def __repr__(self):
        st=''
        st= st+ 'tx :'+str(self.txopt) + '\n'
        st= st+ 'rx :'+str(self.rxopt) + '\n'

    def creategrid(self,full=True,boundary=[]):
        """ create a grid

        Parameters
        ----------
        full : boolean
            default (True) use all the layout area
        boundary : (xmin,ymin,xmax,ymax)
            if full is False the boundary is used

        """
        if full:
            mi=np.min(self.L.Gs.pos.values(),axis=0)+0.01
            ma=np.max(self.L.Gs.pos.values(),axis=0)-0.01
        else:
            mi = np.array([boundary[0],boundary[1]])
            ma = np.array([boundary[2],boundary[3]])

        x=np.linspace(mi[0],ma[0],self.nx)
        y=np.linspace(mi[1],ma[1],self.ny)
        self.grid=np.array((list(np.broadcast(*np.ix_(x, y)))))




    def cover(self):
        """ start the coverage calculation

        Examples
        --------
        .. plot::
            :include-source:

            >>> from pylayers.antprop.coverage import *
            >>> C = Coverage()
            >>> C.cover()
            >>> C.showPr()

        """
        self.Lwo,self.Lwp,self.Edo,self.Edp = Loss0_v2(self.L,self.grid,self.model.f,self.tx)
        self.freespace = PL(self.grid,self.model.f,self.tx)
        self.prdbmo = self.ptdbm - self.freespace - self.Lwo
        self.prdbmp = self.ptdbm - self.freespace - self.Lwp
        self.snro = self.prdbmo - self.pndbm
        self.snrp = self.prdbmp - self.pndbm


    def showEd(self,polarization='o'):
        """ show direct path excess of delay map

        Examples
        --------
        .. plot::
            :include-source:

            >>> from pylayers.antprop.coverage import *
            >>> C = Coverage()
            >>> C.cover()
            >>> C.showEdo()
        """

        fig=plt.figure()
        fig,ax=self.L.showGs(fig=fig)
        l=self.grid[0,0]
        r=self.grid[-1,0]
        b=self.grid[0,1]
        t=self.grid[-1,-1]
        if polarization=='o':
            cov=ax.imshow(self.Edo.reshape((self.nx,self.ny)).T,
                      extent=(l,r,b,t),
                      origin='lower')
            titre = "Map of LOS excess delay, polar orthogonal"
        if polarization=='p':
            cov=ax.imshow(self.Edp.reshape((self.nx,self.ny)).T,
                      extent=(l,r,b,t),
                      origin='lower')
            titre = "Map of LOS excess delay, polar parallel"

        ax.scatter(self.tx[0],self.tx[1],linewidth=0)
        ax.set_title(titre)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        clb = fig.colorbar(cov,cax)
        clb.set_label('excess delay (ns)')
        if self.show:
            plt.show()


    def showPower(self,rxsens=True,nfl=True,polarization='o'):
        """ show the map of received power

        Parameters
        ----------

        rxsens : bool
              clip the map with rx sensitivity set in self.rxsens
        nfl : bool
              clip the map with noise floor set in self.pndbm
        polarization : string
            'o'|'p'

        Examples
        --------
        .. plot::
            :include-source:

            >>> from pylayers.antprop.coverage import *
            >>> C = Coverage()
            >>> C.cover()
            >>> C.showPower()

        """

        fig=plt.figure()
        fig,ax=self.L.showGs(fig=fig)
        l=self.grid[0,0]
        r=self.grid[-1,0]
        b=self.grid[0,1]
        t=self.grid[-1,-1]

        if polarization=='o':
            prdbm=self.prdbmo
        if polarization=='p':
            prdbm=self.prdbmp

#        tCM = plt.cm.get_cmap('jet')
#        tCM._init()
#        alphas = np.abs(np.linspace(.0,1.0, tCM.N))
#        tCM._lut[:-3,-1] = alphas
        title='Map of received power - Pt = '+str(self.ptdbm)+' dBm'

        cdict = {
        'red'  :  ((0., 0.5, 0.5), (1., 1., 1.)),
        'green':  ((0., 0.5, 0.5), (1., 1., 1.)),
        'blue' :  ((0., 0.5, 0.5), (1., 1., 1.))
        }
        #generate the colormap with 1024 interpolated values
        my_cmap = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)



        if rxsens :

            ### values between the rx sensitivity and noise floor
            mcPrf = np.ma.masked_where((prdbm > self.rxsens) 
                                     & (prdbm < self.pndbm),prdbm)
            cov1 = ax.imshow(mcPrf.reshape((self.nx,self.ny)).T,
                             extent=(l,r,b,t),cmap = my_cmap,
                             vmin=self.rxsens,origin='lower')

            ### values above the sensitivity
            mcPrs = np.ma.masked_where(prdbm < self.rxsens,prdbm)
            cov = ax.imshow(mcPrs.reshape((self.nx,self.ny)).T,
                            extent=(l,r,b,t),
                            cmap = 'jet',
                            vmin=self.rxsens,origin='lower')
            title=title + '\n gray : Pr (dBm) < %.2f' % self.rxsens + ' dBm'

        else :
            cov=ax.imshow(prdbm.reshape((self.nx,self.ny)).T,
                          extent=(l,r,b,t),
                          cmap = 'jet',
                          vmin=self.pndbm,origin='lower')

        if nfl:
            ### values under the noise floor 
            ### we first clip the value below he noise floor
            cl = np.nonzero(prdbm<=self.pndbm)
            cPr = prdbm
            cPr[cl] = self.pndbm
            mcPruf = np.ma.masked_where(cPr > self.pndbm,cPr)
            cov2 = ax.imshow(mcPruf.reshape((self.nx,self.ny)).T,
                             extent=(l,r,b,t),cmap = 'binary',
                             vmax=self.pndbm,origin='lower')
            title=title + '\n white : Pr (dBm) < %.2f' % self.pndbm + ' dBm'


        ax.scatter(self.tx[0],self.tx[1],s=10,linewidth=0)

        ax.set_title(title)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        clb = fig.colorbar(cov,cax)
        clb.set_label('Power (dBm)')
        if self.show:
            plt.show()


    def showTransistionRegion(self,polarization='o'):
        """
        Notes
        -----
        See  : Analyzing the Transitional Region in Low Power Wireless Links
                  Marco Zuniga and Bhaskar Krishnamachari

        Examples
        --------
        .. plot::
            :include-source:

            >>> from pylayers.antprop.coverage import *
            >>> C = Coverage()
            >>> C.cover()
            >>> C.showTransitionRegion()

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

        if polarization=='o':
            prdbm=self.prdbmo
        if polarization=='p':
            prdbm=self.prdbmp

        zones = np.zeros(len(prdbm))
        uconnected  = np.nonzero(prdbm>PrU)[0]
        utransition = np.nonzero((prdbm < PrU)&(prdbm > PrL))[0]
        udisconnected = np.nonzero(prdbm < PrL)[0]

        zones[uconnected] = 1
        zones[utransition] = (prdbm[utransition]-PrL)/(PrU-PrL)
        cov = ax.imshow(zones.reshape((self.nx,self.ny)).T,
                             extent=(l,r,b,t),cmap = 'BuGn',origin='lower')

        title='PDR region'
        ax.scatter(self.tx[0],self.tx[1],linewidth=0)

        ax.set_title(title)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(cov,cax)
        if self.show:
            plt.show()

    def showLoss(self,polarization='o'):
        """ show losses map

        Parameters
        ----------
        polarization : string 
            'o'|'p'|'both'

        Examples
        --------
        .. plot::
            :include-source:

            >>> from pylayers.antprop.coverage import *
            >>> C = Coverage()
            >>> C.cover()
            >>> C.showLoss(polarization='o')
            >>> C.showLoss(polarization='p')
        """
        fig = plt.figure()
        fig,ax=self.L.showGs(fig=fig)
        l=self.grid[0,0]
        r=self.grid[-1,0]
        b=self.grid[0,1]
        t=self.grid[-1,-1]

        if polarization=='o':
            cov = ax.imshow(self.Lwo.reshape((self.nx,self.ny)).T,
                            extent=(l,r,b,t),
                            origin='lower')
            title = ('Map of losses, orthogonal (V) polarization') 
        if polarization=='p':
            cov = ax.imshow(self.Lwp.reshape((self.nx,self.ny)).T,
                            extent=(l,r,b,t),
                            origin='lower')
            title = ('Map of losses, parallel (H) polarization') 

        ax.scatter(self.tx[0],self.tx[1],linewidth=0)
        ax.set_title(title)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(cov,cax)

        if self.show:
            plt.show()



if (__name__ == "__main__"):
    plt.ion()
    C=Coverage()
    C.cover()
    C.showPower()
    C.showLoss(polarization='o')
    C.showLoss(polarization='p')
    C.showTransistionRegion(polarization='o')
    C.showEd(polarization='o')
#    C.L.dumpr()
#    sigar,sig=C.L.signature(C.grid[2],C.tx)

