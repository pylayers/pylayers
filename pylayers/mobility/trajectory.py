import numpy as np
import scipy as sp
import pdb
import matplotlib.pyplot as plt
import pylayers.util.pyutil as pyu
from matplotlib.path import Path
from pylayers.util.project import *
from pylayers.gis.layout import Layout
import pandas as pd
import copy
import time
import doctest
from matplotlib.widgets import Slider,CheckButtons



class Trajectory(pd.DataFrame):
    """  Define a trajectory

    This class derives from pandas.DataFrame. It handles a full 3D trajectory
    description.

    A trajectory is time-stamped and contains 3D coordinates of position, velocity and acceleration.

    Attributes
    ----------

    tmin : float
    tmax : float
    tttimr :
    dtot :
    meansp :

    Methods
    -------

    time
    space
    update
    rescale
    generate
    distance
    plot
    replay


    """
    def __init__(self, df = {}):
        """ initialization
        """
        super(Trajectory,self).__init__(df)
        self.has_values=self.update()



    def __repr__(self):
        try:
            dtot = self['s'].values[-1]
            T = self.tmax-self.tmin
            st = ''
            st = st+'t (s) : '+ str(self.tmin)+':'+str(self.tmax)+'\n'
            st = st+'d (m) : '+ str(dtot)+'\n'
            st = st+'Vmoy (m/s) : '+ str(dtot/T)+'\n'
        except:
            st = 'void Trajectory'
        return(st)

    def update(self):
        """ update class member data

        This method updates the following data members

        + tmin (s)
        + tmax (s)
        + ts   time step in second
        + ttime (s) trajectory duration
        + measnsp

        Returns
        -------

        bool :
            True if Trajectory has values, False otherwise

        """

        if len(self.values) != 0:
            N = len(self.index)
            self.tmin = self.index.min().value*1e-9
            self.tmax = self.index.max().value*1e-9
            self.ts = (self.index[1].value*1e-9)-(self.index[0].value*1e-9)

            self.ttime = self.tmax-self.tmin
            self.dtot = self['s'].values[-1]
            self.meansp = self.dtot/self.ttime
            return True
        else :
            return False


    def generate(self,t=np.linspace(0,10,50),pt=np.vstack((np.sin(np.linspace(0,3,50)),np.linspace(0,10,50),np.random.randn(50),)).T,unit='s'):
        """
        Generate a trajectroy from a numpy array

        Parameters
        ----------

        pt : np.ndarray:
            (npt x x x y)
            with
                npt : number of samples
                x : x values
                y : y values
        t = np.ndarray
            (1 x npt)

        unit : str
            time unity ('s'|'ns',...)

        Examples
        --------

        >>> from pylayers.mobility.trajectory import *
        >>> traj = Trajectory()
        >>> traj.generate()
        >>> traj.plot()


        """


        npt = len(t)
        td = pd.to_datetime(t,unit=unit)
        # velocity vector
        v = pt[1:,:]-pt[0:-1,:]
        # acceleration vector
        a = v[1:,:]-v[0:-1,:]
        #
        d = np.sqrt(np.sum(v*v,axis=1))
        s = np.cumsum(d)
        s[-1] = 0
        s = np.roll(s,1)

        df = {'x':pt[:-2,0],
            'y':pt[:-2,1],
            'z':pt[:-2,2],
            'vx':v[:-1,0],
            'vy':v[:-1,1],
            'vz':v[:-1,2],
            'ax':a[:,0],
            'ay':a[:,1],
            'az':a[:,2],
            's':s[:-1]}
        super(Trajectory,self).__init__(df,columns=['x','y','z','vx','vy','vz','ax','ay','az','s'],index=td[:-2])
        self.update()
        return self


    def rescale(self,speedkmph=3):
        """ same length but specified speed

        Parameters
        ----------

        speedkmph : float
            targeted mean speed in km/h

        Returns
        -------

        t : rescaled trajectory

        """
        speedms = speedkmph/3.6
        factor  = speedms/self.meansp
        newtime = self.time()/factor
        pt = self.space(ndim=3)
        t = copy.copy(self)
        t.generate(t=newtime,pt=pt)
        return(t)


    def distance(self,tk):
        """ recover distance at time tk

        Parameters
        ----------

        tk :

        Example
        -------

        >>> from pylayers.mobility.trajectory import *
        >>> bc = body.Body()

        """
        t = self.time()
        u = np.where((t>=tk-self.ts/2.)&(t<=tk+self.ts/2.))[0][0]

        return(self['s'][u])

    def space(self,ndim=2):
        """ extract space information

        Parameters
        ----------

        ndim : int
            number of dimensions (default 2)

        Returns
        -------

        pt : nd.array()

        """
        if ndim==2:
            pt = np.vstack((self['x'].values,self['y'].values)).T
        if ndim==3:
            pt = np.vstack((self['x'].values,self['y'].values,self['z'].values)).T
        return(pt)

    def time(self,unit=0):
        """ extract time

        Parameters
        ----------

        unit : integer
            default 0 (s) - 3 (ms) 6 (mus) 9 (ns)

        Returns
        -------

        t : nd.array
           time in 10**-unit  s

        """
        lt = self.index
        t  = np.array(map(lambda x : x.value,lt))
        conv = 10**(unit-9)
        t = t * conv
        return (t)

    def plot(self,fig=[],ax=[],Nlabels=5,typ='plot',L=[]):
        """ plot trajectory

        Parameters
        ----------

        fig
        ax
        Nlabels : int
        typ : 'plot'|'scatter'
        L : pylayers.gis.layout.Layout object to be displayed

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.mobility.trajectory import *
            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> t = np.arange(0,10,0.01)
            >>> x = 2*t*np.cos(t)
            >>> y = 3*t*np.sin(t)
            >>> z = 0*t
            >>> pt =np.vstack((x,y,z)).T
            >>> traj = Trajectory()
            >>> traj.generate(t,pt)
            >>> f,a = traj.plot()
            >>> plt.show()

        """

        if fig==[]:
            fig = plt.gcf()
        if ax == []:
            ax = plt.gca()

        if L!=[]:
            if isinstance(L,Layout):
                fig,ax=L.showGs(fig=fig,ax=ax)

        if typ == 'plot':
            ax.plot(self['x'],self['y'])
        elif typ == 'scatter':
            ax.scatter(self['x'],self['y'])

        for k in np.linspace(0,len(self),Nlabels,endpoint=False):
            k = int(k)

            ax.text(self['x'][k],self['y'][k],str(self.index[k].strftime("%M:%S")))
            ax.plot(self['x'][k],self['y'][k],'*r')

        plt.xlabel('x (meters)')
        plt.ylabel('y (meters)')

        return fig,ax

    def replay(self,fig=[],ax=[],Nlabels=5,typ='plot',L=[],speed=1,**kwargs):
        """
            replay a trajectory

        Parameters
        ----------

        fig
        ax
        Nlabels : int
            default 5
        typ : string
            'plot'|'scatter'
        L : pylayers.gis.layout.Layout
            Layout for body to be displayed in
        speed : float
            speed ratio

        """

        plt.ion()
        if fig==[]:
            fig = plt.gcf()
        if ax == []:
            ax = plt.gca()

        limkwargs = copy.copy(kwargs)
        if 'c' in kwargs :
            limkwargs.pop('c')
        if 'color'in kwargs:
            limkwargs.pop('c')
        limkwargs['marker']='*'
        limkwargs['s']=20

        if ('m' or 'marker') not in kwargs:
            kwargs['marker']='o'
        if ('c' or 'color') not in kwargs:
            kwargs['color']='b'


        if L!=[]:
            if isinstance(L,Layout):
                fig,ax=L.showGs(fig=fig,ax=ax,**kwargs)

        labels = np.linspace(0,len(self),Nlabels,endpoint=True).tolist()

        for ik,k in enumerate(self.index):
            time.sleep(1/(1.*speed))
            ax.scatter(self['x'][ik],self['y'][ik],**kwargs)
            plt.title(str(self.index[ik].time())[:11].ljust(12),loc='left')

            if ik > labels[0]:
                ax.text(self['x'][ik],self['y'][ik],str(self.index[ik].strftime("%M:%S")))
                ax.scatter(self['x'][ik],self['y'][ik],**limkwargs)
                labels.pop(0)
            plt.draw()


        plt.ioff()
        # for k in :
        #     k = int(k)
        #     ax.text(self['x'][k],self['y'][k],str(self.index[k].strftime("%M:%S")))
        #     ax.plot(self['x'][k],self['y'][k],'*r')
        #     plt.draw()

def importsn(_filename='pos.csv'):
    """
    ****DEPRECATED
    import simulnet csv file

    ****DEPRECATED

    Parameters
    ----------

    filename : string
        default 'pos.csv'

    Returns
    -------

    lt : list of trajectory

    """
    filename = pyu.getlong(_filename,pstruc['DIRNETSAVE'])
    dt = pd.read_csv(filename)
    dtk = dt.keys()
    N = len(dtk)
    Ntraj = (N-1)/3
    lt = []
    for it in range(Ntraj):
        x = dt[dtk[3*it+1]].values
        y = dt[dtk[3*it+2]].values
        z = np.zeros(len(x))
        pt = np.vstack((x,y,z))
        T=Trajectory()
        lt.append(T.generate(t=dt['time'].values,pt=pt.T,unit='s'))

    return(lt)

def importh5(self,_filename='simulnet_TA-Office.h5'):
        """ import simulnet h5 file

        Parameters
        ----------

        filename : string
            default simulnet + Layout_filename . h5

        Returns
        -------

        lt : list of trajectory

        """

        filename = pyu.getlong(_filename,pstruc['DIRNETSAVE'])
        fil = pd.HDFStore(filename)

        lt=[]

        for k in fil.keys():
            df = fil[k]
            df = df.set_index('t')
            v=np.array((df.vx.values,df.vy.values))
            d = np.sqrt(np.sum(v*v,axis=0))
            s = np.cumsum(d)
            df['s'] = s
            lt.append(Trajectory(df))
        fil.close()
        return lt

if __name__ == '__main__':
    plt.ion()
    doctest.testmod()




