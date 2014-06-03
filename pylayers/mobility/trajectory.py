"""

Trajectories Class
==================
.. autosummary::
    :toctree: generated/

    Trajectories.__init__
    Trajectories.__repr__
    Trajectories.loadh5
    Trajectories.ishow

Trajectory Class
================

.. autosummary::
    :toctree: generated/

    Trajectory.__init__
    Trajectory.__repr__
    Trajectory.update
    Trajectory.generate
    Trajectory.resample
    Trajectory.rescale
    Trajectory.replay
    Trajectory.distance
    Trajectory.space
    Trajectory.time
    Trajectory.plot

Utility Functions
=================

.. autosummary::
    :toctree: generated/

    importsn
    importh5

"""
import numpy as np
import scipy as sp
import pdb
import matplotlib.pyplot as plt
import pylayers.util.pyutil as pyu
from pylayers.util.project import *
from pylayers.gis.layout import Layout
import pandas as pd
import copy
import time
import doctest
from matplotlib.widgets import Slider, CheckButtons
import matplotlib.animation as animation



class Trajectories(list):
    """  Define a list of trajectory


    """
    def __init__(self):
        """ initialization
        """
        super(list, self).__init__()
        self.name = []
        self.typ = []
        self.ID = []
        self.t = []

    def __repr__(self):
        """
        """

        try:
            s = 'Trajectories performed in Layout : ' + self.Lfilename + '\n\n'
            for a in self:
                s = s + a.__repr__()
                s = s + '\n'
        except:
            s = 'Issue in Trajectories. Are you sure any Trajectory is loaded ?'
        return s


    def append(self,obj):
        """ overload list.append
        """

        super(Trajectories,self).append(obj)
        self.name.append(obj.name)
        self.typ.append(obj.typ)
        self.ID.append(obj.ID)


    def pop(self,idx=-1):
        """ overloaded list.pop
        """

        super(Trajectories,self).pop(idx)
        self.name.pop(idx)
        self.typ.pop(idx)
        self.ID.pop(idx)
        
        
    def saveh5(self, L=Layout(), filename = ''):
        
        layfile = L.filename.split('.')[0]
        if filename == '':
            store = pd.HDFStore('/private/staff/i/bi/mmhedhbi/Bureau/P1/netsave/traj_'+layfile+'.h5','w') 
        else :
            store = pd.HDFStore('/private/staff/i/bi/mmhedhbi/Bureau/P1/netsave/traj_'+filename+'.h5','w')   
        
        for traj in self:
            
            t = traj.time()*1e9
            td = pd.to_datetime(t)
            
            if traj.typ != 'ap':
                data  = pd.DataFrame(
                    {'t': td, 'x': traj.x.values, 'y': traj.y.values,
                                           'vx': traj.vx.values, 'vy': traj.vy.values,
                                           'ax': traj.ax.values, 'ay': traj.ay.values,
                     }, columns=['t', 'x', 'y',  'vx', 'vy', 'ax', 'ay'], index=td )
            else:
                data  = pd.DataFrame(
                    {'t': td, 'x': traj.x.values, 'y': traj.y.values,'z': traj.z.values, 
                                'vx': traj.vx.values, 'vy': traj.vy.values,
                                 'ax': traj.ax.values, 'ay': traj.ay.values},
                        columns=['t', 'x', 'y',  'z', 'vx', 'vy', 'ax', 'ay' ], index=td )
            
                
            store.put(str(traj.ID),data) 

           

            store.get_storer(str(traj.ID)).attrs.typ = traj.typ
            store.get_storer(str(traj.ID)).attrs.name = traj.name
            store.get_storer(str(traj.ID)).attrs.ID = str(traj.ID)
            store.get_storer(str(traj.ID)).attrs.layout = L.filename
            #saving metadata

        store.close()
                
            
        
           




    def loadh5(self, _filename='simulnet_TA-Office.h5',append =False):
        """ import simulnet h5 file

        Parameters
        ----------

        filename : string
            default simulnet + Layout_filename . h5
        append : boolean
            if True : append new trajectories to preexisting ones

        Returns
        -------

        lt : list of trajectory

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.mobility.trajectory import *
            >>> T = Trajectories()
            >>> T.loadh5()

        """

        filename = pyu.getlong(_filename, pstruc['DIRNETSAVE'])
        if os.path.exists(filename):
            fil = pd.HDFStore(filename)
        else:
            raise NameError(filename + ' not founded')
        if not append:
            [self.pop(0) for i in range(len(self))]
        for k in fil.keys():
           
            df = fil[k]
            
            df = df.set_index('t')
            ID = fil.get_storer(k).attrs.ID
            name = fil.get_storer(k).attrs.name
            typ = fil.get_storer(k).attrs.typ
            layout = fil.get_storer(k).attrs.layout
            v=np.array((df.vx.values,df.vy.values))
            d = np.sqrt(np.sum(v*v,axis=0))
            s = np.cumsum(d)
            df['s'] = s
            self.append(Trajectory(df=df,ID=ID,name=name,typ=typ))
        fil.close()
        self.Lfilename = layout
        self.time()


    def resample(self, sf=2, tstart = -1):
        """ resample trajectories

        Parameters
        ----------

        sf : int
            sampling factor
        tstart : float
            new start time (must be > original start time).
            if tstart = -1 : original start conserved

        Returns
        -------

        T : Trajectories
            new trajectories object updated

        """
        T=Trajectories()
        for t in self:
            if t.typ != 'ap':
                T.append(t.resample(sf=sf, tstart=tstart))
            else:
                T.append(t)
        T.Lfilename = self.Lfilename
        T.time()
        return T

    def time(self,unit='s'):
        """ extract time from a trajectory

        Parameters
        ----------

        unit : integer
            default 0 (s) - 3 (ms) 6 (mus) 9 (ns)

        Returns
        -------

        update self.t


        """

        ut = np.where(np.array(self.typ) == 'ag')[0][0]
        self.t = self[ut].time()

    def replay(self, fig=[], ax=[], **kwargs):
        """ replay a trajectory

        Parameters
        ----------

        fig
        ax

        speed : float
            speed ratio

        """

        # plt.ion()
        if fig==[]:
            fig = plt.gcf()
        if ax == []:
            ax = plt.gca()

        limkwargs = copy.copy(kwargs)
        if 'c' in kwargs:
            limkwargs.pop('c')
        if 'color'in kwargs:
            limkwargs.pop('c')
        limkwargs['marker'] = '*'
        limkwargs['s'] = 20

        if ('m' or 'marker') not in kwargs:
            kwargs['marker'] = 'o'
        if ('c' or 'color') not in kwargs:
            kwargs['color'] = 'b'

        L=Layout(self.Lfilename)
        fig, ax = L.showG('s',fig=fig, ax=ax, **kwargs)

        time=self[0].time()


        line, = ax.plot([], [], 'ob', lw=2)
        time_template = 'time = %.1fs'
        time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

        def init():
            line.set_data([],[])
            time_text.set_text('')
            return line, time_text

        def animate(it):
            X=[]
            Y=[]
            for t in self:
                if t.typ == 'ag':
                    X.append(t['x'].values[it])
                    Y.append(t['y'].values[it])
            line.set_data(X,Y)
            time_text.set_text(time_template%(time[it]))
            return line, time_text

        ani = animation.FuncAnimation(fig, animate, np.arange(1, len(time)),
            interval=25, blit=True, init_func=init)
        plt.show()
       


    def ishow(self):
        """
            interactive show of trajectories

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.mobility.trajectory import *
            >>> T=Trajectories()
            >>> T.loadh5()
            >>> T.ishow()

        """

        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.2, left=0.3)

        t = np.arange(0, len(self[0].index), self[0].ts)
        L = Layout(self.Lfilename)
        fig, ax = L.showG('s', fig=fig, ax=ax)

        valinit = 0
        lines = []
        labels = []
        colors = "bgrcmykw"

        for iT, T in enumerate(self):
            if T.typ == 'ag':
                lines.extend(ax.plot(T['x'][0:valinit],T['y'][0:valinit], 'o',
                             color=colors[iT], visible=True))
                labels.append(T.name + ':' + T.ID)
            else:
                lines.extend(ax.plot(T['x'][0], T['y'][0], '^', ms=12,
                             color=colors[iT], visible=True))
                labels.append(T.name + ':' + T.ID)

        t = self[0].time()

        # init boolean value for visible in checkbutton
        blabels = [True]*len(labels)


        ########
        # slider
        ########
        slider_ax = plt.axes([0.1, 0.1, 0.8, 0.02])
        slider = Slider(slider_ax, "time", self[0].tmin, self[0].tmax,
                        valinit=valinit, color='#AAAAAA')

        def update(val):
            if val >= 1:
                pval=np.where(val>t)[0]
                ax.set_title(str(self[0].index[pval[-1]].time())[:11].ljust(12),
                             loc='left')
                for iT, T in enumerate(self):
                    if T.typ == 'ag':
                        lines[iT].set_xdata(T['x'][pval])
                        lines[iT].set_ydata(T['y'][pval])
                fig.canvas.draw()
        slider.on_changed(update)


        ########
        # choose
        ########
        rax = plt.axes([0.02, 0.5, 0.3, 0.2], aspect='equal')
        # check (ax.object, name of the object , bool value for the obsject)
        check = CheckButtons(rax, labels, tuple(blabels))

        def func(label):
            i = labels.index(label)
            lines[i].set_visible(not lines[i].get_visible())
            fig.canvas.draw()

        check.on_clicked(func)
        fig.canvas.draw()
        plt.show(fig)


class Trajectory(pd.DataFrame):
    """  Define a trajectory

    This class derives from pandas.DataFrame. It handles a full 3D trajectory.

    A trajectory is time-stamped and contains 3D coordinates of p
    position, velocity and acceleration.

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
    def __init__(self, df={}, ID=0, name='', typ=''):
        """ initialization
        """
        super(Trajectory, self).__init__(df)
        self.ID = ID
        self.name = name
        self.typ = typ
        self.has_values = self.update()

    def __repr__(self):
        try:
            # total distance
            dtot = self['s'].values[-1]
            # total time
            T = self.tmax-self.tmin
            st = ''
            typ = self.typ
            if typ == 'ag':
                string ='Trajectory of agent ' + self.name + ' with ID ' + self.ID
            else :
                string ='Access point ' + self.name + ' with ID ' + self.ID
            st = st + string + '\n'
            st = st + '-'*len(string) + '\n'

            if self.typ == 'ag':
                st = st+'t (s) : '+ str("%3.2f" %self.tmin)+" : "+ str("%3.2f" % self.ts) +" : " +str("%3.2f" % self.tmax)+'\n'
                st = st+'dtot (m) : '+ str("%3.2f" %dtot)+'\n'
                st = st+'Vmoy (m/s) : '+ str("%3.2f" % (dtot/T))+'\n'
            else :
                st = st+'t (s) : '+ str("%3.2f" %self.tmin) + '\n'
                st = st+'Vmoy (m/s) : '+ str(self['vx'].values[0]) +'\n'
            st = st + str(self.head(2)) + '\n'
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
            
            self.tmin = self.index.min().value*1e-9
            self.tmax = self.index.max().value*1e-9
            try:
                self.ts = (self.index[-1].value*1e-9)-(self.index[-2].value*1e-9)
            except:
                self.ts = np.nan

            self.ttime = self.tmax-self.tmin
            self.dtot = self['s'].values[-1]
            self.meansp = self.dtot/self.ttime
            return True
        else :
            return False



    #def generate(self,ID = 1, name = '',t=np.linspace(0,10,50),pt=np.vstack((np.sin(np.linspace(0,3,50)),np.linspace(0,10,50),np.random.randn(50),)).T,unit='s', sf = 1):
    def generate(self,**kwargs):
        """
        Generate a trajectory from a numpy array

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

        Id : Agent Id

        name : Agent Name

        unit : str
            time unity ('s'|'ns',...)

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.mobility.trajectory import *
            >>> traj = Trajectory()
            >>> traj.generate()
            >>> traj.plot()



        """
        defaults  = { 'ID': 1,
                     'name': 'MyNameIsNoBody',
                     'typ':'ag',
                     't': np.linspace(0,10,50),
                     'pt': np.vstack((np.sin(np.linspace(0,3,50)),np.linspace(0,10,50),np.random.randn(50),)).T,
                     'unit': 's',
                     'sf': 1
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value
        
        t = kwargs['t']
        if len(t) < 2:
            raise AttributeError('Trajectory.generate requieres at least 3 time stamps')
        pt = kwargs['pt']
        npt = len(t)
        td = pd.to_datetime(t,unit=kwargs['unit'])
        # velocity vector
        v = (pt[1:, :]-pt[0:-1, :])/(t[1]-t[0])
        # acceleration vector
        a = (v[1:, :]-v[0:-1, :])/(t[1]-t[0])
        #
        d = np.sqrt(np.sum(v[:,0:2]*v[:, 0:2], axis=1))
        s = np.cumsum(d)*(t[1]-t[0])
        s[-1] = 0
        s = np.roll(s,1)

        df = {
            'x': pt[:-2, 0],
            'y': pt[:-2, 1],
            'z': pt[:-2, 2],
            'vx': v[:-1, 0],
            'vy': v[:-1, 1],
            'vz': v[:-1, 2],
            'ax': a[:, 0],
            'ay': a[:, 1],
            'az': a[:, 2],
            's': s[:-1]}
        super(Trajectory, self).__init__(df, columns=['x', 'y', 'z', 'vx', 'vy',
                                         'vz', 'ax', 'ay', 'az', 's'],
                                        index=td[:-2])
        self.ID = kwargs['ID']
        self.name = kwargs['name']
        self.typ = kwargs['typ']
        self.update()
        return self

    def resample(self, sf=2, tstart=-1):
        """ resample trajectory

        Parameters
        ----------

        sf : float
            sampling factor
        tstart : float
            new start time (must be > original start time).
            if tstart = -1 : original start conserved

        Returns
        -------

        T : Trajectory
            resampled trajectory

        """
       
        t = self.time()
        x = self.space()[:, 0]
        y = self.space()[:, 1]
        fx = sp.interpolate.interp1d(t, x)
        fy = sp.interpolate.interp1d(t, y)
        if tstart == -1:
            tstart = t[0]
        else:
            if t[0] <= tstart:
                tstart = tstart
            else :
                raise AttributeError('tstart < tmin')

        tstep = (t[1]-t[0])/sf
        # need to add at least 3 values gor generate to estomate acceleration
        tnew = np.arange(tstart, t[-1], tstep)

        # generate requieres 3 measures at least 
        xnew = fx(tnew)
        ynew = fy(tnew)
        T = Trajectory()

        T.generate(ID=self.ID,
                   name=self.name,
                   typ=self.typ,
                   t=tnew,
                   pt=np.vstack((xnew,ynew,np.random.randn(len(tnew)),)).T,
                   unit='s',
                   sf=sf)
        
        T.update()
       
        return T



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
        factor = speedms/self.meansp
        newtime = self.time()/factor
        pt = self.space(ndim=3)
        t = copy.copy(self)
        t.generate(ID=self.ID, name=self.name, t=newtime, pt=pt)
        return(t)

    def distance(self,tk):
        """ recover distance at time tk

        Parameters
        ----------

        tk : float
            time value in seconds

        Example
        -------

        >>> from pylayers.mobility.trajectory import *
        >>> T = Trajectory()
        >>> T.generate()
        >>> T.distance(2)

        """
        t = self.time()
        #pdb.set_trace()
        u = np.where((t >= tk-self.ts/2.) & (t <= tk+self.ts/2.))[0][0]
        

        return(self['s'][u])

    def space(self, ndim=2):
        """ extract space information

        Parameters
        ----------

        ndim : int
            number of dimensions (default 2)

        Returns
        -------

        pt : nd.array()

        """
        if ndim == 2:
            pt = np.vstack((self['x'].values, self['y'].values)).T
        if ndim == 3:
            pt = np.vstack((self['x'].values, self['y'].values,self['z'].values)).T
        return(pt)

    def time(self, unit=0):
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
        t  = np.array(map(lambda x: x.value, lt))
        conv = 10**(unit-9)
        t = t * conv
        return (t)

    def plot(self, fig=[], ax=[], Nlabels=5, typ='plot', L=[]):
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

        if L != []:
            if isinstance(L,Layout):
                fig, ax = L.showGs(fig=fig, ax=ax)

        if typ == 'plot':
            ax.plot(self['x'], self['y'])
        elif typ == 'scatter':
            ax.scatter(self['x'], self['y'])

        for k in np.linspace(0, len(self), Nlabels, endpoint=False):
            k = int(k)

            ax.text(self['x'][k], self['y'][k], str(self.index[k].strftime("%M:%S")))
            ax.plot(self['x'][k], self['y'][k], '*r')

        plt.xlabel('x (meters)')
        plt.ylabel('y (meters)')

        return fig, ax

    def replay(self, fig=[], ax=[], Nlabels=5, typ='plot', L=[], speed=1, **kwargs):
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

        # plt.ion()
        if fig==[]:
            fig = plt.gcf()
        if ax == []:
            ax = plt.gca()

        limkwargs = copy.copy(kwargs)
        if 'c' in kwargs:
            limkwargs.pop('c')
        if 'color'in kwargs:
            limkwargs.pop('c')
        limkwargs['marker'] = '*'
        limkwargs['s'] = 20

        if ('m' or 'marker') not in kwargs:
            kwargs['marker'] = 'o'
        if ('c' or 'color') not in kwargs:
            kwargs['color'] = 'b'

        if L != []:
            if isinstance(L,Layout):
                fig, ax = L.showG('s',fig=fig, ax=ax, **kwargs)

        labels = np.linspace(0, len(self), Nlabels, endpoint=True).tolist()

        line, = ax.plot([], [], 'ob', lw=2)
        time_template = 'time = %.1fs'
        time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

        def init():
            line.set_data([],[])
            time_text.set_text('')
            return line, time_text

        def animate(it):
            thisx = [-100,self['x'].values[it]]
            thisy = [-100,self['y'].values[it]]
            line.set_data(thisx, thisy)
            time_text.set_text(time_template%(self.time()[it]))
            return line, time_text

        ani = animation.FuncAnimation(fig, animate, np.arange(1, len(self.time())),
            interval=25, blit=True, init_func=init)
        plt.show()
        # for ik, k in enumerate(self.index):
        #     time.sleep(1/(1.*speed))
        #     ax.scatter(,**kwargs)
        #     plt.title(str(self.index[ik].time())[:11].ljust(12), loc='left')

        #     if ik > labels[0]:
        #         ax.text(self['x'].values[ik], self['y'].values[ik], str(self.index[ik].strftime("%M:%S")))
        #         ax.scatter(self['x'].values[ik], self['y'].values[ik], **limkwargs)
        #         labels.pop(0)
        #     plt.draw()

        # plt.ioff()
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
    filename = pyu.getlong(_filename, pstruc['DIRNETSAVE'])
    dt = pd.read_csv(filename)
    dtk = dt.keys()
    N = len(dtk)
    Ntraj = (N-1)/3
    lt = []
    for it in range(Ntraj):
        x = dt[dtk[3*it+1]].values
        y = dt[dtk[3*it+2]].values
        z = np.zeros(len(x))
        pt = np.vstack((x, y, z))
        T=Trajectory()
        lt.append(T.generate(t=dt['time'].values, pt=pt.T, unit='s'))

    return(lt)

# def importh5(_filename='simulnet_TA-Office.h5'):
#         """ import simulnet h5 file

#         Parameters
#         ----------

#         filename : string
#             default simulnet + Layout_filename . h5

#         Returns
#         -------

#         lt : list of trajectory

#         """

#         filename = pyu.getlong(_filename,pstruc['DIRNETSAVE'])
#         fil = pd.HDFStore(filename)

#         lt=[]

#         for k in fil.keys():
#             df = fil[k]
#             df = df.set_index('t')
#             v=np.array((df.vx.values,df.vy.values))
#             d = np.sqrt(np.sum(v*v,axis=0))
#             s = np.cumsum(d)
#             df['s'] = s
#             lt.append(Trajectory(df))
#         fil.close()
#         return lt

if __name__ == '__main__':
    plt.ion()
    doctest.testmod()


