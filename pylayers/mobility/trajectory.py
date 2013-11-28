import numpy as np 
import scipy as sp 
import pdb
import matplotlib.pyplot as plt 
import pylayers.util.pyutil as pyu
from matplotlib.path import Path 
import pandas as pd 


class Trajectory(pd.DataFrame):
    """  Define a trajectory

    This class derives from pandas.DataFrame

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
    rescale
    plot 


    """
    def __init__(self,t=[],pt=np.vstack((np.arange(0,10,0.01),np.zeros(1000))).T,unit='s'):
        """ initialization 
        """
        if t ==[]:
            t = np.arange(0,1,0.01)
        td = pd.to_datetime(t,unit=unit)
        #v = np.vstack((pt[1:,:]-pt[0:-1,:],np.array([np.nan,np.nan])))
        #a = np.vstack((v[1:,:]-v[0:-1,:],np.array([np.nan,np.nan])))
        v = pt[1:,:]-pt[0:-1,:]
        a = v[1:,:]-v[0:-1,:]
        d = np.sqrt(np.sum(v*v,axis=1))
        s = np.cumsum(d)
        #vy = self['y'][1:].values-self['y'][0:-1].values
        pd.DataFrame.__init__(self,{'t':td[:-2],
                                    'x':pt[:-2,0],
                                    'y':pt[:-2,1],
                                    'vx':v[:-1,0],
                                    'vy':v[:-1,1],
                                    'ax':a[:,0],
                                    'ay':a[:,1],
                                    's' :s[:-1]},columns=['t','x','y','vx','vy','ax','ay','s'])

        N = len(t) 
        self.tmin = t[0]
        self.tmax = t[N-2]
        #self.tmin = t[0].second + t[0].microsecond/1e6
        #self.tmax = t[N-2].second + t[N-2].microsecond/1e6
        self.ttime = self.tmax-self.tmin
        self.dtot = self['s'].values[-1]
        self.meansp = self.dtot/self.ttime

        if np.shape(pt)[1]>2:
            self['z'] = pt[:,2]

    def __repr__(self):

        dtot = self['s'].values[-1]
        T = self.tmax-self.tmin
        st = ''
        st = st+'t (s) : '+ str(self.tmin)+':'+str(self.tmax)+'\n'
        st = st+'d (m) : '+ str(dtot)+'\n'
        st = st+'Vmoy (m/s) : '+ str(dtot/T)+'\n'
        return(st)
    
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
        pt = self.space()
        t = Trajectory(t=newtime,pt=pt)
        return(t)

        

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
        lt = self['t']
        t  = np.array(map(lambda x : x.value,lt))
        conv = 10**(unit-9)
        t = t * conv
        return (t)

    def plot(self,fig=[],ax=[],Nlabels=5):
        """ plot trajectory

        Parameters
        ----------

        fig 
        ax 
        Nlabels : int 
        
        Examples
        --------
        
        .. plot::
            :include-source:

            >>> t = np.arange(0,10,0.01)
            >>> x = 2*t*np.cos(t)
            >>> y = 3*t*np.sin(t) 
            >>> pt =np.vstack((x,y)).T
            >>> traj = Trajectory(t,pt)
            >>> traj.plot()
            >>> plt.show()

        """

        if fig==[]:
            fig = plt.gcf()
        if ax == []:
            ax = plt.gca()

        ax.plot(self['x'],self['y'])
        for k in np.linspace(0,len(self),Nlabels,endpoint=False):
            k = int(k)
            ax.text(self['x'][k],self['y'][k],str(self['t'][k].strftime("%M:%S")))
            ax.plot(self['x'][k],self['y'][k],'*r')

        plt.xlabel('x (meters')
        plt.ylabel('y (meters')

        return fig,ax 

def importsn(_filename='pos.csv'):
    """ import simulnet csv file
    
    Parameters
    ----------

    filename : string 
        default 'pos.csv'

    Returns
    -------

    lt : list of trajectory

    """
    filename = pyu.getlong(_filename,'save_data')
    dt = pd.read_csv(filename)
    dtk = dt.keys()
    N = len(dtk)
    Ntraj = (N-1)/3
    lt = []
    for it in range(Ntraj):
        pt = np.vstack((dt[dtk[3*it+1]].values,dt[dtk[3*it+2]].values))
        lt.append(Trajectory(dt['time'].values,pt=pt.T,unit='s'))
    return(lt)    

if __name__ == '__main__':

    t = np.arange(0,10,0.01)
    x = 2*t*np.cos(t)
    y = 3*t*np.sin(t) 
    pt =np.vstack((x,y)).T

    traj = Trajectory(t,pt)
    traj.plot()
    plt.show()

   
