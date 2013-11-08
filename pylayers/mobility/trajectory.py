import numpy as np 
import scipy as sp 
import pdb
import matplotlib.pyplot as plt 
from matplotlib.path import Path 
import pandas as pd 


class Trajectory(pd.DataFrame):
    def __init__(self,t=np.arange(0,10,0.01),pt=np.vstack((np.arange(0,10,0.01),np.zeros(1000))).T,unit='s'):
        t = pd.to_datetime(t,unit=unit)
        #v = np.vstack((pt[1:,:]-pt[0:-1,:],np.array([np.nan,np.nan])))
        #a = np.vstack((v[1:,:]-v[0:-1,:],np.array([np.nan,np.nan])))
        v = pt[1:,:]-pt[0:-1,:]
        a = v[1:,:]-v[0:-1,:]
        d = np.sqrt(np.sum(v*v,axis=1))
        s = np.cumsum(d)
        #vy = self['y'][1:].values-self['y'][0:-1].values
        pd.DataFrame.__init__(self,{'t':t[:-2],
                                    'x':pt[:-2,0],
                                    'y':pt[:-2,1],
                                    'vx':v[:-1,0],
                                    'vy':v[:-1,1],
                                    'ax':a[:,0],
                                    'ay':a[:,1],
                                    's' :s[:-1]},columns=['t','x','y','vx','vy','ax','ay','s'])

        N = len(t) 
        self.tmin = t[0].second+t[0].microsecond/1e6
        self.tmax = t[N-2].second+t[N-2].microsecond/1e6
        if np.shape(pt)[1]>2:
            self['z'] = pt[:,2]
    
    
    def time(self,unit=0):
        """
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

if __name__ == '__main__':

    t = np.arange(0,10,0.01)
    x = 2*t*np.cos(t)
    y = 3*t*np.sin(t) 
    pt =np.vstack((x,y)).T

    traj = Trajectory(t,pt)
    traj.plot()
    plt.show()
