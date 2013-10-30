import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
from matplotlib.path import Path 
import pandas as pd 


class Trajectory(pd.DataFrame):
    def __init__(self,t,pt,unit='s'):
        index = pd.to_datetime(t,unit=unit)
        pd.DataFrame.__init__(self,{'x':pt[:,0],'y':pt[:,1]},index=index)
        if np.shape(pt)[1]>2:
            self['z'] = pt[:,2]

    def plot(self,fig=[],ax=[],Nlabels=5):

        if fig==[]:
            fig = plt.gcf()
        if ax == []:
            ax = plt.gca()

        ax.plot(self['x'],self['y'])
        for k in np.linspace(0,len(self),Nlabels,endpoint=False):
            k = int(k)
            ax.text(self['x'][k],self['y'][k],str(self.index[k].strftime("%H:%M:%S")))

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
