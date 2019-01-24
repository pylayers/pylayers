import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio
from mayavi import mlab
import pdb
deg_to_rad = np.pi/180.
rad_to_deg = 180./np.pi

class OTA(object):
    """ Over The Air Simulator

    config = 0 : spherical distribution of probes

    """
    def __init__(self,**kwargs):
        typ = kwargs.pop('config','spherical')
        theta=kwargs.pop('theta',np.array([90,60,30]))
        phi = kwargs.pop('phi',np.array([45,135,225,315]))
        R = kwargs.pop('R',0.66)
        H = kwargs.pop('H',1.5)
        if typ=='spherical':
            theta = theta*deg_to_rad
            phi = phi*deg_to_rad
            x = (R*np.sin(theta[:, None])*np.cos(phi[None, :])).ravel()
            y = (R*np.sin(theta[:, None])*np.sin(phi[None, :])).ravel()
            z = (H + R*np.cos(theta[:, None])*np.ones(len(phi))[None, :] ).ravel()

        self.p = np.vstack((x, y, z))
        self.th = (theta[:, None]*np.ones(len(phi))[None, :]).ravel()
        self.ph = (np.ones(len(theta))[:,None]*phi[None,: ]).ravel()

    def __repr__(self):
        st = ''
        for k in range(self.th.shape[0]):
            st = st + str(k+1)+' '+str(self.th[k]*rad_to_deg)+' '+str(self.ph[k]*rad_to_deg)+'\n'
        return(st)

    def load(self, filename):
        """ load an OTA file


        self.M: np.array
            (f x Nx x Ny x Ant) 
        self.vmax
        self.vmaxdB
        self.vmindB

        """
        U = spio.loadmat(filename)
        #print U.keys()
        self.fGHz = U['freq'].ravel()/1e9
        for k in range(12):
            key = 'E'+str(k+1)
            X = U[key][:,:,:,None]
            try:
                self.M = np.concatenate((self.M,X),axis=3)
            except:
                self.M = X
        self.vmax = np.max(np.abs(self.M))
        self.vmaxdB = 20*np.log10(np.max(np.abs(self.M)))
        self.vmindB = self.vmaxdB-30


    def set_grid(self,rx,ry,rz):
        """ set grid

        rx : range along x
        ry : range along y
        rz : range along z

        """

        try:
            del self.pg
        except:
            pass
        for x in rx:
            for y in ry:
                for z in rz:
                    pt = np.array([x,y,z])[None,:]
                    try:
                        self.pg = np.vstack((self.pg,pt))
                    except:
                        self.pg = pt


    def show(self, **kwargs):
        ax = kwargs.pop('ax', [])
        fGHz = kwargs.pop('fGHz', 2)
        ka = kwargs.pop('ka', 1)
        kind = kwargs.pop('kind', 'ad')
        config = kwargs.pop('config', True)
        grid = kwargs.pop('grid', True)
        label = kwargs.pop('label', False)
        colorbar = kwargs.pop('colorbar', False)
        alpha = kwargs.pop('alpha', 1)
        s = kwargs.pop('s', 10)


        if ax == []:
            ax = plt.gca()
        else:
            pass
        # determine frequency index
        if hasattr(self, 'fGHz'):
            abdf = np.abs(self.fGHz-fGHz)
            kf = np.where(abdf==np.min(abdf))
            fGHz = self.fGHz[kf]
            #print('Freq (GHz) : ',fGHz)

        # display config
        if config:
            if ka> 0:
                ax.plot(self.p[0,ka-1],self.p[1,ka-1],'or')
                if label:
                    ax.annotate(str(ka),xy=(self.p[0,ka-1],self.p[1,ka-1]),fontsize=18)
            else:
                ax.plot(self.p[0,:],self.p[1,:],'or')
                if label:
                    for k in range(self.p.shape[1]):
                        ax.annotate(str(k+1),xy=(self.p[0,k],self.p[1,k]),fontsize=18)

            for k in range(self.p.shape[1]):
                r = np.sqrt(self.p[0,k]**2+self.p[1,k]**2)
                t = np.linspace(0,2*np.pi,100)
                u = r*np.cos(t)
                v = r*np.sin(t)
                ax.plot(u,v,linewidth=0.5,color='blue')
        if hasattr(self,'fGHz'):
            if kind=='m':
                val = np.abs(self.M[kf,:,:,ka-1])
                vmax = self.vmax
                vmin = 0 
            if kind=='l20':
                val = 20*np.log10(np.abs(self.M[kf,:,:,ka-1]))
                vmax = self.vmaxdB
                vmin = self.vmindB
            if kind=='ar':
                val = np.angle(self.M[kf,:,:,ka-1])
                vmax = np.pi
                vmin = -np.pi
            if kind=='ad':
                val = np.angle(self.M[kf,:,:,ka-1])*rad_to_deg
                vmax = 180
                vmin = -180
            #sca=ax.scatter(self.pg[:,0],self.pg[:,1],c=val.T[::-1,::-1],s=s,alpha=alpha,linewidth=0,vmin=vmin,vmax=vmax)
            sca=ax.scatter(self.pg[:,0],self.pg[:,1],c=val.ravel(),s=s,alpha=alpha,linewidth=0,vmin=vmin,vmax=vmax)
        elif hasattr(self,'pg'):
            sca=ax.scatter(self.pg[:,0],self.pg[:,1],c='k',s=3,alpha=0.5)
        else:
            pass
        if grid:
            ax.grid()
        ax.axis('equal')
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.axis('off')
        if colorbar:
            plt.colorbar(sca,ax=ax)
        fig = plt.gcf()
        return(fig,ax)



if __name__=="__main__":
    # specify the config
    theta = np.array([30,60,90])
    phi = np.array([45,135,225,315])
    ota = OTA(theta,phi)
    # specify the grid
    rx = np.arange(-0.07,0.07,0.01)
    ry = np.arange(-0.07,0.07,0.01)
    rz = np.array([1.5])
    ota.set_grid(rx,ry,rz)
