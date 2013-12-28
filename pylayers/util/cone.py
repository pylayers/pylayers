#-*- coding:Utf-8 -*-
"""
Module : Cone

"""
import numpy as np
import shapely as shp
import matplotlib.pyplot as plt
import pylayers.util.geomutil as geu
from matplotlib.path import Path
import matplotlib.patches as patches
import pdb


class Cone(object):

    def __init__(self, bv=np.array([[1,0],[0,1]]), apex=np.array([0, 0])):
        """
        bv : np.array(Ndim,Ndim)
                basis vector
        """
        self.apex = apex
        self.bv = bv / np.sqrt(np.sum(bv * bv, axis=0))
        self.gamma = np.dot(self.bv[:,0],self.bv[:,1])
        self.degenerated = False

    def belong(self, p):
        """ test if p belongs to Cone

        Parameters
        ----------

        p  : np.array (Ndim x Npoints)

        Returns
        -------

        b : np.array boolean (1xNpoints)

        """

        # Ndim x Npoints
        if not self.degenerated:
            pt  = p - self.apex[:,np.newaxis]
            puv = np.sum(self.bv[:,:,np.newaxis]*pt[:,np.newaxis,:],axis=0)
        
            alpha = puv[0,:]-self.gamma*puv[1,:] 
            beta  = puv[1,:]-self.gamma*puv[0,:]

            b = (beta>0)&(alpha>0)
        else:
            a0 = self.seg0[:,0]
            b0 = self.seg0[:,1]

            if self.bv[0,0]<>0:
                slope = self.bv[1,0]/self.bv[0,0]
                y0 = a0[1]-slope*a0[0]
                y1 = b0[1]-slope*b0[0]
                b = (p[1,:] > slope*p[0,:] + min(y0,y1) ) & (p[1,:]<slope*p[0,:]+max(y0,y1) ) 
            else:    
                b = (p[0,:] >  min(a0[0],b0[0]) ) & (p[0,:]< max(a0[0],b0[0]) )
        return(b)

    def above(self, p):
        bo1 = self.belong(p)
        pb = p[:,bo1]
        if self.bv[0,1]<>0:
            slope1 = self.bv[1,1]/self.bv[0,1]
            b1 = self.bv[1,1] - slope1*self.bv[0,1]
            bo2 = pb[1,:] > slope1*pb[0,:]+b  
        else:
            bo2 = pb[1,:] > self.seg1[1,0]   
           
        return(bo1,bo2)    
        
    def from2segs(self,seg0,seg1):
        """ creates a Cone from 2 segments

        Parameters
        ----------

        seg0 : 2 x 2  (Ndim x Npoints)
        seg1 : 2 x 2 

        """ 
        # bv : (4,1)
        self.seg0 = seg0
        self.seg1 = seg1

        a0 = seg0[:,0]
        b0 = seg0[:,1]
        a1 = seg1[:,0]
        b1 = seg1[:,1]

        # check segment orientation (crossing)
        if not (geu.ccw(a0,b0,b1) ^
            geu.ccw(b0,b1,a1) ):
            v0 = (b1 - a0)[:,np.newaxis]
            v1 = (a1 - b0)[:,np.newaxis]
        else:    
            v0 = (a1 - a0)[:,np.newaxis]
            v1 = (b1 - b0)[:,np.newaxis]
        
        bv = np.hstack((v0,v1))
        #pdb.set_trace()
        self.bv = bv / np.sqrt(np.sum(bv * bv, axis=0))
        self.gamma = np.dot(self.bv[:,0],self.bv[:,1])

        if abs(self.gamma)>1-1e-15:
            self.degenerated=True
        else:    
            a0u = np.dot(a0,self.bv[:,0])
            a0v = np.dot(a0,self.bv[:,1])
            b0u = np.dot(b0,self.bv[:,0])
            b0v = np.dot(b0,self.bv[:,1])

            kb  = ((b0v-a0v)-self.gamma*(b0u-a0u))/(self.gamma*self.gamma-1)
            self.apex = b0 + kb*self.bv[:,1]

    def show(self, **kwargs):

        defaults = {'length': 25}
        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if 'seg1' not in self.__dict__:    
            verts = [tuple(self.apex),
                     tuple(self.apex + kwargs['length'] * self.bv[:, 0]),
                     tuple(self.apex + kwargs['length'] * self.bv[:, 1]),
                     tuple(self.apex)
                    ]
            codes = [Path.MOVETO,
                Path.LINETO,
                Path.LINETO,
                Path.CLOSEPOLY,
            ]

        else:
            a0 = self.seg0[:,0]
            b0 = self.seg0[:,1]
            a1 = self.seg1[:,0]
            b1 = self.seg1[:,1]

            if not(self.degenerated):
                #verts = [tuple(self.apex),
                #         tuple(a1),
                #         tuple(b1),
                #         tuple(self.apex)
                #        ]
                verts = [tuple(self.apex),
                     tuple(self.apex + kwargs['length'] * self.bv[:, 0]),
                     tuple(self.apex + kwargs['length'] * self.bv[:, 1]),
                     tuple(self.apex)
                    ]
                codes = [Path.MOVETO,
                Path.LINETO,
                Path.LINETO,
                Path.CLOSEPOLY,
                ]
            else:
                if (geu.ccw(a0,b0,b1) ^
                    geu.ccw(b0,b1,a1) ):
                    verts = [tuple(b0),
                             tuple(a1),
                             tuple(b1),
                             tuple(a0),
                             tuple(b0)
                        ]

                else:
                    verts = [tuple(b0),
                             tuple(b1),
                             tuple(a1),
                             tuple(a0),
                             tuple(b0)
                        ]

                codes = [Path.MOVETO,
                Path.LINETO,
                Path.LINETO,
                Path.LINETO,
                Path.CLOSEPOLY,
                ]


        

        path = Path(verts, codes)

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        if 'seg1' in self.__dict__:
            ax.plot([a0[0],b0[0]],[a0[1],b0[1]],lw=2,color='b')
            ax.plot([a1[0],b1[0]],[a1[1],b1[1]],lw=2,color='r')
        patch = patches.PathPatch(path, facecolor='orange', lw=2, alpha=0.3)
        ax.add_patch(patch)
        ax.axis('equal')
        # ax.set_xlim(-2,2)
        # ax.set_ylim(-2,2)

        return(fig, ax)
