#-*- coding:Utf-8 -*-
r"""

Class Cone
==========

The following conventions are adopted

+ A Cone has an **apex** which is a point in the plane.
+ A cone has two vectors which define the cone aperture. Those two vectors can
always been distinguished as a starting vector (u) and a ending vector (v).

The cone region is defined by the convex angular sector going from starting
vector  :math:`\mathbf{u}` to ending vector :math:`\mathbf{v}`
rotating in the plane in folllowing the trigonometric convention.
The modulus of the cross product between :math:`\mathbf{u}` and :math:`\mathbf{v}` is positive.

:math:`\mathbf{u} \times \mathbf{v} = \alpha \mathbf{z} \;\; \textrm{with} \;\;\alpha > 0`

.. autosummary::
    :toctree:

    Cone.__init__
    Cone.upd_angle
    Cone.belong_seg
    Cone.aboveseg
    Cone.outside_point
    Cone.belong_point2
    Cone.belong_point
    Cone.above
    Cone.fromptseg
    Cone.from2segs
    Cone.from2csegs
    Cone.show


"""
import numpy as np
import shapely as shp
import matplotlib.pyplot as plt
import pylayers.util.geomutil as geu
from pylayers.util.project import *
from matplotlib.path import Path
import matplotlib.patches as patches
import pdb
import logging


class Cone(PyLayers):

    def __init__(self, a=np.array([1,0]), b = np.array([0,1]), apex=np.array([0, 0])):
        """
        a : np.array (,2)
                basis vector
        b : np.array (,2)        
        apex : np.array (,2)
        """

        self.apex = apex
        # normalizing cone vectors 
        an = a/np.sqrt(np.dot(a,a))
        bn = b/np.sqrt(np.dot(b,b))

        if np.cross(an,bn) > 0:
            self.u = an
            self.v = bn
        else:  
            self.u = bn
            self.v = an
        
        # -1 < gamma < 1
        self.dot = np.dot(self.u,self.v)
        self.cross = np.cross(self.u,self.v)

        if self.cross<>0:
            self.degenerated = False
        else:
            self.degenerated = True    

        # update cone angle and probability
        self.upd_angle()
        

    def upd_angle(self):
        """update cone angle attribute 
           and associated probability of the Cone object
        """
        
        self.angle = np.arccos(self.dot)
        self.pcone = self.angle/(1.0*np.pi)


    def belong_seg(self,pta,phe,prob=True):
        """ test if segment belong to cone

        Parameters
        ----------

        pta : np.array (2xNseg)
        phe : np.array (2xNseg)

        Returns
        -------

        typ   : int 
            0 : no visibility 
            1 : full visibility 
            2 : he.v 
            3 : ta.v
            4 : ta.u  
            5 : he.u  
            6 : inside 
        proba : float
            geometric probability 

            
        Notes
        -----

        A segment belongs to the cone if not all termination points 
        lie in the same side outside the cone.

        """
        vc  = (self.u+self.v)/2
        #vcn = vc/np.sqrt(np.dot(vc,vc))
        w = vc/np.sqrt(np.dot(vc,vc))
        w = w.reshape(2,1)
        #w  = np.array([vcn[1],-vcn[0]])

        ptama = pta - self.apex.reshape(2,1)
        phema = phe - self.apex.reshape(2,1)

        dtaw = np.sum(ptama*w,axis=0)
        dhew = np.sum(phema*w,axis=0)

        blta = (dtaw>0)
        blhe = (dhew>0)
        #if 'seg1' in self.__dict__:
        #    pa =  self.seg1[:,0].reshape(2,1)
        #    pb = (self.seg1[:,0]+w).reshape(2,1)
        #else:
        #    pa = self.apex.reshape(2,1)
        #    pb = pa+w.reshape(2,1)
        #blta = geu.isleft(pa,pb,pta)
        #blhe = geu.isleft(pa,pb,phe)
        # segment candidate for being above segment 1 (,Nseg)
        boup = blta & blhe
        # type of segment
        if prob:
            proba = np.zeros(np.shape(pta)[1])
        else :
            proba =[]
        typ   = np.zeros(np.shape(pta)[1])
        # is tail out ? bo1 | bo2  
        # btaol : boolean tail out left 
        # btaor : boolean tail out right 
        # bheol : boolean head out left 
        # bheor : boolean head out right #
        # among upper segment check position wrt cone
        #btaol,btaor = self.outside_point(pta)
        #bheol,bheor = self.outside_point(phe)
        btaor,btaol = self.outside_point(pta)
        bheor,bheol = self.outside_point(phe)
        # tail and head are they out cone on the same side ? 
        # if the two termination points are not on the same side of the cone
        # --> segment is in.
        # boin = (~((btaol&bheol)|(btaor&bheor)))&boup
        # full interception (proba to reach = 1) 
        bfull = ((btaol&bheor)|(btaor&bheol))&boup
        if prob :
            proba[bfull] = 1
        typ[bfull] = 1

        #(he-apex).v
        btalhein  = (btaol & ~bheol & ~bheor)&boup
        if (prob and not (btalhein==False).all()):
            v2  = phe[:,btalhein]-self.apex.reshape(2,1)
            vn2 = v2/np.sqrt(np.sum(v2*v2,axis=0))
            vvn2 = np.dot(self.v,vn2)
            vvn2 = np.minimum(vvn2,np.ones(len(vvn2)))
            vvn2 = np.maximum(vvn2,-np.ones(len(vvn2)))
            pr2 = np.arccos(vvn2)/np.arccos(self.dot)
            proba[btalhein] = pr2
        typ[btalhein] = 2

        #(tai-apex).v
        bheltain  = (bheol & ~btaol & ~btaor)&boup
        if (prob and not (bheltain==False).all()):
            v3  = pta[:,bheltain]-self.apex.reshape(2,1)
            vn3 = v3/np.sqrt(np.sum(v3*v3,axis=0))
            vvn3 = np.dot(self.v,vn3)
            vvn3 = np.minimum(vvn3,np.ones(len(vvn3)))
            vvn3 = np.maximum(vvn3,-np.ones(len(vvn3)))
            pr3 = np.arccos(vvn3)/np.arccos(self.dot)
            proba[bheltain] = pr3
        typ[bheltain] = 3

        #ta.u
        bhertain  = (bheor & ~btaol & ~btaor)&boup
        if (prob and not(bhertain==False).all()):
            v4  = pta[:,bhertain]-self.apex.reshape(2,1)
            vn4 = v4/np.sqrt(np.sum(v4*v4,axis=0))
            vvn4 = np.dot(self.v,vn4)
            vvn4 = np.minimum(vvn4,np.ones(len(vvn4)))
            vvn4 = np.maximum(vvn4,-np.ones(len(vvn4)))
            pr4 = np.arccos(vvn4)/np.arccos(self.dot)
            proba[bhertain] = pr4
        typ[bhertain] = 4

        #he.u
        btarhein  = (btaor & ~bheol & ~bheor)&boup
        if (prob and not(btarhein==False).all()):
            v5  = phe[:,btarhein]-self.apex.reshape(2,1)
            vn5 = v5/np.sqrt(np.sum(v5*v5,axis=0))
            vvn5 = np.dot(self.v,vn5)
            vvn5 = np.minimum(vvn5,np.ones(len(vvn5)))
            vvn5 = np.maximum(vvn5,-np.ones(len(vvn5)))
            pr5 = np.arccos(vvn5)/np.arccos(self.dot)
            proba[btarhein] = pr5
        typ[btarhein] = 5

        #ta.he
        btainhein  = (~btaol & ~btaor & ~bheol & ~bheor)&boup
        if (prob and not (btainhein==0).all()):
            va  = pta[:,btainhein]-self.apex.reshape(2,1)
            vb  = phe[:,btainhein]-self.apex.reshape(2,1)
            vna = va/np.sqrt(np.sum(va*va,axis=0))
            vnb = vb/np.sqrt(np.sum(vb*vb,axis=0))
            vnab = np.sum(vna*vnb,axis=0)
            vnab = np.minimum(vnab,np.ones(len(vnab)))
            vnab = np.maximum(vnab,-np.ones(len(vnab)))
            pr6 = np.arccos(vnab)/np.arccos(self.dot)
            proba[btainhein] = pr6
        typ[btainhein] = 6

        return(typ,proba)

    def aboveseg(self):
        """
        """
        vc  = (self.u+self.v)/2
        vcn = vc/np.sqrt(dot(vc,vc))
        w  = np.array([vcn[1],-vcn[0]])
        self.pa =  self.seg1[:,0].reshape(2,1)  
        self.pb = (self.seg1[:,0]+w).reshape(2,1)

    def outside_point(self,p):
        """ check if p is outside cone

        Parameters
        ----------

        p : np.array  (2xNp)

        Returns
        -------

        ~b1 & ~b2 : boolean (outside on the left)  (,Np)
        b1 & b2 : boolean (outside on the right)  (,Np) 

        Examples
        --------

        Notes
        -----

        If one of the two output booleans is True the point is outside 
        There are 2 output bits but only 3 states due to (uv) convention.

            v    u 
        p    \  /       lv & lu
              \/

             \p /
              \/        ~lv & lu

             \  /  p
              \/        ~lu & ~lv


        """

        a = self.apex[:,np.newaxis]
        b = a + self.u.reshape(2,1)
        c = a + self.v.reshape(2,1)

        p0a0 = p[0,:]-a[0,:]
        p1a1 = p[1,:]-a[1,:]
        lu = ((b[0,:]-a[0,:])* p1a1 - ((b[1,:]-a[1,:])* p0a0 ))>0
        lv = ((c[0,:]-a[0,:])* p1a1 - ((c[1,:]-a[1,:])* p0a0 ))>0

        return(~lu & ~lv , lu & lv)

    def belong_point2(self,p):
        """
        """

        a = self.apex[:,np.newaxis]
        b = a + self.u.reshape(2,1)
        c = a + self.v.reshape(2,1)

        p1a1 = p[1,:]-a[1,:]
        p0a0 = p[0,:]-a[0,:]
        b1 = ((b[0,:]-a[0,:])* p1a1 - ((b[1,:]-a[1,:])* p0a0 ))>0
        b2 = ((c[0,:]-a[0,:])* p1a1 - ((c[1,:]-a[1,:])* p0a0 ))>0

        return(b1^b2)

    def belong_point(self, p):
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
            #puv = np.sum(self.bv[:,:,np.newaxis]*pt[:,np.newaxis,:],axis=0)

            #alpha = puv[0,:]-self.gamma*puv[1,:]
            #beta  = puv[1,:]-self.gamma*puv[0,:]
            pu = np.sum(self.u[:,np.newaxis]*pt,axis=0)
            pv = np.sum(self.v[:,np.newaxis]*pt,axis=0)

            alpha = pu-self.dot*pv
            beta  = pv-self.dot*pu

            b = (beta>0)&(alpha>0)
        else:
            a0 = self.seg0[:,0]
            b0 = self.seg0[:,1]

            if self.u[0]<>0:
                slope = self.u[1]/self.u[0]
                y0 = a0[1]-slope*a0[0]
                y1 = b0[1]-slope*b0[0]
                b = (p[1,:] > slope*p[0,:] + min(y0,y1) ) & (p[1,:]<slope*p[0,:]+max(y0,y1) )
            else:
                b = (p[0,:] >  min(a0[0],b0[0]) ) & (p[0,:]< max(a0[0],b0[0]) )
        return(b)

    def above(self, p):
        """ check if above
        Parameters
        ----------

        p :
        """
        bo1 = self.belong(p)
        pb = p[:,bo1]
        if self.v[1]<>0:
            slope1 = self.v[1]/self.v[0]
            b1 = self.v[1] - slope1*self.v[0]
            bo2 = pb[1,:] > slope1*pb[0,:]+b
        else:
            bo2 = pb[1,:] > self.seg1[1,0]

        return(bo1,bo2)

    def fromptseg(self,pt,seg):
        """ creates a Cone from one point and one segment

        Parameters
        ----------

        pt : nd.array (,2)

        seg : nd.array (2,2)

        """

        self.apex = pt 
        a = seg[:,0]
        b = seg[:,1]
        v0 = b - pt
        v1 = a - pt
        v0n = v0/np.sqrt(np.dot(v0,v0))
        v1n = v1/np.sqrt(np.dot(v1,v1))

        if np.cross(v0n,v1n) > 0:
            self.u = v0n
            self.v = v1n
            self.seg1 = seg
        else:
            self.u = v1n
            self.v = v0n
            self.seg1 = seg[:,::-1]

        self.dot = np.dot(self.u,self.v)
        self.cross = np.cross(self.u,self.v)


        if self.cross < 1e-15:
            self.degenerated=True
        self.upd_angle()

    def from2segs(self,seg0,seg1):
        """ creates a Cone from 2 segments

        Parameters
        ----------

        seg0 : 2 x 2  (Ndim x Npoints)
        seg1 : 2 x 2

        Notes
        -----

        The only way for the cone to be degenerated is when the two segments are on the same line.

        See Also
        --------

        pylayers.gis.layout.Layout.buildGi



        """
        # bv : (4,1)


        self.seg0 = seg0
        self.seg1 = seg1

        a0 = seg0[:,0]
        b0 = seg0[:,1]
        a1 = seg1[:,0]
        b1 = seg1[:,1]

        # check for connected segments (This could be determined earlier) 
        # a0 = a1 | b1
        # b0 = a1 | b1 

        # check segment orientation (crossing)
      
        if not (geu.ccw(a0,b0,b1) ^
                geu.ccw(b0,b1,a1) ):
            v0 = (b1 - a0)
            v1 = (a1 - b0)
            twisted = True
        else:    
            v0 = (a1 - a0)
            v1 = (b1 - b0)
            twisted = False
        
        v0n = v0/np.sqrt(np.dot(v0,v0))
        v1n = v1/np.sqrt(np.dot(v1,v1))

        if np.cross(v0n,v1n) > 0:
            self.u = v0n
            self.v = v1n
            inversion = False
        else:  
            self.u = v1n
            self.v = v0n
            inversion = True 
        
        if  (not twisted) & (not inversion) :
            #reverse seg1
            #print "reverse seg1"
            self.seg1 = self.seg1[:,::-1]
        if (inversion) & (not twisted):
            #reverse seg0 
            #print "reverse seg0"
            self.seg0 = self.seg0[:,::-1]
        if twisted & inversion:
            #reverse seg0 and seg1   
            #print "reverse seg0"
            #print "reverse seg1"
            self.seg0 = self.seg0[:,::-1]      
            self.seg1 = self.seg1[:,::-1]

        self.dot = np.dot(self.u,self.v)
        self.cross = np.cross(self.u,self.v)

        if self.cross < 1e-15:
            self.degenerated=True
        else:    
            a0u = np.dot(self.seg0[:,0],self.u)
            a0v = np.dot(self.seg0[:,0],self.v)
            b0u = np.dot(self.seg0[:,1],self.u)
            b0v = np.dot(self.seg0[:,1],self.v)

            kb  = ((b0v-a0v)-self.dot*(b0u-a0u))/(self.dot*self.dot-1)
            self.apex = self.seg0[:,1] + kb*self.v
        self.upd_angle()

    def from2csegs(self,seg0,seg1):
        """ creates a Cone from 2 connected segments

        Parameters
        ----------

        seg0 : 2 x 2  (Ndim x Npoints)
        seg1 : 2 x 2 

        Notes
        -----

        The only way for the cone to be degenerated is when the two segments are on the same line.


        """ 
        # bv : (4,1)


        self.seg0 = seg0
        self.seg1 = seg1

        a0 = seg0[:,0]
        b0 = seg0[:,1]
        a1 = seg1[:,0]
        b1 = seg1[:,1]
        
        # determine common point 
        if (np.dot(a0-a1,a0-a1)<1e-8):
            p = a0 
            u = b1-p
            v = p-b0
        elif (np.dot(a0-b1,a0-b1)<1e-8):
            p = a0 
            u = a1-p
            v = p-b0
            self.seg1 = self.seg1[:,::-1]
        elif (np.dot(b0-a1,b0-a1)<1e-8):
            p = b0 
            self.seg0 = self.seg0[:,::-1]
            u = b1-p
            v = p-a0
        elif (np.dot(b0-b1,b0-b1)<1e-8):
            self.seg0 = self.seg0[:,::-1]
            self.seg1 = self.seg1[:,::-1]
            p = b0 
            u = a1-p
            v = p-a0
        else:
            logging.critical('segment are not connected')
            pdb.set_trace()
       
        self.apex = p
        self.v = v/np.sqrt(np.dot(v,v))
        self.u = u/np.sqrt(np.dot(u,u))

        self.dot = np.dot(self.u,self.v)
        self.cross = np.cross(self.u,self.v)

        if self.cross < 1e-15:
            self.degenerated=True

        self.upd_angle()


    def show(self, **kwargs):
        """ show cone 

        Parameters
        ----------

        length : float 

        """

        defaults = {'length': 15.}
        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if 'seg1' not in self.__dict__:    
            verts = [tuple(self.apex),
                     tuple(self.apex + kwargs['length'] * self.u),
                     tuple(self.apex + kwargs['length'] * self.v),
                     tuple(self.apex)
                    ]
            codes = [Path.MOVETO,
                Path.LINETO,
                Path.LINETO,
                Path.CLOSEPOLY,
            ]

        else:
            a1 = self.seg1[:,0]
            b1 = self.seg1[:,1]
            if 'seg0' not in self.__dict__:    
                a0 = self.apex
                b0 = self.apex
            else:
                a0 = self.seg0[:,0]
                b0 = self.seg0[:,1]

            if not(self.degenerated):
                #verts = [tuple(self.apex),
                #         tuple(a1),
                #         tuple(b1),
                #         tuple(self.apex)
                #        ]
                verts = [tuple(self.apex),
                     tuple(self.apex + kwargs['length'] * self.u),
                     tuple(self.apex + kwargs['length'] * self.v),
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

        if 'fig' not in kwargs:
            fig = plt.figure(figsize=(10,10))
        else:
            fig = kwargs['fig']
        if 'ax' not in kwargs:    
            ax = fig.add_subplot(111)
        else:
            ax = kwargs['ax']

        ax.plot([self.apex[0],self.apex[0]+kwargs['length']*self.u[0]],
                [self.apex[1],self.apex[1]+kwargs['length']*self.u[1]],lw=1,color='b')
        ax.plot([self.apex[0],self.apex[0]+kwargs['length']*self.v[0]],
                [self.apex[1],self.apex[1]+kwargs['length']*self.v[1]],lw=1,color='r')
        if 'seg0' in self.__dict__:
            ax.plot([a0[0],b0[0]],[a0[1],b0[1]],lw=2,color='b')
        if 'seg1' in self.__dict__:
            ax.plot([a1[0],b1[0]],[a1[1],b1[1]],lw=2,color='r')
        patch = patches.PathPatch(path, facecolor='orange', lw=2, alpha=0.3)
        ax.add_patch(patch)
        ax.axis('equal')
        # ax.set_xlim(-2,2)
        # ax.set_ylim(-2,2)

        return(fig, ax)
