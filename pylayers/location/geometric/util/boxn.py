# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of RGPA.

#Foobar is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#Foobar is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
#authors :
#Nicolas AMIOT        : nicolas.amiot@univ-rennes1.fr
#Bernard UGUEN        : bernard.uguen@univ-rennes1.fr
#Mohamed LAARAIEDH    : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
from pylayers.util.project import *
import numpy as np
import os
import pdb
import copy
import time



#GeomNetType = np.dtype([('Id',np.uint64), 
#        ('time',np.uint64), 
#    ('p',float,(3)),
#    ('v',float,(3)),
#    ('a',float,(3))])


class LBoxN(PyLayers):
    """
    class LBoxN 
    List of BoxN 


    Atributes
    ---------


    box    : np.array
        array containing BoxN object. default void
    vol    : list
        volume of each boxes from self.box. default void
    bd    : np.arrays    2*len(box) x ndim
        contains all the boundaries of boxes from self.box. default void
    ndim    : int
        dimension of the nbox
    parmsh    : dictionnary
        keys ['display']    =True
             ['interactive']    =False

    LB : LBoxN
        Create a BoxN object from another one. default : None, the LBoxN obebject created is void.

    Methods
    -------

    mesure(self):
         measure intervals of box
    append(self,b):
        append a box 'b' to a Lbox
    append_l(self,lb):
        append a lbox 'lb' to a lbox
    info(self):
        display LBoxN information 
    bd2coord(self,Mapping = False):
        convert boundaries of Lbox to their vertexes coordintates
    octant(self):
        quadtree on Lboxes
    volume(self):
        estimate volume of LBoxes
    intersect(self,lb):
        EXPERIMENTAL find the intersection of LBoxN
    show3(self,col='b',Id=0):        
        required file generation for geomview display

    """
    #__slots__=('box','vol','bd','ndim','ctr','grav','parmsh')

    def __init__(self,Lb=None,ndim=3):
        self.ctr= []
        if Lb==None:
            self.box   = np.array([])
            self.vol   = []
            self.bd=[]
            self.ndim  = ndim   #   !!! TO BE DONE

        #    self.bnum = []


        else :

            self.box   = np.array([])
            self.vol   = []
            self.bd=[]
            for b in Lb:
                self.append(b)
                self.ndim=b.ndim

        self.mesure()
        self.parmsh={}
        self.parmsh['display']=True
        self.parmsh['interactive']=False


    def mesure(self):
        """ LMeasure BoxN

        Obtain measure of :
        - size of each interval from each dimension for each boxes
        - center of each interval from each dimension for each boxes
        -  Volume of the BoxN for each boxes (NOT WORKING)

        """

        if len(self.bd ) != 0:
            lbd=len(self.bd)/2
            self.ctr=np.zeros((lbd,self.ndim))
            for i in xrange(lbd):self.ctr[i,:]=(self.bd[2*i]+self.bd[(2*i)+1])/2.
#########################FONTIONNE MAIS TROP LOURD/TENT quand bcp de box
#            C  =np.array(([[1/2.,1/2.]]))
#            #M  =np.array(([[-1.,1.]]))    
#            I  = np.identity(len(self.bd)/2)
#            CTR = np.kron(I,C)
#            #MES = np.kron(I,M)
#            #self.mes  = np.dot(MES,self.bd)#self.bd[1,:]-self.bd[0,:]
#            self.ctr  = np.dot(CTR,self.bd)#(self.bd[1,:]+self.bd[0,:])/2.0
            #self.vol  = np.prod(self.mes,axis=1)
        else :
            self.ctr  = []

    def append(self,b):
        """append : Append a box to LboxN

        Parameters
        ----------
        b    : BoxN
            box to added

        Returns
        -------
        Nothing but update self.box status

        """
        self.box=np.append(self.box,b)


        try:
            self.bd=np.vstack((self.bd,b.bd[0]))
            self.bd=np.vstack((self.bd,b.bd[1]))
        except:
            self.bd=b.bd[0]
            self.bd=np.vstack((self.bd,b.bd[1]))

        V1  = sum(self.vol)
        V2  = b.vol
        # uptate center of gravity 
        try:
            self.grav = (V1*self.grav+V2*b.ctr)/(V1+V2)

        except:
            self.grav = b.ctr

        self.vol.append(b.vol)
        self.ctr.append(b.ctr)

    def append_l(self,lb):    
        """Append  LBoxN to LBoxN

        Parameters
        ----------

        lb    : LBoxN
            lbox to be added

        Returns
        -------

        Nothing but update self.box status

        """
#        for i in xrange(len(lb.box)):self.append(lb.box[i])
        self.box=np.append(self.box,lb.box)

        try :    

            self.bd=np.vstack((self.bd,lb.bd))    
            self.ctr=np.vstack((self.ctr,lb.ctr))
            self.vol=self.vol+lb.vol

        except:
            self.bd=lb.bd
            self.ctr=lb.ctr
            self.vol=lb.vol

#            try:
#                self.bd=np.vstack((self.bd,lb.bd[i][0]))
#                self.bd=np.vstack((self.bd,lb.bd[i][1]))
#            except:
#                self.bd=lb.bd[i][0]
#                self.bd=lb.bd[i][1]


#        self.box = self.box + lb.box
#        V1  = sum(self.vol)
#        V2  = lb.box[0].vol*len(lb.box)
#        # uptate center of gravity 
#        try:
#            self.grav = (V1*self.grav+V2*lb.grav)/(V1+V2)
#        except:
#            self.grav = lb.grav
#        self.vol = self.vol + lb.vol


    def info(self):
        """ display LBoxN information 
        """
        Vtot = 0 
        for k in range(len(self.box)):
            #print "Box : ",k," Volume :",self.vol[k]
            #print "ndim :",self.ndim
            print "------------"
            self.box[k].info()
            Vtot = Vtot+self.box[k].vol
        print "Volume : ",Vtot    

    def bd2coord(self,Mapping = False):
        """Boundaru to coordinates   

        Convert boundaries of Lbox to their vertexes coordintates

        in :
        [xmin ymin zmin]       
        [xmax ymax zmax]

        out :

        [xmin ymin zmin]
        [xmin ymax zmin]
        [xmax ymin zmin]
        [xmax ymax zmin]
        [xmin ymin zmax]
        [xmin ymax zmax]
        [xmax ymin zmax]
        [xmax ymax zmax]


        Parameters
        ----------

        Mapping : Boolean
            return a mapping of the vertex in regards of Lbox. Default = False

        Returns
        -------        

        P : array 2^ndim x ndim
            coordinates of box vertex 

        """

        lbd=len(self.bd)

        dimm1 = pow(2,self.ndim-1)
        dim = pow(2,self.ndim)

        P=np.zeros((lbd*dimm1,self.ndim)) #P=(len self.bd/2)*4
        # organisation de P
        if self.ndim == 3:
            R=np.repeat(range(0,dimm1*lbd,dimm1),2)
            R[range(1,len(R),2)]=R[range(1,len(R),2)]+1
        if self.ndim == 2:
            R=np.repeat(range(0,dimm1*lbd,dim),2)
            R[range(1,len(R),2)]=R[range(1,len(R),2)]+1

        if self.ndim == 3:
            RZ=np.repeat(range(0,dimm1*lbd,dim),dimm1)+(lbd/2)*range(0,dimm1,1)
        # aller chercher dans self.bd
        R2a=np.repeat(self.bd[range(0,lbd,2),:],dimm1,axis=0)
        R2b=np.repeat(self.bd[range(1,lbd,2),:],dimm1,axis=0)

#        # X
#        P[R,0]=R2a[:,0]#np.repeat(L.bd[range(0,lbd,2),0],4)
#        P[R+2,0]=R2b[:,0]#np.repeat(L.bd[range(1,lbd,2),0],4)
#        # Y
#        P[np.sort(np.mod(R+3,4*lbd)),1]=R2a[:,1]#np.repeat(L.bd[range(0,lbd,2),1],4)
#        P[R+1,1]=R2b[:,1]#np.repeat(L.bd[range(1,lbd,2),1],4)


        # X
        P[np.sort(np.mod(R+3,dimm1*lbd)),0]=R2a[:,0]#np.repeat(L.bd[range(0,lbd,2),1],4)
        P[R+1,0]=R2b[:,0]#np.repeat(L.bd[range(1,lbd,2),1],4)

        # Y
        P[R,1]=R2a[:,1]#np.repeat(L.bd[range(0,lbd,2),0],4)
        P[R+2,1]=R2b[:,1]#np.repeat(L.bd[range(1,lbd,2),0],4)


        if self.ndim == 3:
            # Z
            P[RZ,2]=R2a[:,2]#np.repeat(L.bd[range(0,lbd,2),2],4)
            P[RZ+4,2]=R2b[:,2]#np.repeat(L.bd[range(1,lbd,2),2],4)


        # mapping coresponding box
        if Mapping == True:
            Map=np.repeat(range(0,lbd/2,1),dim)
        #    Map=10*np.repeat(self.bnum,8)+range(0,8,1)*(len(self.bnum))
            return(P,Map)
#        P=np.array((
#            [self.bd[0,0],self.bd[0,1],self.bd[0,2]],
#            [self.bd[0,0],self.bd[1,1],self.bd[0,2]],
#            [self.bd[1,0],self.bd[1,1],self.bd[0,2]],
#            [self.bd[1,0],self.bd[0,1],self.bd[0,2]],
#            [self.bd[0,0],self.bd[0,1],self.bd[1,2]],
#            [self.bd[0,0],self.bd[1,1],self.bd[1,2]],
#            [self.bd[1,0],self.bd[1,1],self.bd[1,2]],
#            [self.bd[1,0],self.bd[0,1],self.bd[1,2]]
#               ))
        else:
            return(P)



    def octant(self):
        """ quadtree on boxes 

        Divide each Lboxes into 2^ndim equal parts
        aka Split each interval from each dimension into 2 equal part

        Returns
        -------

        lb    : LBoxN 2^ndim*self.box x ndim
            return theLBoxn from quadtree process
        """

        if self.ndim == 3:

            tlb = []
            lbox = LBoxN([],ndim=self.ndim) 
            C=np.array(([[1,1/2.,0],[0,1/2.,1]])).T
            I=np.identity(len(self.bd)/2)
            CC=np.kron(I,C)
            BD=np.dot(CC,self.bd)

            M=np.repeat(3*np.arange(0,len(self.bd)/2),16)  # groupe de 16 boundaries equivaut a 8 sous boites

            X = (range(0,2)+range(1,3))*4*(len(self.bd)/2)
            Y =(range(0,2)*2+range(1,3)*2)*2*(len(self.bd)/2)
            Z = (range(0,2)*4+range(1,3)*4)*(len(self.bd)/2)

            Rx = BD[X+M,0]
            Ry = BD[Y+M,1]
            Rz = BD[Z+M,2]
            lbd = np.array((Rx,Ry,Rz )).T
            llbd=len(lbd)
            xr=xrange(0,llbd,2)
            lb=LBoxN([BoxN(lbd[i:i+2],ndim=self.ndim) for i in xr])
            lb.mesure()



        if self.ndim == 2:

            tlb = []
            lbox = LBoxN(ndim=self.ndim) 
            C=np.array(([[1,1/2.,0],[0,1/2.,1]])).T
            I=np.identity(len(self.bd)/2)
            CC=np.kron(I,C)
            BD=np.dot(CC,self.bd)

            M=np.repeat(3*np.arange(0,len(self.bd)/2),8)  # groupe de 16 boundaries equivaut a 8 sous boites

            X = (range(0,2)+range(1,3))*2*(len(self.bd)/2)
            Y =(range(0,2)*2+range(1,3)*2)*(len(self.bd)/2)


            Rx = BD[X+M,0]
            Ry = BD[Y+M,1]

            lbd = np.array((Rx,Ry )).T
            llbd=len(lbd)


            lb=LBoxN(np.array([BoxN(lbd[i:i+2],ndim=self.ndim) for i in xrange(0,llbd,2)]))



        return(lb)


    def volume(self):
        """ Evaluate Boxes volume

        Compute the volume on each LBoxes
        """
        self.vol=[]
        for b in self.box:
            vol = b.vol
            if vol>=0:
                self.vol.append(vol)
            else:
                pdb.set_trace()



    def intersect(self,lb):
        """ EXPERIMENTAL

        Intersection of 2 LBOXN 

        """
        new_lb=LBoxN(ndim=self.ndim)
        for k in range(len(self.box)):
            for l in range(len(lb.box)):
                b=self.box[k].intersect(lb.box[l])
                if b.vol>0:   # si intersection non vide
                    new_lb.append(b)
        return(new_lb)        


    def show3(self,col='b',Id=0):
        """Show box into geomview

        generate a geomview file which allow to represent a box.


        Parameters
        ----------

        col    : string
            choose box color. compliant with matplotlib colorConverter. default 'r'
        Id    : list
            Identity of boxes to show.Default : [0]

        Returns
        -------

        filename    : string
            name of the generated file

        """
        filename  = "lbox"+str(Id)+".list"
        filename2 = basename +"/geom/"+filename
        fd = open(filename2,"w")
        fd.write("LIST\n")
        for k in range(len(self.box)):
            b = self.box[k]
            b.parmsh['display']=False
            filebox = b.show3(col=col,Id=[Id,k])
#            filebox = b.show3(col=col[k],Id=[Id,k])
            chaine = "{<"+filebox+"}\n"
            fd.write(chaine)
        #chaine = "{<cloud.list}\n"
        #fd.write(chaine)
        fd.close()

        if self.parmsh['display']:
            chaine = "geomview  -nopanel -b 1 1 1 " + filename2 + " 2>/dev/null &"
            os.system(chaine)
        return(filename)


class BoxN(PyLayers):
    """BoxN Class

    A box is determined by its boundary interval along each dimension



    Attributes
    ----------

    bd    : numpy array 2 x ndim  
        box boundary 

    ndim    : int
        dimension of the box (2D, 3D,...)

    self.parmsh : dictionnary
        display dictionnary for show 3 method TO BE CHANGE !!!
        keys    :['display']=True
            ['interactive']=False


    OBTAIN FROM mesure()
    self.mes    : array 1 x ndim 
        size of intervals of each dimension
    self.ctr    : array 1 x ndim 
        center of intervals of each dimension 
    self.vol    : float
        Volume of box


    Methods
    -------

    info()         : info about class 
    def mesure(self):
         measure intervals of box
    volume()       : evaluate volume 
    inbox(p)       : is p in box ?
    intersect(box) : intersection of two boxes               
    show3()        : geomview vizualization 
    cut()          : cut a box along given direction 



    TODO
    ----

    Remove parmsh and replace it by a ini file


    """
#    __slots__=('bd','ndim','mes','ctr','vol','parmsh')
    def __init__(self,bd=None,ndim=3):
#        if bd==None:
#            self.bd  = np.array([]).astype('float')
#        else:
#            for i in range(np.shape(bd)[1]):
#                assert bd[1,i]>=bd[0,i] , pdb.set_trace()
#                self.bd = bd.astype('float') 
        self.bd=bd
        self.ndim = ndim#np.shape(bd)[1]
        self.mesure()
        self.parmsh={}
        self.parmsh['display']=True
        self.parmsh['interactive']=False


#        print "%s from %s" % (inspect.stack()[1][3],inspect.stack()[1][1])


    def mesure(self):
        """ Measure BoxN
        Obtain measure of :
        - size of each interval from each dimension 
        - center of each interval from each dimension 
        - Volume of the BoxN

        """
        self.mes  = self.bd[1,:]-self.bd[0,:]
        self.ctr  = (self.bd[1,:]+self.bd[0,:])/2.0
        self.vol  = np.prod(self.mes)


    def setbd(self,vmin,vmax,axis=0):
        """
        setbd : set boundary value on axis
        """

        assert vmin<=vmax, "Incorrect bound"
        self.bd[0,axis]= vmin.astype('float')
        self.bd[1,axis]= vmax.astype('float')
        self.mesure()

    def void(self):
        """ return True if box is void
        """

        b = False
        if self.meas==[]:
            b = True
        else:
            pmes = np.prod(self.meas)
            if pmes==0:
                b = True

        return(b)

    def info(self):
        """ Information on BoxN
        """

        print "Volume (.vol) :",self.vol
        print "Center (.ctr)  :",self.ctr


    def bd2coord(self):
        """Boundary to coordinates

        Return an array containing of vertex from a box

        3D case :
        in :
        [xmin ymin zmin]
        [xmax ymax zmax]

        out :

        [xmin ymin zmin]
        [xmin ymax zmin]
        [xmax ymin zmin]
        [xmax ymax zmin]
        [xmin ymin zmax]
        [xmin ymax zmax]
        [xmax ymin zmax]
        [xmax ymax zmax]




        Returns
        -------

        P : array 2^ndim x ndim
            coordinates of box vertex 

        """    

        if self.ndim == 3:
            P=np.array(([self.bd[0,0],self.bd[0,1],self.bd[0,2]],
                [self.bd[0,0],self.bd[1,1],self.bd[0,2]],
                [self.bd[1,0],self.bd[1,1],self.bd[0,2]],
                [self.bd[1,0],self.bd[0,1],self.bd[0,2]],
                [self.bd[0,0],self.bd[0,1],self.bd[1,2]],
                [self.bd[0,0],self.bd[1,1],self.bd[1,2]],
                [self.bd[1,0],self.bd[1,1],self.bd[1,2]],
                [self.bd[1,0],self.bd[0,1],self.bd[1,2]]))
            return(P)

        if self.ndim == 2:
            P=np.array(([self.bd[0,0],self.bd[0,1]],
                [self.bd[0,0],self.bd[1,1]],
                [self.bd[1,0],self.bd[1,1]],
                [self.bd[1,0],self.bd[0,1]]))
            return(P)

    def coord2bd(self,coord):
        """

        Coordinates to boundary
        update boundary array from numpy array  of coordinates

        Parameters
        ----------

        coord : array 2^ndim x ndim
            vertexes coordinates of a boxN

        Returns
        -------

        Nothing but fills self.bd

        """
        self.bd[0,:]=np.min(coord,axis=0)
        self.bd[1,:]=np.max(coord,axis=0)


#    def octant(self):
#        tlb = []
#        lbox = LBoxN([])
#        BD = np.array((self.bd[0],(self.bd[0]+self.bd[1])/2,self.bd[1] ))
##        Rx = BD[(range(0,2)+range(1,3))*4,0]
##        Ry = BD[(range(0,2)*2+range(1,3)*2)*2,1]
##        Rz = BD[range(0,2)*4+range(1,3)*4,2]
##        O = np.array((Rx,Ry,Rz )).T
##        LB = LBoxN([BoxN(O[0:2]),BoxN(O[2:4]),BoxN(O[4:6]),BoxN(O[6:8]),BoxN(O[8:10]),BoxN(O[10:12]),BoxN(O[12:14]),BoxN(O[14:16])])
##    #    LB.bnum = range(00,010,1)
##        return(LB)
##        

    def octant(self):
        """ quadtree on boxes OBSOLETE

        Divide boxes into 2^ndim equal parts
        aka Split each interval from each dimension into 2 equal part

        """
        tlb = []
        lbox = LBoxN([])
        for k in range(self.ndim):
            tlb.append(self.cut(self.ctr[k],axis=k))

        lbm = tlb[0]
        for l in range(len(tlb)-1):
            lbp=lbm.intersect(tlb[l+1])
            lbm=lbp

        return(lbm)


    def intersect(self,box):
        """ Find intersection between current box and a given one

        Parameters
        ----------

        box    : BoxN
            a BoxN object

        Returns
        -------

        new_box    : BoxN
            a BoxN object

        """

        new_box = BoxN(np.zeros((2,self.ndim)),ndim=self.ndim)
        for k in range(self.ndim):    
            newmin = max(self.bd[0,k],box.bd[0,k])
            newmax = min(self.bd[1,k],box.bd[1,k])
            if (newmax>newmin):
                new_box.bd[0,k]= newmin
                new_box.bd[1,k]= newmax
        new_box.mesure()    
        return(new_box)    

    def bdiff(self,box):
        """ OBSOLETE

        USE self.intersect instead !!!

        """

        new_box = BoxN(np.zeros((2,self.ndim)),ndim=self.ndim)
        for k in range(self.ndim):    
            newmin = max(self.bd[0,k],box.bd[0,k])
            newmax = min(self.bd[1,k],box.bd[1,k])
            if (newmax>newmin):
                new_box.bd[0,k]= newmin
                new_box.bd[1,k]= newmax
        new_box.mesure()    
        return(new_box)    




    def show3(self,dim=(0,1,2),col='r',Id=[0],H_Id=0,alpha=0.2):
        """Show box into geomview

        generate a geomview file which allow to represent a box.


        Parameters
        ----------

        dim    : tuple
            chose dimension to display. default : (0,1,2)
        col    : string
            choose box color. compliant with matplotlib colorConverter. default 'r'
        Id    : list
            Identity of boxes to show.Default : [0]
        alpha   : float
            set transparency. Default 0.2

        Returns
        -------

        fname    : string
            name of the generated file

        """
        #for k in range(len(self.bb)):
        b        = np.zeros((6,4,3)).astype('int')
        b[0,:,:] = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]])
        b[1,:,:] = np.array([[1,0,0],[1,0,1],[1,1,1],[1,1,0]])
        b[2,:,:] = np.array([[0,0,0],[0,0,1],[0,1,1],[0,1,0]])
        b[3,:,:] = np.array([[0,1,0],[1,1,0],[1,1,1],[0,1,1]])
        b[4,:,:] = np.array([[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
        b[5,:,:] = np.array([[0,0,0],[0,0,1],[1,0,1],[1,0,0]])

        trans = -np.sort(-np.linspace(0.1,0.5,H_Id),-1)

        if self.parmsh['interactive']==True:
            print '(read geometry {define foo \n'
            print "LIST\n"
            print "{appearance  {-edge material {alpha "+str(alpha)+ " diffuse 1 0 0" +" }}"
            # 6 faces du cube
            for k in range(6):
                print "{QUAD "
                # 4 sommets de la face
                for li in range(4):
                    x = str(self.bd[b[k,i,0],dim[0]])
                    y = str(self.bd[b[k,i,1],dim[1]])
                    z = str(self.bd[b[k,i,2],dim[2]])
                    print x+" "+y+" "+z+" "
                print "} "
            print '}})'
        else:
            ch ='_'
            for c in Id:
                ch = ch+str(c)+'_'
            fname = "box"+ch+".list"
            filename = basename+"/geom/"+fname
            fd = open(filename,"w")
            fd.write("LIST\n")
            if self.ndim==3 :
                for k in range(6):
                    # 6 faces du cube
                    #filebbk = "face" +str(Id)+str(o)+str(l)+str(k)+".quad"
                    #fdk = open("./geom/"+filebbk,"w")
                    if col=='r':
                         col = " 1 0 0 "
                    elif col=='b':
                         col = " 0 0 1 "
                    elif col=='m':
                         col = " 1 0 1 "
                    elif col=='y':
                         col = " 1 1 0 "
                    elif col=='c':
                         col = " 0 1 1 "
                    elif col=='g':
                         col = " 0 1 0 "
                    elif col=='k':
                         col = " 0 0 0 "
                    elif col=='orange':
                        col = " 1 1 0.125 "
                    elif col=='skyblue':
                        col = " 0.91 1 1 "


                    fd.write("{appearance  {-edge material {alpha "+str(alpha)+" diffuse "+col +" }}")
                    fd.write("{ QUAD\n")
                    # 4 points de la face


                    for i in range(4):
                        x = str(self.bd[b[k,i,0],dim[0]])
                        y = str(self.bd[b[k,i,1],dim[1]])
                        z = str(self.bd[b[k,i,2],dim[2]])
                        fd.write(x+" "+y+" "+z+"\n")
#                    if self.ndim==2 :
#                        for i in range(4):
#                            x = str(self.bd[b[k,i,0],dim[0]])
#                            y = str(self.bd[b[k,i,1],dim[1]])
#                            z = str(5.0)
#                            fd.write(x+" "+y+" "+z+"\n")


                    fd.write(" }}\n")
                fd.close()

            if self.ndim==2 :

                Z=1.0*np.array(range(0,2)*len(self.bd))
                for k in range(6):        
                    # 6 faces du cube
                    #filebbk = "face" +str(Id)+str(o)+str(l)+str(k)+".quad"
                    #fdk = open("./geom/"+filebbk,"w")
                    if col=='r':
                         col = " 1 0 0 "
                    elif col=='b':
                         col = " 0 0 1 "
                    elif col=='m':
                         col = " 1 0 1 "
                    elif col=='y':
                         col = " 1 1 0 "
                    elif col=='c':
                         col = " 0 1 1 "
                    elif col=='g':
                         col = " 0 1 0 "
                    elif col=='k':
                         col = " 0 0 0 "
                    elif col=='orange':
                        col = " 1 0.59 0.125 "
                    elif col=='skyblue':
                        col = " 0.61 1 1 "


                    fd.write("{appearance  {-edge material {alpha "+str(alpha)+" diffuse "+col +" }}")
                    fd.write("{ QUAD\n")
                    # 4 points de la face


                    for i in range(4):
                        x = str(self.bd[b[k,i,0],dim[0]])
                        y = str(self.bd[b[k,i,1],dim[1]])
#                        z = str(self.bd[b[k,i,2],dim[2]])
                        z = str(Z[b[k,i,2]])
                        fd.write(x+" "+y+" "+z+"\n")
#                    if self.ndim==2 :
#                        for i in range(4):
#                            x = str(self.bd[b[k,i,0],dim[0]])
#                            y = str(self.bd[b[k,i,1],dim[1]])
#                            z = str(5.0)
#                            fd.write(x+" "+y+" "+z+"\n")


                    fd.write(" }}\n")
                fd.close()


            if self.parmsh['display']:
                chaine = "geomview  -nopanel -b 1 1 1 " + filename + " 2>/dev/null &"
                os.system(chaine)

            return(fname)

