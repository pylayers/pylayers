import pdb
import numpy as np
from xml.etree.ElementTree import Element
from xml.etree import ElementTree

class Point(Element):
    def __init__(self,name,p):
        Element.__init__(self,name,X=str(p[0]),Y=str(p[1]),Z=str(p[2]))

class Vertex(Element):
    def __init__(self,x=0,y=0,z=[]):
        Element.__init__(self,'Vertex')
        if z==[]:
            #2D
            self.attrib['X1'] = str(x)
            self.attrib['X2'] = str(y)
        else:
            #3D
            self.text = str(x)+','+str(y)+','+str(z)
        pass

class Face(Element):
    def __init__(self, x=0, y=0, z=0):
        Element.__init__(self, 'Face')
        self.text = str(x)+','+str(y)+','+str(z)

class Cylinder(Element):
    def __init__(self,P1,P2,Radius=50,Priority=10):
        Element.__init__(self,'Cylinder')
        self['Priority'] = Priority
        self['P1'] = Point(P1)
        self['P2'] = Point(P2)
        self['Radius']  = Radius

class Sphere(Element):
    """ a Sphere is a Primitives of a Properties element

    Use set  of CSX object
    """
    def __init__(self,P,R=50,Pr=10):
        Element.__init__(self,'Sphere',Priority=str(Pr),Radius=str(R))
        self.append(Point('Center',P))

class Box(Element):
    def __init__(self,P1,P2,Pr):
        Element.__init__(self,'Box',Priority=str(Pr))
        self.append(Point('P1',P1))
        self.append(Point('P2',P2))

class Polygon(Element):
    def __init__(self,LP,Priority=1,elevation=254,normdir=2,coordsystem=0):
        """
        Parameters
        ----------

        Priority
        Elevation

        """
        Element.__init__(self,'Polygon',
                         Priority=str(Pr),
                         Elevation=str(elevation),
                         NormDir=str(normdir),
                         CoordSystem=str(coordsystem))
        for P in LP:
            V = Vertex(x=P[0],y=P[1])
            self.append(V)

class ConductingSheet(Element):
  pass

class Polyhedron(Element):
    """
        AddPolyhedron.m
    """
    def __init__(self,LP= [[0,0,0],
                          [1,0,0],
                          [0,1,0],
                          [0,0,1]]
                     ,LV= [[0,1,2],
                           [1,2,3],
                           [0,1,3],
                           [0,2,3]],
                      Pr=0):
        """
        Parameters
        ----------

        LP : List of Points
        LV : List of vertex indexes
        Pr : Priority

        """
        Element.__init__(self,'Polyhedron',Priority=str(Pr))
        for P in LP:
            self.append(Vertex(x=P[0],y=P[1],z=P[2]))
            self.append(Face(x=P[0], y=P[1],z=P[2]))

class Curve(Element):
   def __init__(self,point ,Pr):
        Element.__init__(self,'Curve',Priority=str(Pr))
        for P in LP:
            V = Vertex(x=P[0],y=P[1],z=P[2])
            self.append(V)

class XLines(Element):
    def __init__(self,X=np.arange(-10,11,1)):
        Element.__init__(self,'XLines')
        c=reduce(lambda a,b: str(a)+','+str(b),X)
        self.text = c

class YLines(Element):
    def __init__(self,X=np.arange(-10,11,1)):
        Element.__init__(self,'YLines')
        c=reduce(lambda a,b: str(a)+','+str(b),X)
        self.text = c


class ZLines(Element):
    def __init__(self,X=np.arange(-10,31,1)):
        Element.__init__(self,'ZLines')
        c=reduce(lambda a,b: str(a)+','+str(b),X)
        self.text = c

class RectilinearGrid(Element):
    def __init__(self,X,Y,Z,CoordSystem=0,DeltaUnit="1"):
        """
        Parameters
        ----------

        X  : np.array
        Y  : np.array
        Z  : np.array
        CoordSystem : int
            default 0 : cartesian
        DeltaUnit : string
            default "1"

        """
        if CoordSystem==0:
            Element.__init__(self,'RectilinearGrid')
            self.attrib['DeltaUnit']=DeltaUnit
            self.attrib['CoordSystem']='0'
            self.append(XLines(X))
            self.append(YLines(Y))
            self.append(ZLines(Z))

class LinPoly(Element):
    """ Planar surface defined by a polygon

     Parameters
     ----------

     prio:      priority
     normDir:   normal direction of the polygon,
                e.g. 'x', 'y' or 'z', or numeric 0..2
     elevation: position of the polygon plane
     points:    list of points (2xN) 
     length:    linear extrution in normal direction, starting at elevation

    """
    def __init__(self,**kwargs):
        defaults = {'elevation' : -1,
                    'pr':10,
                    'lp':np.array([[0,10],
                                 [50,27.1010071662834],
                                 [50,-27.1010071662834],
                                 [0,-10]]),
                    'normdir':1,
                    'length':2}
        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        Element.__init__(self,'LinPoly', Priority=str(kwargs['pr']),Elevation=str(kwargs['elevation']),Length=str(kwargs['length']),NormDir=str(kwargs['normdir']))
        for p in kwargs['lp']:
            self.append(Vertex(x=p[0],y=p[1],z=[]))

class Transformation(Element):
    """docstring for Transformation"""
    def __init__(self):
        Element.__init__(self,'Transformation')

class Translate(Element):
    def __init__(self,p):
        Element.__init__(self,'Translate',Argument=str(p[0])+","+str(p[1])+","+str(p[2]))

class Rotate_X(Element):
    def __init__(self,ang):
        Element.__init__(self,'Rotate_X',Argument=str(ang))

class Rotate_Y(Element):
    def __init__(self,ang):
        Element.__init__(self,'Rotate_Y',Argument=str(ang))

class Rotate_Z(Element):
    def __init__(self,ang):
        Element.__init__(self,'Rotate_Z',Argument=str(ang))

def SmoothMeshLine( plines, max_res, ratio=1.3,check=True,max_ratio=1.25):
	"""
	Parameters
	----------

	plines : np.array

	Returns
	-------

	lines : np.array

	"""

	dlines = np.diff(plines)
	Npoints = np.ceil(dlines/max_res)
	Npoints = Npoints.astype(int)

	for k,N in enumerate(Npoints):
		if k!=len(dlines)-1:
			l = np.linspace(plines[k],plines[k+1],N,endpoint=False)
		else:
			l = np.linspace(plines[k],plines[k+1],N+1,endpoint=True)

		try:
			lines = np.hstack((lines,l))
		except:
			lines = l
	#max_ratio = ratio*allowed_max_ratio


	if check:
		EC,pos,Etype=CheckMesh(lines,0,max_res,max_ratio,0);

	return lines

def CheckMesh(lines, min_res, max_res, ratio, verbose=False):
	""" Check if mesh lines are valid

   Parameters
   ----------

   lines   : np.array()
   min_res: minimal allowed mesh-diff
   max_res: maximal allowed mesh-diff
   ratio:   maximal allowed mesh-diff ratio
   be_quiet: disable warnings

   Returns
   -------

   EC:     error code (number of found errors)
   pos:    line positions with error
   E_type: error type

   From Matlab code of Thorsten Liebig

   See Also
   --------

	"""


	diff_lines = np.diff(lines)
	EC = 0
	E_type = []

	pos = []
	max_err = np.where(diff_lines>max_res)[0]

	if (len(max_err) >0 & verbose):
		print('CheckMesh : found resolution larger than max_res')

		pos = np.hstack((pos,max_err))
		EC = EC + len(max_err)
		E_type.append(1)

	min_err = np.where(diff_lines<min_res)[0]
	if (len(min_err)>0 & verbose):
		warning('CheckMesh : found resolution smaller than min_res')
		pos = np.hstack((pos,min_err))
		EC = EC + len(min_err)
		E_type.append(2)

	r = diff_lines[1:]/diff_lines[0:-1]

	if (r>ratio*1.01).any():
		u = np.where(r>ratio*1.01)
		pos = np.hstack((pos,u))
		E_type.append(3)

	if (r<(1/ratio)*1.01).any():
		u = np.where(r<(1/ratio)*1.01)
		pos = np.hstack((pos,u))
		E_type.append(4)


	return EC,pos,E_type