from xml.etree import ElementTree
from xml.etree.ElementTree import Element
from xml.etree.ElementTree import SubElement
import numpy as np
import scipy as sp
import io
import pdb
import os
from xml.etree.ElementTree import Element
from xml.etree import ElementTree
# temporary import
from geometry import *
from matter import *
from excitation import *
from probe import *
from dumpbox import *


class OpenEMS(Element):
    """
    Main Class decribing the openEMS simulation
    """
    def __init__(self,FDTD,CSX):
        Element.__init__(self,'openEMS')
        self.append(FDTD)
        self.append(CSX)


    def __repr__(self):
        st = ElementTree.tostring(self)
        return(st)

    def save(self,filename='openEMS.xml'):
        """ save the whole openems xml configuration file

        Parameters
        ----------

        filename : string
            default : "openEMS.xml"

        Notes
        -----
        We here use xmlllint script for pretty formatting

        """
        self.filename = filename
        output_file = open(filename, 'w' )
        output_file.write( '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>' )
        output_file.write(ElementTree.tostring( self ))
        output_file.close()
        os.system("xmllint --format "+filename+"> tmpf")
        os.system("mv tmpf "+filename)


    def geomplot(self):
        os.system('~/Apps/openEMS/bin/AppCSXCAD '+self.filename)

    def run(self):
        os.system('~/Apps/openEMS/bin/openEMS.sh '+self.filename)

class CSX(Element):
    """ Continuous Structure Class

    Methods
    -------

    add
    set
    save

    """
    def __init__(self,CoordSystem=0):
        if CoordSystem==0:
            Element.__init__(self,'ContinuousStructure',CoordSystem=str(CoordSystem))
            # CSX has always a Properties Section
            P = Element('Properties')
            self.append(P)

    def save(self,filename='CSX.xml'):
        output_file = open( filename, 'w' )
        output_file.write( '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>' )
        output_file.write( ElementTree.tostring( self ) )
        output_file.close()


    def primitive(self,name,elp):
        """ append a primitive to a property name

        Parameters
        ----------

        name : property name to be set
        elp  : element primitive

        """
        P = self.find('Properties')
        for ch in P.getchildren():
            if ch.attrib['Name']==name:
                Prim = ch.find('Primitives')
                Prim.append(elp)

    def add(self, el, p=[], **kwargs):
        """ add property:name

        Parameters
        ----------

        el : element

        All properties should have a distinct name

        Examples
        --------

        >>> C = CSX()
        >>> Me = Element('Metal',Name='gnd')
        >>> Ma = Element('Material',Name='substrate')
        >>> b1 = Box(P1=[0,0,0], P2=[100,100,200],Pr=0)
        >>> b2 = Box(P1=[0,0,0], P2=[100,100,200],Pr=0)
        >>> C1 = Cylinder(P1=[0,0,0], P2=[0,0,200],P=10,R=10)
        >>> M = RectilinearGrid(L1,L2,L3,DeltaUnit=0.001,CoordSyst="0")
        >>> C.add(Me,b1,offset=)
        >>> C.add(Ma,b2,offset=)
        >>> C.add(M)

        """
        #if 'name' not in self.proper
        if el.attrib.has_key('Name'):
            nameExist = False
            P = self.find('Properties')
            for ch in P.getchildren():
                if el.attrib['Name']==ch.attrib['Name']:
                    nameExist = True
                    raise NameError(el.tag + ' already exists')
            # if no error raised
            P.append(el)
            if p!=[]:
                self.primitive(el.attrib['Name'],p)
        else:
            # el is not a Property
            for ch in self.getchildren():
                if el.tag==ch.tag:
                    raise AttributeError('Element already exists')
            # if no error raised
            self.append(el)
            if p!=[]:
                self.primitive(el.attrib['Name'],p)

class FDTD(Element):
    """

     Inititalize the FDTD data-structure.

     optional field arguments for usage with openEMS:

       NrTS:           max. number of timesteps to simulate (e.g. default=1e9)
       EndCriteria:    end criteria, e.g. 1e-5, simulations stops if energy has
                       decayed by this value (<1e-4 is recommended, default=1e-5)
       MaxTime:        max. real time in seconds to simulate
       OverSampling:   nyquist oversampling of time domain dumps
       CoordSystem:    choose coordinate system (0 Cartesian, 1 Cylindrical)
       MultiGrid:      define a cylindrical sub-grid radius
       TimeStep:       force to use a given timestep (dangerous!)
       TimeStepFactor: reduce the timestep by a given factor (>0 to <=1)
       TimeStepMethod: 1 or 3 chose timestep method (1=CFL, 3=Rennigs (default))
       CellConstantMaterial: set to 1 to assume a material is constant inside
                             a cell (material probing in cell center)
       f_max

    Examples
    --------

    default init with 1e9 max. timesteps and -50dB end-criteria
        >>> F = FDTD();

    init with 1e6 max. timesteps and -60dB end-criteria
        >>> F = FDTD(NrTs=1e6, EndCriteria=1e-6);

    cylindrical FDTD simulation
        >>> F = FDTD('CoordSystem', 1);

    See also
    --------


    openEMS matlab interface
    --------------------------

    author: Thorsten Liebig (c) 2010-2013
    python version : B.Uguen, Meriem Asghar, Amath Ndiaye, Adil AIT SOUDANE

    """

    def __init__(self,**kwargs):
        """
        
        Parameters
        ----------

        'NumberOfTimestepsrTs' : 1e9,
        'endCriteria': 1e-5,
        'MaxTime': 0 ,
        'OverSampling': 0,
        'CoordSystem':0,
        'MultiGrid':0,
        'f_max' : 

        """

        Element.__init__(self,'FDTD')
        d = {}
        for k in kwargs:
            if type(kwargs[k])!=str:
                d[k]=str(kwargs[k])
            else:
                d[k]=kwargs[k]

        self.attrib.update(d)

    def __repr__(self):
        return(ElementTree.tostring( self ) )

    def add(self,el):
        self.append(el)

    def save(self,filename='FDTD.xml'):
        output_file = open( filename, 'w' )
        output_file.write( '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>' )
        output_file.write( ElementTree.tostring( self ) )
        output_file.close()

class Exc(Element):
    """

    Parameters
    ----------

    typ : Gaussian (0)
            f0 : center frequency
            fc :
          Sinus (1)
            f0 : frequency
          Dirac (2)
          Step (3)
          Custom (10)
            f0 : nyquist rate
            funcStr : string describing the excitation function e(t)
    """
    def __init__(self,typ='Gaussian',**kwargs):
        dtrans={'Gaussian':0,
                'Sinus':1,
                'Dirac':2,
                'Step':3,
                'Custom':10}
        Element.__init__(self,'Excitation',Type=str(dtrans[typ]))
        if typ=='Gaussian':
            self.attrib['f0']=str(kwargs['f0'])
            self.attrib['fc']=str(kwargs['fc'])
            self.attrib['f_max']=str(kwargs['f0']+kwargs['fc'])

        if typ=='Sinus':
            self.attrib['f0']=str(kwargs['f0'])
            if 'f_max' in kwargs:
              self.attrib['f_max']=str(kwargs['f_max'])

        if typ=='Custom':
            self.attrib['f0']=str(kwargs['f0'])
            self.attrib['Function']=str(kwargs['funcStr'])

class BoundaryCond(Element):
  """

   BC = [xmin xmax ymin ymax zmin zmax]

   0 = PEC      or  'PEC'
   1 = PMC      or  'PMC'
   2 = MUR-ABC  or  'MUR'
   3 = PML-ABC  or  'PML_x' with pml size x => 4..50

   Examples
   --------

   BC = [ 1 , 1  , 0  , 0  , 2   , 3     ]
   BC = ['PMC' 'PMC' 'PEC' 'PEC' 'MUR' 'PML_8'}

  mur-abc definitions
  define a phase-velocity to be used by the mur-abc
  useful e.g. for dispersive waveguides
  FDTD = SetBoundaryCond(FDTD,BC,'MUR_PhaseVelocity',299792457.93272);


 pml definitions
   arguments:  'PML_Grading','gradFunction'
     Define the pml grading grading function.
     Predefined variables in this grading function are:
       D  = depth in the pml in meter
       dl = mesh delta inside the pml in meter
       W  = width (length) of the pml in meter
       N  = number of cells for the pml
       Z  = wave impedance at the current depth and position

  """
  def __init__(self,BC=['MUR','MUR','MUR','MUR','MUR','MUR']):
    Element.__init__(self,'BoundaryCond')
    if type(BC[0])==str:
      self.attrib['xmin']=BC[0]
      self.attrib['xmax']=BC[1]
      self.attrib['ymin']=BC[2]
      self.attrib['ymax']=BC[3]
      self.attrib['zmin']=BC[4]
      self.attrib['zmax']=BC[5]
    else:
      dtrans={0:"PEC",1:"PMC",2:"MUR",3:"PML"}
      self.attrib['xmin']=dtrans[BC[0]]
      self.attrib['xmax']=dtrans[BC[1]]
      self.attrib['ymin']=dtrans[BC[2]]
      self.attrib['ymax']=dtrans[BC[3]]
      self.attrib['zmin']=dtrans[BC[4]]
      self.attrib['zmax']=dtrans[BC[5]]


