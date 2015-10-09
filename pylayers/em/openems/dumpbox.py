import numpy as np
import scipy as sp
import pdb

from xml.etree import ElementTree
from xml.etree.ElementTree import Element
from xml.etree.ElementTree import SubElement
from geometry import *

class DumpBox(Element):
    """
    Add a dump property to CSX with the given name.

   Frequency:  specify a frequency vector (required for dump types >=10)

   SubSampling:   field domain sub-sampling, e.g. '2,2,4'
   OptResolution: field domain dump resolution, e.g. '10' or '10,20,5'

   MultiGridLevel: Request to dump from a multigrid level (default is 0)
                   Note: This only takes effect if the method supports and
                   uses multiple grids!

   Warning:
       FDTD Interpolation abnormalities:
         - no-interpolation: fields are located in the mesh by the
           Yee-scheme, the mesh only specifies E- or H-Yee-nodes 
             --> use node- or cell-interpolation or be aware of the offset
         - E-field dump & node-interpolation: normal electric fields on
           boundaries will have false amplitude due to forward/backward
           interpolation  in case of (strong) changes in material
           permittivity or on metal surfaces
             --> use no- or cell-interpolation
         - H-field dump & cell-interpolation: normal magnetic fields on
           boundaries will have false amplitude due to forward/backward
           interpolation in case of (strong) changes in material permeability
             --> use no- or node-interpolation

 e.g. AddDump(CSX,'Et');
      CSX = AddBox(CSX,'Et',10,[0 0 0],[100 100 200]); %assign box

 or   AddDump(CSX,'Ef',DumpType, 10, 'Frequency',[1e9 2e9]);
      CSX = AddBox(CSX,'Ef',10,[0 0 0],[100 100 200]); %assign box

 or   AddDump(CSX,'Ht','SubSampling','2,2,4','DumpType',1);
      CSX = AddBox(CSX,'Ht',10,[0 0 0],[100 100 200]); %assign box

     See also AddMaterial, AddExcitation, AddProbe, AddMetal, AddBox

    """
    def __init__(self,name="Et",DumpMode='nointerp',save='vtk'):
        """

        Parameters
        ----------

        name : string
            name of the DumpBox (def = 'Et')

        DumpMode : string
            'none' : No Interpolation
            'node' : Node Interpolation
            'cell' : Cell Interpolation

        DumpType:  string
                'et' :  E-field time-domain dump (default)  (0)
                'ht' :  H-field time-domain dump (1)
                'jt' :  electric current time-domain (2)
                'tt' :  total current density (rot(H)) (3)

                 'ef' : E-field frequency-domain (10)
                 'hf' : H-field frequency-domain (11)
                 'jf' : for electric current frequency-domain (12)
                 'tf' : total current density (rot(H)) (13)

                 'sar' :  local SAR  (20)
                 'sarg' :  1g averaging SAR (21)
                 'sarggg' : 10g averaging SAR (22)

                 'raw' :  raw data needed for SAR calculations
                 (electric field FD,cell volume, conductivity and density) (29)

        save : string
                'vtk'  : 0
                'hdf5' : 1

        """

        dumptype = {'et':0,'ht':1,'jt':2,'tt':2,'ef':10,'hf':11,'jf':12,'tf':13,'sar':20,'sarg':21,'sargg':22,'raw':29}
        savetype = {'vtk':0,'hdf5':1}

        Element.__init__(self,'DumpBox')

        self.attrib['Name']  =name
        if DumpMode=='none':
            self.attrib['DumpMode']="0"
        if DumpMode=='node':
            self.attrib['DumpMode']="1"
        if DumpMode=='cell':
            self.attrib['DumpMode']="2"
        Prim = Element('Primitives')
        self.append(Prim)

