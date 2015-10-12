import numpy as np
import scipy as sp
import pdb

from xml.etree import ElementTree
from xml.etree.ElementTree import Element
from xml.etree.ElementTree import SubElement
from geometry import *
"""

Creates an E-field or H-field excitation.

name: property name for the excitation
type:

    Es=E-field soft excitation
    E1=E-field hard excitation
    2=H-field soft excitation
    3=H-field hard excitation
    10=plane wave excitation

excite: e.g. [2 0 0] for excitation of 2 V/m in x-direction

additional options for openEMS:
 'Delay'  : setup an excitation time delay in seconds
 'PropDir': direction of plane wave propagation (plane wave excite only)

example:

 CSX = AddExcitation( CSX, 'infDipole', 1, [1 0 0] );
start = [-dipole_length/2, 0, 0];
stop  = [+dipole_length/2, 0, 0];
CSX = AddBox( CSX, 'infDipole', 1, start, stop );

CSXCAD matlab interface
-----------------------

author: Thorsten Liebig

See also SetExcitationWeight, AddMetal, AddExcitation, AddProbe,
AddDump, AddBox

CSX = AddProperty(CSX, 'Excitation', name, 'Type', type, 'Excite', excite, varargin{:});

"""

class Excitation(Element):
    def __init__(self,name="exc0",typ='Es', excite="1,0,0"):
        Element.__init__(self,'Excitation')
        self.attrib['Name']  =name
        self.attrib['Excite']=excite
        if typ=='Es':
            self.attrib['Type']="0"
        if typ=='Eh':
            self.attrib['Type']="1"
        if typ=='Hs':
            self.attrib['Type']="2"
        if typ=='Hh':
            self.attrib['Type']="3"
        if typ=='PW':
            self.attrib['Type']="10"

        Prim = Element('Primitives')

        self.append(Prim)



#def AddRectWaveGuidePort( prio, portnr, start, stop, dir, a, b, mode_name, exc_amp, varargin ):
#    """
# function [CSX,port] = AddRectWaveGuidePort( CSX, prio, portnr, start, stop, dir, a, b, mode_name, exc_amp, varargin )
#
# Create a rectangular waveguide port, including an optional excitation and probes
#
# Note: - The excitation will be located at the start position in the given direction
#       - The voltage and current probes at the stop position in the given direction
#
# input:
#   CSX:        complete CSX structure (must contain a mesh)
#   prio:       priority of primitives
#   start:      start coordinates of waveguide port box
#   stop:       stop  coordinates of waveguide port box
#   dir:        direction of port (0/1/2 or 'x'/'y'/'z'-direction)
#   a,b:        rectangular waveguide width and height (in meter)
#   mode_name:  mode name, e.g. 'TE11' or 'TM21'
#   exc_amp:    excitation amplitude (set 0 to be passive)
#
# optional (key/values):
#   varargin:   optional additional excitations options, see also AddExcitation
#   'PortNamePrefix': a prefix to the port name
#
# output:
#   CSX:        modified CSX structure
#   port:       port structure to use with calcPort
#
# example:
#   % create a TE10 circular waveguide mode, using cylindircal coordinates
#   start=[mesh.r(1)   mesh.a(1)   0  ];
#   stop =[mesh.r(end) mesh.a(end) 100];
#   [CSX,port] = AddCircWaveGuidePort( CSX, 99, 1, start, stop, 320e-3, 'TE11', 0, 1);
#
# openEMS matlab interface
# -----------------------
# (c) 2013 Thorsten Liebig (thorsten.liebig@gmx.de)
#
# See also InitCSX, AddExcitation, calcWGPort, calcPort
# """
#    if mode_name[0:2]!='TE':
#        print  'currently only TE type modes are supported'
#        exit()
#
#    m = mode_name[2]
#    n = mode_name[3]
#
#    # values by David M. Pozar, Microwave Engineering, third edition
#    # Pozar 3.83 page 123
#    #
#    kc = np.sqrt((m*np.pi/a)**2 + (n*np.pi/b)**2)
#
#
#    unit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit
#
#    dir = DirChar2Int(dir);
#    dir_names={'x','y','z'};
#
#    dirP = mod((dir+1),3)+1;
#    dirPP = mod((dir+2),3)+1;
#    nameX = ['(' dir_names{dirP}  '-' num2str(start(dirP)) ')'];
#    nameY = ['(' dir_names{dirPP} '-' num2str(start(dirPP)) ')'];
#
#    # convert a&b to drawing units
#    a = a/unit;
#    b = b/unit;
#
#    # functions by David M. Pozar, Microwave Engineering, second edition
#    # electric field mode profile
#    # 3.82 a and 3.82 b page 122
#
#    func_Ex = str( n/b) +'*cos('+str(m*pi/a) +'*'+ nameX + ')*sin(' + str(n*pi/b)+ '*' + nameY + ')'
#    func_Ey = str(-m/a) +'*sin('+str(m*pi/a) +'*'+ nameX + ')*cos(' + str(n*pi/b)+ '*' + nameY + ')'
#
#    # magnetic field mode profile
#
#    func_Hx = str(m/a) +'*sin(' +str(m*pi/a) + '*' nameX + ')*cos(' + str(n*pi/b) + '*' nameY ')'
#    func_Hy = str(n/b) +'*cos(' +str(m*pi/a) + '*' nameX + ')*sin(' + str(n*pi/b) + '*' nameY ')'
#
#
#    func_E{dir+1} = 0;
#    func_E{dirP} = func_Ex;
#    func_E{dirPP} = func_Ey;
#
#    func_H{dir+1} = 0;
#    func_H{dirP} = func_Hx;
#    func_H{dirPP} = func_Hy;
#
#[CSX,port] = AddWaveGuidePort( CSX, prio, portnr, start, stop, dir, func_E, func_H, kc, exc_amp, varargin{:} );

#
#def  AddWaveGuidePort( CSX, prio, portnr, start, stop, dir, E_WG_func, H_WG_func, kc, exc_amp, varargin )
#    """
# function [CSX,port] = AddWaveGuidePort( prio, portnr, start, stop, dir, E_WG_func, H_WG_func, kc, exc_amp, varargin )
#
# Create a waveguide port, including an optional excitation and probes
#
# Note: - The excitation will be located at the start position in the given direction
#       - The voltage and current probes at the stop position in the given direction
#
# parameter:
#   CSX:        complete CSX structure (must contain a mesh)
#   prio:       priority of primitives
#   start:      start coordinates of waveguide port box
#   stop:       stop  coordinates of waveguide port box
#   dir:        direction of port (0/1/2 or 'x'/'y'/'z'-direction)
#   E_WG_func:  electric field mode profile function as a string
#   H_WG_func:  magnetic field mode profile function as a string
#   kc:         cutoff wavenumber (defined by the waveguide dimensions)
#   exc_amp:    excitation amplitude (set 0 to be passive)
#
# optional (key/values):
#   varargin:   optional additional excitations options, see also AddExcitation
#   'PortNamePrefix': a prefix to the port name
#
# output:
#       CSX:        modified CSX structure
#       port:       port structure to use with calcPort
#
#    Examples:
#    ---------
#       % create a TE11 circular waveguide mode, using cylindircal coordinates
#       >>> p11 = 1.841;
#       >>> kc = p11 / radius;  % cutoff wavenumber with radius in meter
#       >>> kc_draw = kc*unit;  % cutoff wavenumber in drawing units
#
#       # electric field mode profile
#
#       >>> func_E{1} = [ num2str(-1/kc_draw^2,15) '/rho*cos(a)*j1('  num2str(kc_draw,15) '*rho)'];
#       .func_E{2} = [ num2str(1/kc_draw,15) '*sin(a)*0.5*(j0('  num2str(kc_draw,15) '*rho)-jn(2,'  num2str(kc_draw,15) '*rho))'];
#   func_E{3} = 0;
#
#   % magnetic field mode profile
#   func_H{1} = [ '-1*' num2str(1/kc_draw,15) '*sin(a)*0.5*(j0('  num2str(kc_draw,15) '*rho)-jn(2,'  num2str(kc_draw,15) '*rho))'];
#   func_H{2} = [ num2str(-1/kc_draw^2,15) '/rho*cos(a)*j1('  num2str(kc_draw,15) '*rho)'];
#   func_H{3} = 0;
#
#   start=[mesh.r(1)   mesh.a(1)   0  ];
#   stop =[mesh.r(end) mesh.a(end) 100];
#   [CSX, port{1}] = AddWaveGuidePort(CSX, 0, 1, start, stop, 2, func_E, func_H, kc, 1);
#
# openEMS matlab interface
# -----------------------
# (c) 2013 Thorsten Liebig (thorsten.liebig@gmx.de)
#    adapted in python B.Uguen
#
# See also InitCSX, AddExcitation, calcWGPort, calcPort
#
#"""
#
#    dir = DirChar2Int(dir);
#
#    port.type='WaveGuide';
#    port.nr=portnr;
#    port.kc = kc;
#    port.dir = dir;
#    port.drawingunit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;
#
#    PortNamePrefix = '';
#
#    varargin_tmp  = varargin;
#    for n=1:2:numel(varargin_tmp):
#        if strcmpi('PortNamePrefix',varargin_tmp{n}):
#        PortNamePrefix = varargin_tmp{n+1};
#        varargin([n n+1]) = [];
#
## matlab adressing
#    dir = dir + 1;
#    dir_sign = sign(stop(dir) - start(dir));
#    if (dir_sign==0):
#        dir_sign = 1
#
#    port.direction = dir_sign;
#
#E_WG_func{dir} = 0;
#H_WG_func{dir} = 0;
#
#port.excite = 0;
#if (exc_amp~=0):
#    if (start(dir)==stop(dir)):
#        error 'if waveguide port is to be excited, the length in propagation direction must not be zero'
#    e_start = start;
#    e_stop = stop;
#    e_stop(dir) = e_start(dir);
#    port.excite = 1;
#    port.excitepos = e_start(dir);
#    e_vec = [1 1 1]*exc_amp;
#    e_vec(dir) = 0;
#    exc_name = [PortNamePrefix 'port_excite_' num2str(portnr)];
#
#    CSX = AddExcitation( CSX, exc_name, 0, e_vec, varargin{:});
#    CSX = SetExcitationWeight(CSX, exc_name, E_WG_func );
#	CSX = AddBox( CSX, exc_name, prio, e_start, e_stop);
#
## voltage/current planes
#m_start = start;
#m_stop = stop;
#m_start(dir) = stop(dir);
#
#port.measplanepos = m_start(dir);
#port.U_filename = [PortNamePrefix 'port_ut' int2str(portnr)];
#CSX = AddProbe(CSX, port.U_filename, 10, 'ModeFunction', E_WG_func);
#CSX = AddBox(CSX, port.U_filename, 0 ,m_start, m_stop);
#
#port.I_filename = [PortNamePrefix 'port_it' int2str(portnr)];
#CSX = AddProbe(CSX, port.I_filename, 11, 'ModeFunction', H_WG_func, 'weight', dir_sign);
#CSX = AddBox(CSX, port.I_filename, 0 ,m_start, m_stop);
#
#
#
#
