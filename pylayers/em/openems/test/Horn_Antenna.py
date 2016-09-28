"""
Tutorials / horn antenna

Description at:
http://openems.de/index.php/Tutorial:_Horn_Antenna

 (C) 2011,2012,2013 Thorsten Liebig <thorsten.liebig@uni-due.de>
 Python Adaptation : ESIR Project 2015

"""
from pylayers.em.openems.openems import *
import scipy.constants as cst
import numpy as np
# setup the simulation

unit = 1e-3 # all length in mm
class HornAntenna(object):
    def __init__(self,**kwargs):
        defaults = {'unit'   : 1e-3,
                    'width'  : 20,
                    'height' : 30 ,
                    'length' : 50 ,
                    'feed_length' : 50 ,
                    'thickness' : 2,
                    'angle' : np.array([20,20])*np.pi/180.
        }
        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self.unit = kwargs['unit']
        self.width = kwargs['width']
        self.height = kwargs['height']
        self.length = kwargs['length']
        self.feed_length = kwargs['feed_length']
        self.thickness = kwargs['thickness']
        self.angle = kwargs['angle']


HA = HornAntenna()
# size of the simulation box
SimBox = np.r_[200,200,200]

# frequency range of interest

f_start =  10e9
f_stop  =  20e9

# frequency of interest

f0 = 15e9

#waveguide TE-mode definition

TE_mode = 'TE10'

a = HA.width
b = HA.height

# setup FDTD parameter & excitation function

F = FDTD(EndCriteria="1e-4")
F.add(Exc(typ='Gaussian',f0=0.5*(f_start+f_stop),fc=0.5*(f_stop-f_start)))


F.add(BoundaryCond(['PML 8','PML 8','PML 8','PML 8','PML 8','PML 8']))


# setup CSXCAD geometry & mesh
# currently, openEMS cannot automatically generate a mesh
max_res = ((cst.c/f_stop)/unit)/15. # cell size: lambda/20
C = CSX()
#
# Warning : It is not the same thing to add a new properties (add) and to add
# a new primitive to an existing property (primitive)
#
C.add(Matter('horn',
             p=Box(
                 P1=[-a/2.-HA.thickness,-b/2.,0],
                 P2=[-a/2.,-b/2.,0],Pr=10)
             ))

#
# Define Mesh
#

linex = [-SimBox[0]/2.,-a/2., a/2., SimBox[0]/2.]
meshx = SmoothMeshLine( linex, max_res, 1.4)

liney = [-SimBox[1]/2., -b/2., b/2., SimBox[1]/2.]
meshy = SmoothMeshLine( liney, max_res, 1.4 )

linez = [-HA.feed_length, 0 ,SimBox[2]-HA.feed_length ]
meshz = SmoothMeshLine( linez, max_res, 1.4 )

C.add(RectilinearGrid(meshx,meshy,meshz))

#
# Waveguide
#
C.primitive('horn',Box(
                 P1=[-a/2.-HA.thickness,-b/2.,meshz[0]],
                 P2=[-a/2.,b/2.,0],Pr=10)
             )
C.primitive('horn',Box(
                 P1=[a/2.+HA.thickness,-b/2.,meshz[0]],
                 P2=[a/2.,b/2.,0],Pr=10)
             )
C.primitive('horn', Box(
                 P1=[-a/2.-HA.thickness,b/2.+HA.thickness,meshz[0]],
                 P2=[a/2.+HA.thickness,b/2.,0],Pr=10)
             )
C.primitive('horn', Box(
                 P1=[-a/2.-HA.thickness,-b/2.-HA.thickness,meshz[0]],
                 P2=[a/2.+HA.thickness,-b/2.,0],Pr=10)
             )

#
# horn opening 4 metallic plates
#

horn_opening1 = np.array([[0, HA.length, HA.length, 0],
          [a/2.,
           a/2 + np.sin(HA.angle[0])*HA.length,
          -a/2 - np.sin(HA.angle[0])*HA.length,
          -a/2.]])

horn_opening2 = np.array([[b/2+HA.thickness,
          b/2+HA.thickness + np.sin(HA.angle[1])*HA.length,
          -b/2-HA.thickness - np.sin(HA.angle[1])*HA.length,
          -b/2-HA.thickness],
          [ 0, HA.length, HA.length, 0]])

L1 = LinPoly(lp=horn_opening1.T,Pr=10)
L2 = LinPoly(lp=horn_opening1.T,Pr=10)
L3 = LinPoly(lp=horn_opening2.T,Pr=10,normdir=0)
L4 = LinPoly(lp=horn_opening2.T,Pr=10,normdir=0)

T1 = Transformation()
T2 = Transformation()
T3 = Transformation()
T4 = Transformation()

# y translate
Tr1 = Translate([0,-b/2-HA.thickness/2,0])
Tr2 = Translate([0,b/2+HA.thickness/2,0])
# x translate
Tr3 = Translate([-a/2-HA.thickness/2,0,0])
Tr4 = Translate([a/2+HA.thickness/2,0,0])

Rx1 = Rotate_X(HA.angle[1])
Rx2 = Rotate_X(-HA.angle[1])
Rx3 = Rotate_Y(-HA.angle[1])
Rx4 = Rotate_Y(HA.angle[1])

T1.append(Rx1)
T1.append(Tr1)

T2.append(Rx2)
T2.append(Tr2)

T3.append(Rx3)
T3.append(Tr3)

T4.append(Rx4)
T4.append(Tr4)

L1.append(T1)
L2.append(T2)
L3.append(T3)
L4.append(T4)

C.primitive('horn',L1)
C.primitive('horn',L2)
C.primitive('horn',L3)
C.primitive('horn',L4)

## first ProbeBox
#C.add(ProbeBox(name='port_ut1', Type='wv', Weight='1'),
#               a=Attributes([(0*cos(0.15708*(x--10))*sin(0*(y--15))),
#                    (-0.05*sin(0.15708*(x--10))*cos(0*(y--15))),0]),
#               p=Box(P1=[-10,-15,-25],P2=[10,15,-25],Pr=0)
#
## second ProbeBox
#
#C.add(ProbeBox(name='port_it1', Type='wc', Weight='1'), a=Attributes([(0.05*sin(0.15708*(x--10))*cos(0*(y--15))),0*cos(0.15708*(x--10))*sin(0*(y--15))),0]), p=Box(P1=[-10,-15,-25],P2=[10,15,-25],Pr=0)
#
#
##
A = (a + 2*np.sin(HA.angle[0])*HA.length)*unit * (b + 2*np.sin(HA.angle[1])*HA.length)*unit;
##
## apply the excitation
start=[-a/2, -b/2 ,meshz[7] ];
stop =[ a/2,  b/2 ,meshz[0]+HA.feed_length/2. ];

C.add(Excitation('port_excite_1',typ="Es",excite="1,1,0"))
# AddRectWaveGuidePort( CSX, 0, 1, start, stop, 2, a*unit, b*unit, TE_mode, 1);
##
##%% nf2ff calc
##start = [mesh.x(9) mesh.y(9) mesh.z(9)];
##stop  = [mesh.x(end-8) mesh.y(end-8) mesh.z(end-8)];
##[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'Directions', [1 1 1 1 0 1]);
##
##%% prepare simulation folder
##Sim_Path = 'tmp_Horn_Antenna';
##Sim_CSX = 'horn_ant.xml';
##
##[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
##[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder
##
##%% write openEMS compatible xml-file
##WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);
##
##%% show the structure
##CSXGeomPlot([Sim_Path '/' Sim_CSX]);
##
##%% run openEMS
##RunOpenEMS(Sim_Path, Sim_CSX);
##
##%% postprocessing & do the plots
##freq = linspace(f_start,f_stop,201);
##
##port = calcPort(port, Sim_Path, freq);
##
##Zin = port.uf.tot ./ port.if.tot;
##s11 = port.uf.ref ./ port.uf.inc;
##
##plot( freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
##ylim([-60 0]);
##grid on
##title( 'reflection coefficient S_{11}' );
##xlabel( 'frequency f / GHz' );
##ylabel( 'reflection coefficient |S_{11}|' );
##
##drawnow
##
##%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
##% calculate the far field at phi=0 degrees and at phi=90 degrees
##thetaRange = (0:2:359) - 180;
##disp( 'calculating far field at phi=[0 90] deg...' );
##nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, [0 90]*pi/180);
##
##Dlog=10*log10(nf2ff.Dmax);
##G_a = 4*pi*A/(c0/f0)^2;
##e_a = nf2ff.Dmax/G_a;
##
##% display some antenna parameter
##disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
##disp( ['directivity: Dmax = ' num2str(Dlog) ' dBi'] );
##disp( ['aperture efficiency: e_a = ' num2str(e_a*100) '%'] );
##
##%%
##% normalized directivity
##figure
##plotFFdB(nf2ff,'xaxis','theta','param',[1 2]);
##drawnow
##%   D_log = 20*log10(nf2ff.E_norm{1}/max(max(nf2ff.E_norm{1})));
##%   D_log = D_log + 10*log10(nf2ff.Dmax);
##%   plot( nf2ff.theta, D_log(:,1) ,'k-', nf2ff.theta, D_log(:,2) ,'r-' );
##
##% polar plot
##figure
##polarFF(nf2ff,'xaxis','theta','param',[1 2],'logscale',[-40 20], 'xtics', 12);
##drawnow
##%   polar( nf2ff.theta, nf2ff.E_norm{1}(:,1) )
##
##%% calculate 3D pattern
##phiRange = sort( unique( [-180:5:-100 -100:2.5:-50 -50:1:50 50:2.5:100 100:5:180] ) );
##thetaRange = sort( unique([ 0:1:50 50:2.:100 100:5:180 ]));
##
##disp( 'calculating 3D far field...' );
##nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, phiRange*pi/180, 'Verbose',2,'Outfile','nf2ff_3D.h5');
##
##figure
##plotFF3D(nf2ff);
##
##%%
##E_far_normalized = nf2ff.E_norm{1}/max(nf2ff.E_norm{1}(:));
##DumpFF2VTK([Sim_Path '/Horn_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,'scale',1e-3);
S = OpenEMS(F,C)
#
S.save(filename='HornAntenna.xml')
