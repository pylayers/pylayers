from openems.openems import *
# A simple simulation

#
# FDTD Simulation Setting
#

F = FDTD()
F.add(Exc(typ='Sinus',f0=100000))
F.add(BoundaryCond(['PMC','PMC','PEC','PEC','MUR','MUR']))

#
# CSX (Geometry setting)
#

C = CSX()

# The Box is added as a property 

C.add(Excitation('excitation'),p=Box(P1=[-10,-10,0],P2=[10,10,0],Pr=0))
C.add(DumpBox('Et'),p=Box(P1=[-10,0,-10],P2=[10,0,30],Pr=0))
C.add(RectilinearGrid(np.arange(-10,11,1),np.arange(-10,11,1),np.arange(-10,11,1)))
C.add(Polyhedron())
S = OpenEMS(F,C)
S.save(filename='RectWaveguide.xml')
#gnd = Matter('gnd')
#sphere = Matter('sphere')
#patch = Matter('patch')
#substrate = Matter('substrate',typ='Ma',Epsilon="3.38",Kappa="0.00046")
#cdgsht = Matter('copper',typ='Cs',conductivity="56e6",thickness="40e-6")
#b1 = Box(P1=[0,0,0],P2=[100,100,200],Pr=0)
#b2 = Box(P1=[0,0,0],P2=[10,20,30],Pr=10)
#b4 = Box(P1=[-10,0,-10],P2=[10,0,30],Pr=0)
#s1 = Sphere(P=[0,0,0],R=100,Pr=50)
#dump = DumpBox()
#C.add(gnd)
#C.add(patch)
#C.add(substrate)
#C.add(sphere)
#C.add(cdgsht)
#C.add(exc)
#C.add(dump)
#C.set('gnd',b1)
#C.set('gnd',b2)
#C.set('sphere',s1)
#C.set('copper',b1)
#C.set('copper',b2)
#C.set('Et',b4)
#C.save(filename='structure.xml')
##C.AddBox(prop='ConductingSheet',name='copper',P1=[0,-50,200],P2=[1000,50,200],Pri=10)
##C.AddCylinder(prop='Metal',name='cyl0',P1=[0,0,0],P2=[0,0,100],Rad=50,Pri=10)
#
