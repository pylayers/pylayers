from pylayers.em.openems.openems import *
# A simple simulation

#
# FDTD Simulation Setting
#

F = FDTD(NumberOfTimesteps="500",OverSampling="1",endCriteria="0",f_max="10000000")
F.add(Exc(typ='Sinus',f0=100000))

# Boundary Condition
F.add(BoundaryCond(['PMC','PMC','PEC','PEC','MUR','MUR']))

#
# CSX (Geometry setting)
#

C = CSX()

# Excitation

C.add(Excitation('excitation'),p=Box(P1=[-10,-10,0],P2=[10,10,0],Pr=0))

# DumpBox

C.add(DumpBox('Et'),p=Box(P1=[-10,0,-10],P2=[10,0,30],Pr=0))

# Calculus Grid

C.add(RectilinearGrid(np.arange(-10,11,1),np.arange(-10,11,1),np.arange(-10,11,1)))

# Create OpenEMS structure

S = OpenEMS(F,C)

# Save xml file

S.save(filename='ParallelPlate.xml')
