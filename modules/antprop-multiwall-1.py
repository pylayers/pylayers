import matplotlib.pyplot as plt 
from pylayers.simul.simulem import * 
from pylayers.measures.mesuwb import *
from pylayers.antprop.multiwall import *
S = Simul()
S.layout('Lstruc.str','matDB.ini','slabDB.ini')
fGHz = 4 
Tx,Rx = ptw1()
Lwo,Lwp,Edo,Edp = Loss0_v2(S.L,Tx,fGHz,Rx[1,0:2])
fig,ax = S.L.showGs()
tit = plt.title('test Loss0_v2')
sc2 = ax.scatter(Rx[1,0],Rx[1,1],s=20,marker='x',c='k')
sc1 = ax.scatter(Tx[:,0],Tx[:,1],s=Edo,c=Edo,linewidth=0)
plt.show()
