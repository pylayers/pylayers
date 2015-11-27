from pylayers.gis.layout import *
from pylayers.mobility.agent import *
from pylayers.network.network import *
from pylayers.network.emsolver import *
import matplotlib.pyplot as plt

# ## Layout Creation
#
# First we need to load a layout structure.
#
# A place where nodes from a network can move and communicate

# <codecell>

L=Layout('TA-Office.ini')
try:
    L.dumpr() # build 's'tructure, 't'opological and 'r'oom graphs
except:
    L.build()
    L.dumpw()

fig,ax=L.showG('s')


# <markdowncell>

# ## Network Creation
#
# Now we must instanciate a void Network.
#
# A Network will gather all radio acces technology RAT from all agents
# The Network class inherits from networkx

# <codecell>

N=Network()

# <markdowncell>

# ## Agents Creation
#
# Then we need to instanciate agents ( nodes ) into that layout.
# An agent is either a moving person or a static acces point which have one or more (RAT).
#
# - For that tutorial we create 3 agents in room **0 ,3 and 5** respectively.
# - We suppose they communicate on the same and unique RAT : **'rat1'**
# - Obviously agent will evolve into the previously created network **N**.
#

# <codecell>

Ag=[]
nb_agent = 3
room_init=[0,3,5]
pdb.set_trace()
for na in range(nb_agent):
    Ag.append(
              Agent(ID = na,
                    Layout = L,
                    net = N,
                    roomId = room_init[na],
                    RAT=['rat1']
                    )
              )

# <markdowncell>

# ## Connection Creation
#
# Now we can create the Network. It means, we connect agents all together, according to their common RAT.

# <codecell>

N.create()


# show the Layout
fig2,ax2=L.showG('')
# show the network layer
N.show(fig = fig2, legend = True)

# <markdowncell>

# ## Nodes and edges

# <markdowncell>

# Because Network inherits of networkx, you can have information about nodes and edges as such

# <markdowncell>

# The node is a dictionnary whi containts the folowing keys :
#
# * 'PN'  : Its Personnal Network ( described in the following)
# * 'RAT' : A list of RAT of which it belongs
# * 'p'   : Its true position
# * 'pe'  : Its estimated position if it has computed it by itself ( cf. location tutorial - IN CONSTRUCTION -)
# * 't'   : A time stamp
# * 'type': Its type ( 'ag' : for agent or 'ap' for access point )
#
# example with node '0'

# <codecell>

N.node[0]

# <markdowncell>

# The edge is a dictionnary of dictionnary.
# Each RAT is a key of the first dictionnary.
# Location dependent parameter (LDP) is the key of the second.
#
# example with edge 0-1

# <codecell>

N.edge[0][1]

# <markdowncell>

# Obviously, LDPs values are void, because no location depend parameter have been still computed.

# <markdowncell>

# # Compute Location dependent parameters (LDPs)
#
# LDPs are radio measurements bewteen agents. It could be either:
#
# 1. Time of Arrival (**TOA**)
# 2. Received Power (**Pr**)

# <markdowncell>

# ## EMS initialization
#
# We first initialize the electromagnetic solver (EMS) for the given Network N.
# We must give it a Layout strucuture, in order to be able to compute accurate LDP.
#
# .. note:: it could have be made during the Network instantiation

# <codecell>

N.EMS=EMSolver(L)

# <markdowncell>

# ## Computation
#
# The we compute *TOA* and *received power*  dependent parameters, and display for link 0-1

# <codecell>

for ldp in ['Pr','TOA']:
    N.compute_LDPs(N.nodes(),RAT='rat1',LDP=ldp)
N.edge[0][1]

# <markdowncell>

# # tests on PN

# <codecell>

N.update_PN()

N.node[0]['PN'].node[0]['pe']=np.array((4,4))
N.node[0]['PN'].node[1]['pe']=np.array((8,8))
N.node[0]['PN'].node[2]['pe']=np.array((30,8))
fig3,ax3=L.showG('')
# show the network layer
N.node[0]['PN'].show(fig = fig3, legend = True)

# <codecell>
plt.show()

