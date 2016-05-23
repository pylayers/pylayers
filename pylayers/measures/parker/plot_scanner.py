import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
from matplotlib.path import Path
import matplotlib.lines as mlines
from matplotlib.collections import PatchCollection
pA = np.array([0,0])
pB = np.array([0,1.2])
pC = np.array([0.7,0])
pD = np.array([0.7,1.2])
pO = (pA+pB+pC+pD)/4

verts = [
         (pA[0], pA[1]), # left, bottom
         (pB[0], pB[1]), # left, top
         (pD[0], pD[1]), # right, top

         (pC[0],pC[1]), # right, bottom
         (pA[0],pA[1]), # ignored
         ]

codes = [Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
        ]

path = Path(verts, codes)
patch = patches.PathPatch(path,facecolor='blue',lw=2,alpha=0.3)
ax = plt.gca()
circle = patches.Circle(pO,0.1,fc='orange',ec=None)
ax.add_patch(circle)
ax.add_patch(patch)
a1 = ax.arrow(pA[0],pA[1],(pB[0]-pA[0]),(pB[1]-pA[1]),head_width=0.05,head_length=0.1,fc='k',ec='k')
a2 = ax.arrow(pC[0],pC[1],(pD[0]-pC[0]),(pD[1]-pC[1]),head_width=0.05,head_length=0.1,fc='k',ec='k')
ax.axis('equal')
plt.show()
