import vtk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams['animation.ffmpeg_path']='/usr/bin/ffmpeg'

from vtk.util.numpy_support import vtk_to_numpy
reader = vtk.vtkXMLRectilinearGridReader()

fig,ax = plt.subplots()

NrTs = 1500
ims = []
for k in np.arange(0,NrTs,1):
    # create the filename from sequence number
    if k<10:
        filename="tmp/Et_000000000"+str(k)+'.vtr'
    elif k<100:
        filename="tmp/Et_00000000"+str(k)+'.vtr'
    elif k<1000:
        filename="tmp/Et_0000000"+str(k)+'.vtr'

    reader.SetFileName(filename)
    reader.Update()
    Efield = reader.GetOutput().GetPointData().GetArray(0)
    u=vtk_to_numpy(Efield)
    im = ax.imshow(u,interpolation='nearest')
    a = im.get_axes()
    # Draw colorbar only once
    if k==0:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right",size="5%",pad=0.05)
        clb = fig.colorbar(im,cax)
        clb.set_label("Efield V/m")
    ax.axis('tight')
    ims.append([im])

ani = animation.ArtistAnimation(fig,ims,interval=50,blit=True,repeat_delay=1000)
#ani.save('Efield.mp4')
