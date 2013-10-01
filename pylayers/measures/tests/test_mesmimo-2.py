from pylayers.measures.mesmimo import *
import matplotlib.pyplot as plt 
import numpy as np 
import numpy.fft as fft
import matplotlib.animation as animation
from pylayers.gis.layout import *

L=Layout('11Dbibli.ini')
#
#V = MIMO('S00_11_V.csv','12-9-2013/S0/',calibration=True);
print "lecture" 
M = MIMO('S00_11.csv','12-9-2013/S0/',calibration=True,Nz=20000);
#D = M-V
#M.mulcplot('t',types=['v'],color='k')
print "calbento" 
M.calBento2(time=True)
#D.calBento2(time=True)
x = np.arange(1,7,0.1)
y = np.arange(1,7,0.1)
Me = np.meshgrid(x,y)
Me = np.dstack((Me[0],Me[1]))
###D.grid(Me)
###L.showG('s')
###fig,ax = M.C.plot(types=['r'],titles=titles,xlabels=xlabels,ylabels=ylabels)
###fig2,ax2 = M.Hcal.plot(types=['m'],titles=titles,xlabels=xlabels,ylabels=ylabels,fig=fig1,ax=ax1,color='r',labels=lB)
###print "calcul h"
#D.calBento2(duR=0.05,duT=0.05,time=True)
#V.calBento2(duR=0.05,duT=0.05,time=True)
#M.calBento2(duR=0.05,duT=0.05,time=True)
##D.calBento2(duR=0.2,duT=0.2,time=True)
##.calBento2(duR=0.2,duT=0.2,time=True)
##M.calBento2(duR=0.2,duT=0.2,time=True)
print "grid"
M.grid(Me)
M.showgrid()
#print "animgrid"
#D.animgrid1(layout=L)
#V.grid(Me)
#M.grid(Me)
###plt.ion()
###print "animation"
###ims = []
###fig = plt.figure()
###ax1 = fig.add_subplot(131)
###ax2 = fig.add_subplot(132)
###ax3 = fig.add_subplot(133)
###im1 = ax1.scatter(D.grid[...,0],D.grid[...,1],c=abs(V.gloc[:,769]),s=100,vmin=0,vmax=0.01)
###im2 = ax2.scatter(D.grid[...,0],D.grid[...,1],c=abs(M.gloc[:,769]),s=100,vmin=0,vmax=0.01)
###im3 = ax3.scatter(D.grid[...,0],D.grid[...,1],c=abs(D.gloc[:,769]),s=100,vmin=0,vmax=0.001)
###im1 = ax1.imshow(abs([:,:,770]),vmin=0,vmax=0.1)
###im2 = ax2.imshow(abs(Mtt[:,:,770]),vmin=0,vmax=0.1)
###im3 = ax3.imshow(abs(Dtt[:,:,770]),vmin=0,vmax=0.001)
###while 1:
###    for k in range(770,850):
###        t = ax1.title
###        t.set_text(str(k))
###        im1.set_array(abs(V.gloc[:,k]))
###        im2.set_array(abs(M.gloc[:,k]))
###        im3.set_array(abs(D.gloc[:,k]))
##        #im1.set_data(abs(V.gloc[:,k]))
##        #im2.set_data(abs(M.gloc[:,k]))
##        #im3.set_data(abs(D.gloc[:,k]))
###        plt.draw()
###
###
###for k in range(750,900):
###    im = ax1.scatter(D.grid[...,0],D.grid[...,1],c=abs(D.gloc[:,769]),s=100,vmin=0,vmax=0.001)
###    #im = plt.imshow(C[:,:,k],vmin=0,vmax=0.1
###    ax = im.get_axes()
###    #t  = ax.title
###    #t.set_text(str(k))
####    ims.append([im])
###    ims.append([im])
###
####fig = plt.gcf()
###anim = animation.ArtistAnimation(fig, ims, interval=60, blit=True, repeat_delay=1000)
##f = np.linspace(1.8,2.2,1601)
##df = f[1]-f[0]
##t = np.linspace(-1./(2*df),1./(2*df),1601)
#fig = plt.figure(figsize=(20,20))
#ax1 = fig.add_subplot(111)
#fig,ax1=L.showG('s',nodes=False,fig=fig,ax=ax1)
#Nframe = D.gloc.y.shape[1]
#scat = ax1.scatter(D.grid[...,0],D.grid[...,1],c=abs(D.gloc.y[:,0]),s=100,vmin=0,vmax=0.032,linewidth=0)
#cb   = plt.colorbar(scat)
#delay_template = '%d : tau = %5.2f (ns) d= %5.2f (m)'
#delay_text  = ax1.text(0.25,0.9,'',transform=ax1.transAxes,fontsize=18)
### initialization function: plot the background of each frame
#def init():
#    delay_text.set_text('')
#    scat.set_array(abs(D.gloc.y[:,0]))
#    return scat,delay_text
##
### animation function.  This is called sequentially
#def animate(i):
#    delay_text.set_text(delay_template%(i,D.gloc.x[i],D.gloc.x[i]*0.3))
#    scat.set_array(abs(D.gloc.y[:,i]))
#    return scat,delay_text
#    #return delay_text,
##
### call the animator.  blit=True means only re-draw the
### parts that have changed.
#anim = animation.FuncAnimation(fig, animate, init_func=init, frames=Nframe, interval=1, blit=True)
##
### save the animation as an mp4.  This requires ffmpeg
### or mencoder to be
### installed.  The extra_args ensure that the x264
### codec is used, so that
### the video can be embedded in html5.  You may need to
### adjust this for
### your system: for more information, see
### http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('S11.mp4', fps=30)
###, extra_args=['-vcodec', 'libx264'])
#plt.show()
