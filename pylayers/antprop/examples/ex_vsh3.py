from pylayers.antprop.antenna import *
from numpy import *
from matplotlib.pylab import *

kf = 30
A = Antenna('S2R3.vsh3','ant')
phi = linspace(0,2*pi,180)
theta=array([1.57])
Fth,Fph = A.pattern(theta,phi)
polar(phi,abs(Fth[kf,0,:]),phi,abs(Fph[kf,0,:]))
B = Antenna('S2R3.mat','ant/UWBAN/Matfile')
polar(B.phi,abs(B.Ftheta[kf,45,:]),B.phi,abs(B.Fphi[kf,45,:]))
legend((u'$F_{\\theta}^{vsh}$',u'$F_{\phi}^{vsh}$',u'$F_{\\theta}^{original}$',u'$F_{\phi}^{original}$'),loc= 'best')
t=title('$\\theta=\\frac{\pi}{2}$'+', f = '+str(A.fa[kf])[0:6]+ ' GHz')
t.set_fontsize(18)
savefig('polarvsh3.png')
show()
