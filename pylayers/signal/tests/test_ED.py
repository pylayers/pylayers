import numpy as np
import matplotlib.pyplot as plt
from pylayers.signal.ED import *
from pylayers.signal.DF import *
import pdb
#
# Definition of Input signal
#
fsGHz = 50
fcGHz = 4
BGHz = 3.5
Tns = 20
#
# creates the energy detector
#
ed= ED(fsGHz=fsGHz, fcGHz=fcGHz, BGHz=BGHz, Tns=Tns)

w = bs.Noise(fsGHz=fsGHz, tf=1000)

pdb.set_trace()
y  = ed.apply(w)

# band limitation
ed.filterB.freqz()
# time limitation
ed.filterT.freqz()
#Hz2 = Var_w*np.real(filter1.h*np.conj(filter1.h))


w.plot(typ='v')

# plot power spectal density
w.ppsd(mask=False)

Rn = ED.filterB.wnxcorr(var=w.var)
Rx = Rn*Rn*2



f,a = Rn.plot(typ=['v'],xlabels=['sample'],figsize=(8,3))


f,a = Rx.plot(typ=['v'],xlabels=['sample'],figsize=(8,3))



Phix = np.real(fft.fft(Rx.y))
HH  = ed.filterT.H.symH(0)
Phiy = Phix * np.real(HH.y * np.conj(HH.y))


plot(fft.fftshift(Phix))
title(u'$\Phi_x(f)$')

plot(fft.fftshift(Phiy))
title(u'$\Phi_y(f)$')


print( 'empirical :', mean(y.y))
print( ' mean xQ : ', mean(ED.xq))
print( 'theoretical :',Rn.max()[0]*ED.beta*ED.bet)



print( 'theoretical : ',mean(Phiy))
print( 'empirical : ', np.var(y.y))
print( 'std : ', np.sqrt(np.var(y.y)))



order_e = (2*mean(y.y)**2)/np.var(y.y)
order_t = (2*(Rn.max()[0]*ED.beta*ED.beta)**2)/mean(Phiy)
scale_e = np.sqrt(np.var(y.y)/(2*order_e))
scale_t = np.sqrt(mean(Phiy)/(2*order_t))
print( "ordre empirical : ",order_e)
print( "ordre theoretical : ",order_t)
print( scale_e)
print (scale_t)


lchi2_e = st.chi2(df=order_e,loc=0,scale=scale_e)
lchi2_t = st.chi2(df=order_t,loc=0,scale=scale_t)



h = hist(y.y,200,normed=True)
t = h[1]
thl_e = lchi2_e.pdf(t)
thl_t = lchi2_t.pdf(t)
plot(t,thl_e,'r',label='empirical')
plot(t,thl_t,'g',label='theoretical')
plt.legend()



order = (2*BGHz*Tns+1)/2
print( order)


ip=bs.EnImpulse(fe=fsGHz)
ip.translate(10)
ip.help(typ='mb')
print( ip.energy())




IP = ip.psd(Tpns=500)
IP.plotdB(mask=True)




IP = ip.psd




flt=DF()
fc = 4 
fs = fsGHz
fN = 50 
B = 0.5
w1 = (fc-B/2)/fN
w2 = (fc+B/2)/fN
flt.butter(w=[w1,w2],typ='bandpass')
flt.freqz()


### Creating a simple channel



#fig = plt.figure(figsize=(10,3))
#f,a=ip.plot(typ=['v'],fig=fig)
tauk = np.array([4,7,8,10,12])
ak = np.array([0.2,0.5,0.8,0.2,0.1735])
sum(ak*ak)
stem(tauk,ak)
try:
    del(s)
except:
    pass
for k in range(len(ak)):
    ip =  bs.EnImpulse(fe=fsGHz)
    ip.translate(tauk[k])
    try:
        s = s + ip*ak[k]
    except:
        v = ak[k]
        s = ip*v
s2 = bs.TUsignal(s.x,s.y)
Etx= ys.y[-1]
EtxdB = 10*np.log10(Etx)
print( EtxdB)
dm = 4
PL = 32.4+20*np.log10(dm)+20*np.log10(fcGHz)
ErxdB = EtxdB-PL
print( ErxdB)
Erx = 10**(ErxdB/10)
print( Erx)
s = s*(10**(-PL/10))




s.plot(typ='v')




s.energy()




ys = ED.apply(s)




s3=bs.TUsignal(s.x,ED.xq)




s3.plot(typ='v')




ys.plot(typ=['v'])














s.y.shape




s.plot(typ='v')
w.plot(typ='v',c='r')




sw = s + w
sw.plot(typ='v')




ysw = ED.apply(sw)




y.plot(typ='v')




ysw.plot(typ='v')




#h1 = hist(y.y,200,normed=True)
#h2 = hist(ysw.y,200,normed=True)
#t = h2[1]
printi( Erx)
lchi2_e = st.chi2(df=order_e,loc=0,scale=scale_e)
lchi2_t = st.chi2(df=order_t,loc=0,scale=scale_t)
chi2nc_e = st.ncx2(df=order_e,nc=100.4,scale=scale_e,)
chi2nc_t = st.ncx2(df=order_t,nc=100.4,scale=scale_t)
t = linspace(0,5e-9,150)
thl_e = lchi2_e.pdf(t)
thl_t = lchi2_t.pdf(t)
thlc_e = chi2nc_e.pdf(t)
thlc_t = chi2nc_t.pdf(t)
plot(thl_e,'r',label='empirical')
plot(thl_t,'g',label='theoretical')
plot(thlc_e,'k',label='empirical')
plot(thlc_t,'m',label='theoretical')
plt.legend()




chi2nc_t = st.ncx2




chi2nc_t = st.ncx2






