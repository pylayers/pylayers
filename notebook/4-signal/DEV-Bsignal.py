# In[1]:
from pylayers.signal.bsignal import *
from pylayers.simul.simulem import *
figsize(8,8)

# # Generation of an Impulse of normalized energy 
# 
# One possible manner to define an ernergy normalized short UWB impulse is as follows. 
# 
# \\[ p(t)= \frac{\sqrt{2\sqrt{2}}}{\tau\sqrt{\pi}} \cos(2\pi f_c t) e^{-(\frac{t}{\tau})^2} \\]
# with
# 
# \\[ \tau = \frac{2}{B\pi}\sqrt{\frac{\gamma_{dB}\ln{10}}{20}}\\]
# 
# where \\(B\\) is the desired bandwidth defined at \\(\gamma_{dB}\\) below the spectrum maximum and \\(f_c\\) is the central frequency of the pulse.
# 
# 
# This waveform is a gaussian windowing
# of a sine wave of frequency $f_c$. The normalization term depends on the exponential scaling factor $\tau$. 
# 

# In[2]:
fc     = 4 
band   = 2
thresh = 10
fe     = 100 
ip =EnImpulse([],fc,band,thresh,fe)

# In[3]:
ip.info()

# Out[3]:
#     TUsignal
#     --------
#     shx :  (343,)
#     shy :  (343,)
#     dx :   0.01
#     xmin : -1.71
#     xmax : 1.71
#     ymin : -1.89545539648
#     ymax : 2.16154131873
# 

# ## Verification of energy normalization in both domains

# In[4]:
E1= sum(ip.y*ip.y)*ip.dx()
print "Integration in time",E1

# Out[4]:
#     Integration in time 1.0
# 

# In[5]:
P = ip.esd()
E2 = sum(P.y)*P.dx()
print "Integration in frequency domain ",E2

# Out[5]:
#     Integration in frequency domain  1.0
# 

# ## Calcul of UWB channel impulse response

# In[6]:
S= Simul()
S.load('where2.ini')

# In[7]:
st = S.wav.st
sf = S.wav.sf
S.wav.info()

# Out[7]:
#     tw  :  30.0
#     band  :  4.0
#     fc  :  4.493
#     thresh  :  3.0
#     fe  :  50.0
#     Np  :  1500.0
#     te  :  0.02
#     type  :  file
# 

# Here the time domain waveform is measured and the anticausal part of the signal is artificially set to 0.
# 
# To handle properly the time domain wavefom it is required to center the signal in the middle of the array. 
# 
# `st` stands for signal in time domain

### Ploting the waveform

### In time domain

# In[33]:
S.wav.st.plot()


# Out[33]:
# image file:

### in frequency domain 

# The frequency domain version of the signal is embedded in the same object. 
# 
# `sf` stands for signal in frequency domain. 

# In[8]:
f,ax=S.wav.sf.plot()

# Out[8]:
# image file:

### Construction of the propagation channel 

# The link between Txid = 1 and Rxid =1 is simply loaded as

# In[9]:
vc = S.VC(1,1)

# Out[9]:
#     nray :  500
#     nfreq :  181
#     nb rayons dans .tauk :  500
#     nb rayons 2:  500
# 

# The following representation shows the spatial spreading of the propagation channel. 
# On the left are scattered the intensity of rays wrt to angles of departure (in azimut and elevation). 
# On the right is the intensity of rays wrt to angles of arrival. It misses the application between the 2
# planes as well as the delay dimension of the propagation channel.

# In[37]:
vc.doadod()

# Out[37]:
# image file:

# ##  Construction of the transmission channel

# The transmission channel is obtain from the combianation of the propagation channel and the vector antenna pattern at bot side of the radio link

# In[10]:
sc = vc.vec2scal()

# The ScalChannel object contains all the information about the ray transfer functions. 
# The transmission channel is obtained by applying a vector radiation pattern using an antenna file.
# In the presented case, it comes from a real antenna which has been used during the FP7 WHERE1 measurement campaign 
# M1.

# In[11]:
S.tx.A.info()

# Out[11]:
#     defant.vsh3
#     type :  vsh3
#     --------------------------
#     fmin (GHz) : 2.0
#     fmax (GHz) : 8.0
#     Nf   : 121
#     Br
#     -------------
#     Nf   :  121
#     fmin (GHz) :  2.0
#     fmax (GHz) :  8.0
#     Ncoeff s3 :  18
#     Bi
#     -------------
#     Nf   :  121
#     fmin (GHz) :  2.0
#     fmax (GHz) :  8.0
#     Ncoeff s3 :  18
#     Cr
#     -------------
#     Nf   :  121
#     fmin (GHz) :  2.0
#     fmax (GHz) :  8.0
#     Ncoeff s3 :  18
#     Ci
#     -------------
#     Nf   :  121
#     fmin (GHz) :  2.0
#     fmax (GHz) :  8.0
#     Ncoeff s3 :  18
# 

# In[13]:
f,ax=sc.H.plot()


# Out[13]:
# image file:

# The antenna can also been taken into account

# In[14]:
alpha = 1./sqrt(30)  # scaling constant depends on how are stored the antenna data
sca = vc.vec2scalA(S.tx.A,S.rx.A,alpha)
sca.H.plot()

# Out[14]:
#     (<matplotlib.figure.Figure at 0xbc650cc>,
#      array([Axes(0.125,0.547727;0.775x0.352273),
#            Axes(0.125,0.125;0.775x0.352273)], dtype=object))

# image file:

# ## Calculate UWB Channel Impulse Response 

# In[15]:
cir = sc.applywavB(S.wav.sfg)

# In[16]:
cir.plot()

# Out[16]:
# image file:

# In[17]:
#CIR=cir.esd(mode='unilateral')
#CIR.plot()

# Out[17]:
### Hermitian symetry enforcment 

# If the number of point for the transmission channel and the waveform were the same the mathematical operation is an Hadamrd-Shur product between 
# $\mathbf{Y}$ and $\mathbf{W}$. 
# 
# $\mathbf{Y} = \mathbf{S} \odot \mathbf{W}$
# 
# In practice this is what is done after a resampling of the time base with the greater time step. 

# The process whic consist in going time domain to frequency domain is delegate to a specialized class which maintain the proper 
# binding between signal samples and their indexation either in time or in frequency domain.

# In[21]:
wgam = S.wav.sfg
Y    = sc.apply(wgam)
tau  = Y.tau0

# The transmission channel has a member data which is the time delay of each path

# In[25]:
print 'tau =', tau

# Out[25]:
#     tau = [ 23.86713221  25.90353456  26.71456706  31.5656444   36.76975426
#       33.13220501  29.70841516  33.77007551  38.12306298  38.67871406
#       27.8861973   27.64609677  40.8039797   36.1848012   44.63339207
#       40.97791791  26.00589677  56.03851137  56.50109556  28.40727606
#       45.75470239  46.21869515  28.64099626  31.45147939  56.93564665
#       57.39100019  57.30918643  57.76159531  49.6141571   31.23879355
#       45.89051822  25.74826855  34.70312186  34.5104815   48.01134886
#       59.20577447  58.76448492  39.86471567  35.1232117   45.28595588
#       40.03159863  35.31250959  46.67807592  55.69281322  33.21104234
#       65.58905058  65.1909849   33.00969555  41.92951758  63.97682332
#       64.44036847  43.69434237  42.08821424  43.38811921  42.01423039
#       37.6279447   53.22774303  37.45035232  41.69566991  69.27701651
#       69.74102761  46.85683446  54.41707085  49.77517159  64.76408746
#       65.22203762  51.0450064   65.09271861  65.54837284  70.4639201
#       70.0047007   70.76608669  70.30884091  42.25683377  42.09877275
#       48.44775243  38.19069113  39.32593151  50.54234748  66.37758547
#       66.8244792   50.67407841  71.94971066  71.5000346   68.36337908
#       67.98156041  46.89321202  45.5860115   59.98403491  47.21435752
#       38.01572774  39.49509112  69.29340929  58.934858    59.54859827
#       45.29257972  47.03516415  60.09507282  60.01118612  59.6604466
#       47.35534699  48.9702404   58.11322282  58.55933363  72.12893787
#       72.54040836  49.10618879  77.28751962  76.86907554  58.2278284
#       58.67306783  51.16678171  66.29124293  58.31307211  43.20334349
#       50.48955888  65.89749616  66.39173309  77.23599638  77.70053523
#       56.52016356  55.17961088  56.97789045  42.89361517  65.99858584
#       45.809703    56.63799275  88.72342731  89.65357417  89.18849218
#       77.88935896  78.35002411  56.86076562  78.16282511  49.04279311
#       78.62188797  55.26377154  46.30185981  54.0932939   45.66394152
#       46.44561934  46.76033457  57.99665686  55.38427374  50.24950477
#       54.21639768  47.80371906  89.29277494  90.21705196  89.75488978
#       89.5314171   90.45325563  72.3931858   89.99230651  79.65543756
#       47.66405587  79.24949727  72.0327308   50.70038352  75.90628341
#       80.45504151  75.51315577  57.43282066  64.89406137  50.83170489
#       67.51358711  70.66550806  67.0713716   58.19389048  51.94736653
#       63.56523608  57.74495649  79.2360339   79.68891458  46.90268886
#       69.03683397  67.61226056  68.6588345   67.17069467  52.07554342
#       68.1027834   69.13333333  68.7558644   69.95773327  91.38223729
#       90.46986507  57.31662547  90.92599776  70.05296409  66.25100419
#       52.13525146  84.11271541  84.53947558  47.98186891  59.71670155
#       82.66784744  77.88221877  82.35237675  66.35155529  77.46706827
#       60.77898028  73.17464649  95.64158676  94.77022903  95.20576082
#       57.86029151  58.3083375   59.82823562  72.76684227  77.96777112
#       77.5530786   65.80030395  60.88856853  72.14362219  73.26569608
#       71.41643135  70.96346478  75.05823665  72.58922479  62.36672794
#       72.85840148  61.27669898  64.79979424  70.03675226  72.23597134
#       90.45565753  90.44995498  89.97290449  65.90154272  83.13970461]
# 

# symHz force the Hermitian symetry of Y with zero padding of 500 

# In[26]:
UH   = Y.symHz(500)
uh   = UH.ifft(1)
UH.plot()
plt.figure()
uh.plot()

# Out[26]:
# image file:

# image file:

# In[100]:
ips  = Y.ift(500,1)
t    = ips.x 
ip0  = TUsignal(t,ips.y[0,:])

# In[101]:
plot(UH.x,real(UH.y[0,:]),UH.x,imag(UH.y[0,:]))
U0 = FHsignal(UH.x,UH.y[0,:])
u0 = U0.ifft(1)
u1 = ifft(U0.y)
plt.figure()
plot(uh.x,uh.y[0,:]*1000+3)
S.wav.st.plot()

# Out[101]:
#     False False

# image file:

# image file:

# In[102]:
U0.plot()

# Out[102]:
# image file:

# # Here is the problem 
# 
# For some reason the Hermitian symmetry forcing is not working here

# In[103]:
U1=u0.fft()
g = fft(u1)
plot(abs(g))
plt.figure()
s  = fftshift(u1)
plot(abs(g))
#plt.figure()
#plot(uh.x,uh.y[0,:])
#plot(uh.x,s*50+0.003)

# Out[103]:
#     [&lt;matplotlib.lines.Line2D at 0xb3b3234c&gt;]

# image file:

# image file:

# In[104]:
plot(abs(fft(s)),'r')
plot(abs(fft(uh.y[0,:])),'g')

# Out[104]:
#     [&lt;matplotlib.lines.Line2D at 0xb2d78e2c&gt;]

# image file:

# In[105]:
wgam.plot()

# Out[105]:
# image file:

# In[106]:
S.wav.sf.plot()

# Out[106]:
# image file:

# In[107]:
print uh.y[0,:]

# Out[107]:
#     [  6.97818719e-19   2.43627115e-10   6.40438282e-10 ...,   1.08726676e-09
#        5.07489405e-10   8.53183024e-11]

# In[108]:
plot(imag(s))

# Out[108]:
#     [&lt;matplotlib.lines.Line2D at 0xb4f71aac&gt;]

# image file:

# Problem $s$ is not real 

# In[109]:
u0

# Out[109]:
#     &lt;Bsignal.TUsignal at 0xb3b36eec&gt;

# In[110]:
plot(real(u0.y))

# Out[110]:
#     [&lt;matplotlib.lines.Line2D at 0xb03296cc&gt;]

# image file:

# In[111]:
plot(imag(s))

# Out[111]:
#     [&lt;matplotlib.lines.Line2D at 0xad6ab24c&gt;]

# image file:

# In[112]:
U0.y

# Out[112]:
#     array([ 0.+0.j,  0.+0.j,  0.+0.j, ...,  0.+0.j,  0.+0.j,  0.+0.j])

# In[126]:
plot(real(U0.y))

# Out[126]:
#     [&lt;matplotlib.lines.Line2D at 0xb51b7dcc&gt;]

# image file:

# In[114]:
U0.y[0]

# Out[114]:
#     0j

# In[115]:
U0.y[50]

# Out[115]:
#     0j

# In[116]:
U0.y[-50]

# Out[116]:
#     0j

# In[117]:
UH.y[0,2]

# Out[117]:
#     0j

# In[118]:
UH.y[0,-2]

# Out[118]:
#     0j

# In[119]:
N = len(UH.y)


# In[120]:
v1 = UH.y[1:(N-1)/2.]
v2 = UH.y[N:-1:(N-1)/2.]

# In[121]:

len(v1)

# Out[121]:
#     150

# In[122]:
len(v2)

# Out[122]:
#     0

# In[123]:
UH.y[0,-1]

# Out[123]:
#     0j

# In[124]:
UH.y[0,1]

# Out[124]:
#     0j

# In[125]:
plot(real(UH.y[0,:]))
plot(imag(UH.y[1,:]))

# Out[125]:
#     [&lt;matplotlib.lines.Line2D at 0xad6bb10c&gt;]

# image file:

# In[49]:
cir.x


# Out[49]:
#     array([  -9.99307479,   -9.97922438,   -9.96537396, ...,  105.58864266,
#             105.60249307,  105.61634349])

# In[50]:
cir.y

# Out[50]:
#     array([  0.00000000e+00,   0.00000000e+00,   0.00000000e+00, ...,
#             -1.98587668e-12,  -1.43517473e-12,  -7.48555146e-13])

# In[51]:
cir.plot()

# Out[51]:
# image file:

# In[ ]:

