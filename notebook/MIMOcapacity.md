## Uniform capacity and waterfilling capacity

```python
>>> from pylayers.measures.mesmimo import *
>>> from pylayers.gis.layout import *
>>> from matplotlib.pylab import *
>>> from pylayers.signal.bsignal import *
>>> %matplotlib inline
```

We read a MIMO matrix from the measurement dataset and invoke the capacity method.

```python
>>> #Hp = Hc.Hcal.y.swapaxes(1,2)
... #Hpp = Hp.swapaxes(0,1)
... #Hd = Hpp.swapaxes(1,2)
```

```python
>>> Hc =MIMO('S00_11_V.csv','/MIMO-FB/12-9-2013/S00/')
>>> Pt = np.logspace(-9,-3,100)
>>> rho,CB = Hc.Bcapacity(Pt=Pt)
>>> fGHz = Hc.Hcal.x
>>> dfGHz = (fGHz[1]-fGHz[0])
>>> rho=np.sum(rho*dfGHz,axis=0)
>>> rho=np.sum(rho,axis=0)
>>> HdH,U,S,V = Hc.transfer()
>>> HdH = HdH.swapaxes(0,2)
>>> rg = np.real(np.sqrt(np.trace(HdH)/(32)))
```

```python
>>> H2.x
array([ 1.8    ,  1.80025,  1.8005 , ...,  2.1995 ,  2.19975,  2.2    ])
```

```python
>>> H2 = FUsignal(x=fGHz,y=rg)
```

```python
>>> C,SNR=H2.capacity(Pt=Pt)
```

```python
>>> Pt
array([  1.00000000e-09,   1.14975700e-09,   1.32194115e-09,
         1.51991108e-09,   1.74752840e-09,   2.00923300e-09,
         2.31012970e-09,   2.65608778e-09,   3.05385551e-09,
         3.51119173e-09,   4.03701726e-09,   4.64158883e-09,
         5.33669923e-09,   6.13590727e-09,   7.05480231e-09,
         8.11130831e-09,   9.32603347e-09,   1.07226722e-08,
         1.23284674e-08,   1.41747416e-08,   1.62975083e-08,
         1.87381742e-08,   2.15443469e-08,   2.47707636e-08,
         2.84803587e-08,   3.27454916e-08,   3.76493581e-08,
         4.32876128e-08,   4.97702356e-08,   5.72236766e-08,
         6.57933225e-08,   7.56463328e-08,   8.69749003e-08,
         1.00000000e-07,   1.14975700e-07,   1.32194115e-07,
         1.51991108e-07,   1.74752840e-07,   2.00923300e-07,
         2.31012970e-07,   2.65608778e-07,   3.05385551e-07,
         3.51119173e-07,   4.03701726e-07,   4.64158883e-07,
         5.33669923e-07,   6.13590727e-07,   7.05480231e-07,
         8.11130831e-07,   9.32603347e-07,   1.07226722e-06,
         1.23284674e-06,   1.41747416e-06,   1.62975083e-06,
         1.87381742e-06,   2.15443469e-06,   2.47707636e-06,
         2.84803587e-06,   3.27454916e-06,   3.76493581e-06,
         4.32876128e-06,   4.97702356e-06,   5.72236766e-06,
         6.57933225e-06,   7.56463328e-06,   8.69749003e-06,
         1.00000000e-05,   1.14975700e-05,   1.32194115e-05,
         1.51991108e-05,   1.74752840e-05,   2.00923300e-05,
         2.31012970e-05,   2.65608778e-05,   3.05385551e-05,
         3.51119173e-05,   4.03701726e-05,   4.64158883e-05,
         5.33669923e-05,   6.13590727e-05,   7.05480231e-05,
         8.11130831e-05,   9.32603347e-05,   1.07226722e-04,
         1.23284674e-04,   1.41747416e-04,   1.62975083e-04,
         1.87381742e-04,   2.15443469e-04,   2.47707636e-04,
         2.84803587e-04,   3.27454916e-04,   3.76493581e-04,
         4.32876128e-04,   4.97702356e-04,   5.72236766e-04,
         6.57933225e-04,   7.56463328e-04,   8.69749003e-04,
         1.00000000e-03])
```

```python
>>> plt.semilogy(10*log10(SNR),C*1e9/400e6)
[<matplotlib.lines.Line2D at 0x7f9c92782a10>]
```

```python
>>> plt.plot(rg,linewidth=4)
>>> plt.plot(abs(Hc.Hcal.y[0,0,:]))
>>> plt.plot(abs(Hc.Hcal.y[0,1,:]))
>>> plt.plot(abs(Hc.Hcal.y[0,2,:]))
>>> plt.plot(abs(Hc.Hcal.y[1,0,:]))
>>> plt.plot(abs(Hc.Hcal.y[1,1,:]))
>>> plt.plot(abs(Hc.Hcal.y[1,2,:]))
>>> plt.plot(abs(Hc.Hcal.y[2,0,:]))
>>> plt.plot(abs(Hc.Hcal.y[2,1,:]))
>>> plt.plot(abs(Hc.Hcal.y[2,2,:]))
[<matplotlib.lines.Line2D at 0x7f9c924caa10>]
```

```python
>>> rho2,CWf = Hc.WFcapacity(Pt=Pt)
>>> rho2=np.sum(rho2*dfGHz,axis=0)
>>> rho2=np.sum(rho2,axis=0)
```

```python
>>> CB = np.sum(CB,axis=0)
>>> CWf = np.sum(CWf,axis=0)
```

```python
>>> CB.shape
(100,)
```

```python
>>> kB = 1.3806488e-23
>>> N0 = kB*273
>>> B  = 400e6
>>> Pb = N0*B
>>> #H =  rg[None,:]
... H = Hc.Hcal.y[0,0,:][None,:]
>>> SNR = Pt[:,None]*H*np.conj(H)/(N0*B)
>>> Csiso = np.log(1+SNR)/np.log(2)
```

```python
>>> Csiso=np.sum(Csiso,axis=1)*dfGHz
```

```python
>>> SNR = np.real(np.sum(SNR,axis=1)*dfGHz)
```

```python
>>> plt.semilogx(Pt,CB*1e9/400e6,label='BLAST',linewidth=4,alpha=0.5)
>>> plt.semilogx(Pt,CWf*1e9/400e6,label='WF',linewidth=4,alpha=0.5)
>>> plt.semilogx(Pt,Csiso*1e9/400e6,label='SISO',linewidth=4,alpha=0.5)
>>> plt.semilogx(Pt,C*1e9/400e6,label='SISO function',linewidth=1,alpha=0.5)
>>> plt.xlabel('Pt')
>>> plt.ylabel('Spectral Efficiency (bit/s/Hz)')
>>> plt.legend(loc='best')
/home/uguen/anaconda/lib/python2.7/site-packages/numpy/core/numeric.py:462: ComplexWarning: Casting complex values to real discards the imaginary part
  return array(a, dtype, copy=False, order=order)
```

```python
>>> alpha1 = CWf/Csiso
>>> alpha2 = CB/Csiso
```

```python
>>> alpha1[-1]
(3.4454854741006731-2.7796736825509248e-19j)
```

```python
>>> plt.semilogx(Pt,alpha1,label='Water Filling')
>>> plt.semilogx(Pt,alpha2,label='Blast')
>>> plt.legend(loc='best')
```

```python
>>> Hc.Hcal.y.shape
(8, 4, 1601)
```

```python
>>> CBF = Hc.BFcapacity(Pt=Pt)
```

```python
>>> CB2.shape
```

```python
>>> CBF[1].shape2
(1601, 4)
```

```python
>>> HdH1,U1,S1,V1=Hc.transfer()
```

```python
>>> H1 = Hc.Hcal.y[0,0,:]
```

```python
>>> Hc.normalize()
```

```python
>>> Pt=1e-3
>>> Tp = 273
>>> kB = 1.03806488e-23
>>> # N0 ~ J ~ W/Hz ~ W.s
... N0 = kB*Tp
```

```python
>>> H1 = Hc.Hcal.y
>>> H1c = np.conj(Hc.Hcal.y)
```

```python
>>> #Hsiso = Hc.Hcal[0,:,:]
... HH = H1*H1c
>>> CSISO1 = np.log(1+(Pt*np.real(HH)/N0))/np.log(2)
>>> CSISO3 = np.log(1+np.sum(np.sum(Pt*np.real(HH)/N0,axis=0),axis=0))/np.log(2)
```

```python
>>> CSISO3.shape
(1601,)
```

```python
>>> fGHz = Hc.Hcal.x
>>> df = fGHz[1]-fGHz[0]
```

```python
>>> C1=sum(CSISO*df,axis=2)
>>> C3=sum(CSISO3*df)
```

```python
>>> sum(C3)
```

```python
>>> plt.imshow(C,interpolation='nearest')
>>> plt.colorbar()
```

```python
>>> mean(C)*1e9/400e6
44.017021311673858
```

```python
>>> CSISO2 = np.log(1+(Pt*np.real(Hc.rg*Hc.rg)/N0))/np.log(2)
```

```python
>>> C2=sum(CSISO2*df)
>>> print C2
17.6918916897
```

```python
>>> plt.plot(Hc.rg)
[<matplotlib.lines.Line2D at 0x7fb7c85a5810>]
```

```python
>>> Hc.normalize
True
```

```python
>>> HdH2,U2,S2,V2=Hc.transfer()
```

```python
>>> plt.subplot(211)
>>> plt.semilogy(S1[:,0])
>>> plt.semilogy(S1[:,1])
>>> plt.semilogy(S1[:,2])
>>> plt.semilogy(S1[:,3])
>>> plt.subplot(212)
>>> plt.semilogy(S2[:,0])
>>> plt.semilogy(S2[:,1])
>>> plt.semilogy(S2[:,2])
>>> plt.semilogy(S2[:,3])
>>> plt.semilogy(sum(S2,axis=1))
[<matplotlib.lines.Line2D at 0x7fb7c9187410>]
```

```python
>>> plt.subplot(211)
>>> plt.plot(S1[:,0])
>>> plt.plot(S1[:,1])
>>> plt.plot(S1[:,2])
>>> plt.plot(S1[:,3])
>>> plt.subplot(212)
>>> plt.plot(S2[:,0])
>>> plt.plot(S2[:,1])
>>> plt.plot(S2[:,2])
>>> plt.plot(S2[:,3])
>>> plt.plot(sum(S2,axis=1))
[<matplotlib.lines.Line2D at 0x7fb7c8b47190>]
```

```python
>>> plt.plot(S2[:,0])
[<matplotlib.lines.Line2D at 0x7fb7c88462d0>]
```

```python
>>> W1 = U1[:,0,:]
>>> W2 = U2[:,0,:]
```

plt.imshow(abs(W2),interpolation='nearest')
plt.axis('tight')
plt.colorbar()

```python
>>> plt.imshow(abs(W1),interpolation='nearest')
>>> plt.axis('tight')W2 = U2[:,0,:]
>>> plt.colorbar()
```

```python
>>> plt.imshow(abs(W2),interpolation='nearest')
>>> plt.axis('tight')
>>> plt.colorbar()
```

```python
>>> W1-W2
array([[  3.33066907e-16 +2.77555756e-17j,
         -1.99840144e-15 +0.00000000e+00j,
         -8.32667268e-15 -1.11022302e-15j,
          2.60902411e-15 -2.30926389e-14j],
       [  0.00000000e+00 -6.93889390e-18j,
         -5.55111512e-16 -2.22044605e-16j,
          1.72084569e-15 -4.82947016e-15j,
         -2.33146835e-15 -6.99440506e-15j],
       [  4.44089210e-16 +1.73472348e-17j,
         -2.22044605e-16 +3.33066907e-16j,
         -3.96904731e-15 -5.44009282e-15j,
          2.44249065e-15 +2.40918396e-14j],
       ..., 
       [ -7.77156117e-16 +1.38777878e-17j,
         -1.77635684e-15 +5.55111512e-16j,
         -1.13242749e-14 +1.22124533e-15j,
         -8.99280650e-15 +5.88418203e-15j],
       [  0.00000000e+00 +2.77555756e-17j,
          7.77156117e-16 +1.22124533e-15j,
         -1.20459198e-14 -9.32587341e-15j,
         -5.55111512e-15 +3.27515792e-15j],
       [  5.55111512e-16 -4.16333634e-17j,
          3.10862447e-15 -7.77156117e-16j,
         -4.38538095e-15 -7.88258347e-15j,
          8.88178420e-16 -1.77635684e-15j]])
```

```python
>>> try: 
...     del tBlast
...     del twf1
...     del twf2
...     del ta.shape
>>> except:
...     pass
>>> tBlast=[]
>>> twf = []
>>> ta = []
>>> Hc =MIMO('S00_11_V.csv','/MIMO-FB/12-9-2013/S00/')
>>> for alpha in range(50):
...     print alpha
...     Hc.Hcal.y = Hc.Hcal.y*0.8
...     ta.append(0.8**(alpha+1))
...     M1,detM1,ldM1,C1,C2,sv1 = Hc.Bcapacity(Pt=1e-3)
...     CBlast = sum(C1)
...     tBlast.append(CBlast)
...     Cw,Q=Hc.WFcapacity(Pt=1e-3)
...     Cwf = sum(Cw)
...     twf.append(Cwf)
>>> ta = np.array(ta)
>>> twf = np.array(twf)
0
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
```

```python
>>> plt.loglog(1/ta[::-1],tBlast[::-1],label='Blast')
>>> plt.loglog(1/ta[::-1],twf[::-1],label='WF sum vp')
>>> plt.legend(loc='best')
>>> plt.xlabel('1/alpha')
>>> plt.ylabel('Capacity(GB/s')
```

```python
>>> sum(C2)
19.56534523138852
```

```python
>>> Hc.Hcal.y = Hc.Hcal.y
```

```python
>>> M1,detM1,ldM1,C1,C2,sv1 = Hc.Bcapacity(Pt=1e-3)
>>> print "BLAST capacity (Gb/s) : ",np.sum(C1)
BLAST capacity (Gb/s) :  19.5653452314
```

```python
>>> C1,C2,Q=Hc.WFcapacity(Pt=1e-3)
0.001
```

```python
>>> sum(C1)
19.563895575988624
```

```python
>>> sum(C2)
19.563911597666468
```

```python
>>> Q.shape
(1601, 4)
```

```python
>>> plt.imshow(Q,interpolation='nearest')
>>> plt.axis('tight')
>>> plt.colorbar()
```

```python
>>> plt.imshow(abs(u[0][0,:,:]))
```

```python
>>> plt.semilogy(sv1[:,0])
>>> plt.semilogy(sv1[:,1])
>>> plt.semilogy(sv1[:,2])
>>> plt.semilogy(sv1[:,3])
[<matplotlib.lines.Line2D at 0x7f9318ada550>]
```

```python
>>> kB = 1.3806488e-23
>>> N0 = kB*273
>>> Pt = 1e-3
>>> B  = 400e6
>>> Pb = N0
```

$\mathbf{H}^{\dagger}\mathbf{H}$ is a square Hermitian matrix which can be factorized as
$\mathbf{H}^{\dagger}\mathbf{H} = \mathbf{U}(f)\mathbf{\Lambda}(f)\mathbf{U}^{\dagger}(f)$

# Channel matrix power normalization

Let $\mathbf{H}_m(f)$ be the $(n_R \times n_T)$ complex matrix which corresponds a possibly measured matrix between transmitted and received symbols at frequency $f$.

$$\mathbf{y}(f) = \mathbf{\tilde{H}}(f) \mathbf{x}(f) + \mathbf{n}(f)$$

The channel matrix can be normalized as follows :

$$\mathbf{\tilde{H}}(f) = \sqrt{g(f)} \mathbf{H}(f)$$

This matrix normalization is common and can be fould in \cite{verdu12}. 
This normalization is interesting because it allows to distinguish between 2 different time variability which should be different for practical situations ( to be verified experimentally) . The partial CSI which could be used for the beamforming could be obtained only from the normalized channel matrix $\mathbf{H}$ and could have a much longer lifetime w.r.t to the Coherence time.

An idea would be to define two different coherence time based either on time variability of $g(f,t)$ or time variability of $\mathbf{H}(f,t)$

Where

$$g(f) = \frac{\mathrm{Tr}(\mathbf{\tilde{H}}(f)^{\dagger}\mathbf{\tilde{H}}(f))}{n_T n_R}$$

$$\mathrm{Tr}[\mathbf{H}(f)\mathbf{H}^{\dagger}(f)] 
= \frac{\mathrm{Tr}(\mathbf{\tilde{H}}\mathbf{\tilde{H}}^{\dagger})}{g(f)}= 
n_T n_R$$

#MIMO channel capacity

## No CSI at the transmitter

We assume that we have at our disposal a finite amount of power $P_t$ whose the best usage is seek it term of maximising the channel capacity.

The channel capacity in that case is determined as. This expression is assuming that $\mathbf{Q}=\mathbf{U}(f)\mathbf{D}(f)\mathbf{U}^{\dagger}(f)$:

$$C_{BLAST} = \int_{B} \log_2 \textrm{det}(\mathbf{I}_t + \frac{Pt }{N_0 B N_t} \mathbf{H}^{\dagger}(f)\mathbf{H}(f))df \;\; \textrm{(bit/s)}$$

$$C = \int_{B} \log_2 \textrm{det}(\mathbf{I}_t + \frac{E(f) }{N_0 N_t} \mathbf{H}^{\dagger}(f)\mathbf{H}(f))df \;\; \textrm{(bit/s)}$$

But if the power spectral density (homogeneous to an Energy) of the transmitted signal is assumed constant $E(f) = E$ and thus $ P_t = EB$

$$C_{BLAST} = \int_{B} \log_2 \textrm{det}(\mathbf{I}_t + \frac{Pt }{N_0 B N_t} \mathbf{H}^{\dagger}(f)\mathbf{H}(f))df \;\; \textrm{(bit/s)}$$

If we assume now that the channel tensor $\mathbf{H}(f)$ is known at the transmitter and in particular the eigenvalues of $\mathbf{H}^{\dagger}(f)\mathbf{H}(f)$

$\mathbf{\Lambda}(f)$ is a tensor $(f,N_t,N_t)$
$$\mathbf{\Lambda}(f) = \mathbf{diag}([\lambda_0(f) \ldots,\lambda_{N_t}(f)])$$

$$\mathbf{P}(f,\mu) = \mathbf{diag}([max(\mu-\frac{1}{\lambda_0(f)},0),\ldots, max(\mu-\frac{1}{\lambda_{N_t}(f)},0)])$$

or in a more compact form
$$\mathbf{P}(f,\mu) = \max(\mu-\mathbf{\Lambda}^{-1}(f),0)$$

$$C_w(\mu) = \int_{B} \log_2 \textrm{det}(\mathbf{I}_t + \frac{1}{N_0} \mathbf{U}(f)\mathbf{E}(f,\mu)\mathbf{\Lambda}(f)\mathbf{U}^T(f)df \;\; \textrm{(bit/s)}$$

We can shift the $\mathbf{P}(f,\mu)$ matrix at the right of the unitary matrix $\mathbf{U}$ because we can show that
$\mathbf{diag}(\mathbf{P}(f,\mu)\mathbf{U}(f))=\mathbf{diag}(\mathbf{U}(f)\mathbf{P}(f,\mu))$ and as the next matrix on the right \mathbf{\Lambda}(f) is also diagonal this doesn't modify the whole expression.

$$C_w(\mu) = \int_{B} \log_2 \textrm{det}(\mathbf{I}_t + \frac{1}{N_0} \mathbf{U}(f)\mathbf{P}(f,\mu)\mathbf{\Lambda}(f)\mathbf{U}^{\dagger}(f)df \;\; \textrm{(bit/s)}$$

## SVD and eigenvalue decomposition of $\mathbf{H}^{\dagger}\mathbf{H}$

```python
>>> HHT,U,S,V=Hc.transfer()
```

```python
>>> HHT.shape
(1601, 4, 4)
```

```python
>>> U[0,:,:]
array([[-0.39471602 +1.38777878e-17j,  0.76088052 -1.66533454e-16j,
        -0.26801053 -3.33066907e-15j, -0.43980729 -1.66533454e-14j],
       [-0.24319345 -3.58525932e-01j, -0.16418990 +2.78181914e-01j,
         0.36435070 +4.83905988e-01j, -0.28782235 +5.08148074e-01j],
       [-0.01070260 -5.37556157e-01j,  0.01184301 -8.32378671e-02j,
        -0.54633748 -2.21886374e-01j,  0.36302217 +4.73652542e-01j],
       [ 0.59745206 -1.02397423e-01j,  0.54627334 +1.05990478e-01j,
         0.12828912 +4.43934249e-01j,  0.33069484 +4.74057677e-03j]])
```

```python
>>> U[0,:,:]-np.conj(V[0,:,:].T)
array([[  1.44328993e-15 +1.38777878e-17j,
         -1.11022302e-15 -1.66533454e-16j,
         -2.72004641e-15 -3.33066907e-15j,
          2.22044605e-16 -1.66533454e-14j],
       [ -3.88578059e-16 +3.88578059e-16j,
          1.38777878e-15 +1.83186799e-15j,
         -4.88498131e-15 +3.88578059e-16j,
         -2.00395256e-14 -9.54791801e-15j],
       [ -1.73472348e-18 +0.00000000e+00j,
         -3.05311332e-16 +6.80011603e-16j,
          0.00000000e+00 +2.96984659e-15j,
         -1.23789867e-14 +9.76996262e-15j],
       [  2.22044605e-16 +2.08166817e-16j,
          2.22044605e-16 +2.22044605e-16j,
         -7.21644966e-16 +1.16573418e-15j,
         -3.83026943e-15 +9.57567359e-15j]])
```

### Reconstruction from SVD decomposition

```python
>>> H = Hc.Hcal.y
>>> Hd  = np.conj(Hc.Hcal.y.swapaxes(0,1))
```

```python
>>> Hd.shape
(4, 8, 1601)
```

```python
>>> HdH=np.einsum('ijk,jlk->ilk',Hd,H)
>>> HHd=np.einsum('ijk,jlk->ilk',H,Hd)
```

```python
>>> HdH.shape
(4, 4, 1601)
```

```python
>>> HHd.shape
(8, 8, 1601)
```

```python
>>> T1 = np.real(trace(HdH))
>>> T2 = np.real(trace(HHd))
>>> Hc.Hcal.y.shape
(8, 4, 1601)
```

```python
>>> HdH.shape
(4, 4, 1601)
```

```python
>>> T1.shape
(1601,)
```

```python
>>> rho0 = sum(Trace)*(fGHz[1]-fGHz[0])*1e9
```

```python
>>> 10*log10(rho0)
58.694040524393159
```

```python
>>> plt.plot(np.sqrt(T1/32.))
>>> plt.plot(np.abs(Hc.Hcal.y[0,0,:]))
[<matplotlib.lines.Line2D at 0x7f56903a67d0>]
```

```python
>>> Hn = H/np.sqrt(Trace[None,None,:])
```

```python
>>> SNRMIMO = 10*np.log10(Pt*Trace/(N0*B))
```

```python
>>> plt.plot(fGHz,SNRMIMO)
[<matplotlib.lines.Line2D at 0x7f5690267c10>]
```

```python
>>> Hn.shape
(8, 4, 1601)
```

L'expérimentation numérique ci-dessous est surprenante elle montre que la trace de la matice de transfert n'est autre que l'accumulation des modules carrés des 32 canaux SISO de la matice $\mathbf{H}$

```python
>>> del TY
>>> for m in range(8):
...     for n in range(4): 
...         Y = np.abs(Hc.Hcal.y[m,n,:])**2
...         try:
...             TY = TY+Y
...         except:
...             TY = Y
>>> 
>>> plt.figure(figsize=(10,5))
>>> plt.plot(TY/32,'b',alpha=0.6,linewidth=4,label='Average of the 32 Squared Channel Transfer function')
>>> plt.plot(T1/32,'r',label='Trace of HHd/32')
>>> plt.plot(np.abs(Hc.Hcal.y[3,1,:])**2,'g',label='One among 32 channel')
>>> plt.plot(np.abs(Hc.Hcal.y[1,1,:])**2,'m',label='One among 32 channel')
>>> plt.legend(loc='best')
```

```python
>>> S.shape
(1601, 4)
```

expanding S

```python
>>> I4 = np.eye(4)
>>> Se=S[:,:,None]*I4[None,:,:]
```

```python
>>> Se.shape
(1601, 4, 4)
```

```python
>>> US=np.einsum('ijk,ikl->ijl',U,Se)
>>> USV=np.einsum('ijk,ikl->ijl',US,V)
```

So we have been able to reconstruct the transfer matrix from its decompostion in singular value decomposition in multi dimensional arrays.

```python
>>> sum(USV-HHT)
(-3.8164148082061271e-16+3.2451664976200339e-15j)
```

```python
>>> plt.imshow(abs(M1[:,2,:]))
>>> plt.axis('tight')
(-0.5, 3.5, 1600.5, -0.5)
```

This method also return the positive eigenvalues $\lambda_k(f)$ of the transfer matrix $\mathbf{H}^{\dagger}\mathbf{H}$

In [MIMO Capacity with Channel State information at the Transmitter]

```python
>>> SNR = 10*np.log10(Pt*T1/(32*N0*B))
```

```python
>>> plt.plot(SNR)
[<matplotlib.lines.Line2D at 0x7f568ff82890>]
```

```python
>>> I4 = np.eye(4)
>>> G = sv1[:,:,None]*I4[None,:,:]
>>> mu = 1e6
```

```python
>>> sv1.shape
(1601, 4)
```

The noise spectral density is

```python
>>> B
400000000.0
```

```python
>>> 10*log10(N0*B)+30
-88.216941418646954
```

```python
>>> pb = N0*(fGHz[1]-fGHz[0])*1e9*ones(1601)[:,None]*ones(4)
```

```python
>>> pb.shape
(1601, 4)
```

```python
>>> sum(pb)
6.0344431296260165e-12
```

```python
>>> plt.imshow(log(sv1),interpolation='nearest')
>>> plt.axis('tight')
(-0.5, 3.5, 1600.5, -0.5)
```

```python
>>> G = pb/sv1
```

```python
>>> G
array([[  7.35522492e-13,   2.48215410e-11,   2.23535218e-10,
          1.08863931e-09],
       [  7.30662847e-13,   2.49989407e-11,   2.23749069e-10,
          1.16993486e-09],
       [  7.26031992e-13,   2.50890958e-11,   2.15571584e-10,
          1.17219008e-09],
       ..., 
       [  4.15553815e-13,   1.76386009e-11,   2.19735813e-10,
          3.59773153e-10],
       [  4.15595713e-13,   1.76525194e-11,   2.12327994e-10,
          3.59960768e-10],
       [  4.15290639e-13,   1.73907427e-11,   2.12974193e-10,
          3.60119479e-10]])
```

```python
>>> #pt = Pt/(4*1601*(fGHz[1]-fGHz[0])*1e9)*ones((1601,4))
... pt = (Pt/(4*1601))*ones((1601,4))
```

```python
>>> sum(pt)
0.001
```

```python
>>> sum(pb/sv1)
2.2013333492922998e-06
```

```python
>>> SNR=10*log10((pt/(pb/sv1)))
```

```python
>>> plt.imshow(SNR,interpolation='nearest')
>>> plt.axis('tight')
>>> plt.colorbar()
```

```python
>>> #mu = 1.002*Pt/(1601*4*(fGHz[1]-fGHz[0])*1e9)
... mu = Pt/(1601*4)
>>> Q = np.maximum(0,mu-pb/sv1)
```

```python
>>> Peff = sum(Q)
>>> print Peff
0.000997798666651
```

```python
>>> DeltaP=Pt-Peff
>>> print DeltaP
2.20133334929e-06
```

```python
>>> Deltamu = DeltaP/(1601*4)
```

```python
>>> muopt = mu + Deltamu
```

```python
>>> print muopt
1.56496148243e-07
```

```python
>>> Qopt = np.maximum(0,muopt-pb/sv1)
```

```python
>>> Qopt.shape
(1601, 4)
```

```python
>>> sum(Qopt)
0.001
```

```python
>>> Q-Qopt
array([[ -1.37497398e-15,  -1.37497398e-15,  -1.37497398e-15,
         -1.37497398e-15],
       [ -1.37497398e-15,  -1.37497398e-15,  -1.37497398e-15,
         -1.37497398e-15],
       [ -1.37497398e-15,  -1.37497398e-15,  -1.37497398e-15,
         -1.37497398e-15],
       ..., 
       [ -1.37497398e-15,  -1.37497398e-15,  -1.37497398e-15,
         -1.37497398e-15],
       [ -1.37497398e-15,  -1.37497398e-15,  -1.37497398e-15,
         -1.37497398e-15],
       [ -1.37497398e-15,  -1.37497398e-15,  -1.37497398e-15,
         -1.37497398e-15]])
```

```python
>>> plt.imshow(Q-Qopt,interpolation='nearest')
>>> plt.axis('tight')
>>> plt.colorbar()
```

```python
>>> g=np.maximum(0,mu-pb/sv1)
>>> #M2,detM2,ldM2,C2,sv2 = Hc.capacity(Pt=Pt)
```

```python
>>> pb/sv1
array([[  2.94208997e-18,   9.92861639e-17,   8.94140871e-16,
          4.35455723e-15],
       [  2.92265139e-18,   9.99957630e-17,   8.94996277e-16,
          4.67973942e-15],
       [  2.90412797e-18,   1.00356383e-16,   8.62286335e-16,
          4.68876030e-15],
       ..., 
       [  1.66221526e-18,   7.05544035e-17,   8.78943252e-16,
          1.43909261e-15],
       [  1.66238285e-18,   7.06100778e-17,   8.49311976e-16,
          1.43984307e-15],
       [  1.66116256e-18,   6.95629709e-17,   8.51896770e-16,
          1.44047792e-15]])
```

```python
>>> Gn=G/sum(G)
```

```python
>>> tmu = np.logspace(5,7,100)
```

```python
>>> tC = []
>>> for mu in tmu:
...     g=np.maximum(0,mu-1/G)
...     Pt = g*1e-3/sum(g)
...     M,detM,ldM,Cf,sv = Hc.capacity(Pt=Pt)
...     C = sum(Cf)
...     tC.append(C)
```

```python
>>> 1e-3/(1601*4)
1.561524047470331e-07
```

```python
>>> Es,Cf,sv = Hc.capacity(Pt=1e-3)
```

```python
>>> sum(Cf)
```

```python
>>> ref=1e-3/(1601*4)*np.ones(1601)[:,None,None]*I4[None,:,:]
```

```python
>>> sum(ref)
0.001
```

```python
>>> Pt.shape
(1601, 4, 4)
```

```python
>>> E=Pt-ref
```

```python
>>> sum(E)
1.0587911840678754e-19
```

```python
>>> Es,Cf,sv = Hc.capacity(Pt=ref)
```

```python
>>> sum(Cf)
```

```python
>>> semilogx(tmu,tC)
>>> xlabel(r"$\mu$")
>>> ylabel("Mimo capacity")
```

```python
>>> 1e-3/(4*1601)
```

```python
>>> np.where(tC==max(tC))
```

```python
>>> tC[99]-sum(C1)
```

```python
>>> Es,Cf,sv = Hc.capacity(Pt=Pt)
```

```python
>>> sum(Cf)
```

```python
>>> plt.plot(tmu,tC)
```

```python
>>> Es2,C2,sv2 = Hc.capacity(Pt=Pt)
```

```python
>>> sum(C2)
```

```python
>>> sum(C1)
```

```python
>>> sum(Pt)*1000
```

```python
>>> plt.figure(figsize=(10,10))
>>> plot(Hc.Hcal.x,Es2)
>>> xlabel('Frequency (Ghz)')
```

```python
>>> A = np.random.randn(8,4)+1j*np.random.randn(8,4)
```

```python
>>> A.shape
(8, 4)
```

```python
>>> H1=np.dot(A,np.conj(A.T))
>>> H2=np.dot(np.conj(A.T),A)
```

```python
>>> trace(H1)
(61.468784899215152-5.4470317145671743e-16j)
```

```python
>>> trace(H2)
(61.468784899215152+4.0939474033052647e-16j)
```

```python
>>> sum(abs(A)**2)
61.468784899215152
```

```python
>>> import numpy.linalg as la
```

```python
>>> l,U=la.eig(H1)
```

```python
>>> L=np.diag(l)
```

```python
>>> np.dot(U,L),np.c
array([[  3.06655843e+01 +7.65741511e-16j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j],
       [  0.00000000e+00 +0.00000000e+00j,
          1.77837422e+01 +6.58989659e-16j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j],
       [  0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          9.54272285e+00 +3.05849165e-16j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j],
       [  0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          3.47673549e+00 +1.10895497e-16j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j],
       [  0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          2.70548202e-15 +1.81282286e-16j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j],
       [  0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          6.84704496e-16 +1.18374590e-15j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j],
       [  0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          2.66331914e-17 -1.17889970e-15j,
          0.00000000e+00 +0.00000000e+00j],
       [  0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
         -1.15276621e-15 +2.59285372e-16j]])
```

```python
>>> ones((2,2))
array([[ 1.,  1.],
       [ 1.,  1.]])
```

```python
>>> It  = np.eye(4)
>>> Se  = S[:,:,None]*It[None,:,:]
```

```python
>>> Se.shape
(1601, 4, 4)
```

```python
>>> la.det(Se)
array([  1.77455244e-19,   1.64885305e-19,   1.71282574e-19, ...,
         1.36057503e-18,   1.40605826e-18,   1.42331044e-18])
```

```python
>>> la.det(Se[0,:,:])
1.7745524407309276e-19
```

```python
>>> S.shape
(1601, 4)
```

```python
>>> import scipy as sp
```

```python
>>> import numpy as np
```

```python
>>> np.trace?
```

```python
>>> np.trace?
```

```python
>>> H3 = Hc.Hcal.y
>>> H3 = H3.swapaxes(0,2)
```

```python
>>> U,S,V=la.svd(H3)
```

```python
>>> U.shape
(1601, 4, 4)
```

```python
>>> S.shape
(1601, 4)
```

```python
>>> V.shape
(1601, 8, 8)
```

```python
>>> Pt=np.random.rand(100)
```

```python
>>> np.maximum(0.5,Pt)
array([ 0.85851988,  0.96979777,  0.5       ,  0.52189233,  0.5       ,
        0.80216662,  0.5       ,  0.5       ,  0.5       ,  0.67235463,
        0.5       ,  0.5       ,  0.5       ,  0.5       ,  0.5       ,
        0.5       ,  0.5       ,  0.5       ,  0.77906121,  0.5       ,
        0.55608479,  0.97391911,  0.5       ,  0.70418788,  0.5       ,
        0.5       ,  0.59350424,  0.5       ,  0.59303855,  0.5       ,
        0.73461344,  0.5       ,  0.5       ,  0.5       ,  0.5       ,
        0.5       ,  0.75436453,  0.71910043,  0.98449101,  0.53388298,
        0.5       ,  0.86289408,  0.5       ,  0.68152341,  0.5       ,
        0.5       ,  0.7062983 ,  0.5       ,  0.74917636,  0.5       ,
        0.89625822,  0.54048752,  0.99481911,  0.74024208,  0.696227  ,
        0.98301098,  0.65457003,  0.5       ,  0.5       ,  0.50439848,
        0.5       ,  0.9307032 ,  0.5       ,  0.5       ,  0.5       ,
        0.5       ,  0.66533338,  0.5       ,  0.5       ,  0.5       ,
        0.56757035,  0.96646346,  0.87779665,  0.5       ,  0.78236851,
        0.5       ,  0.5       ,  0.57152478,  0.5       ,  0.5       ,
        0.66138117,  0.64482097,  0.65604395,  0.90231806,  0.6511418 ,
        0.5       ,  0.92777895,  0.86937716,  0.78514714,  0.5       ,
        0.751361  ,  0.5       ,  0.5       ,  0.5       ,  0.5       ,
        0.7353576 ,  0.5       ,  0.72399649,  0.98889226,  0.59960217])
```

```python
>>> a = r_[1,5,6,0]
```

```python
>>> b = r_[5,0,7,-1]
```

```python
>>> usup=where(a>b)[0]
```

```python
>>> usup
array([1, 3])
```

```python
>>> a[usup]
array([5, 0])
```

```python

```
