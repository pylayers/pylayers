### WHERE1 UWB Measurement campaign M1

```python
>>> %matplotlib inline
>>> from pylayers.measures.mesuwb import *
>>> from pylayers.gis.layout import *
>>> from pylayers.simul.link import *
>>> from pylayers.signal.waveform import *
```

First of all, we load the Layout of the environment. If the Layout associated graphs have already been built, one can load them with the `dumpr()` method.

```python
>>> L=Layout('WHERE1.ini')
>>> L.dumpr()
```

```python
>>> try:
...     del td1
...     del td2
...     del td3
...     del td4
...     del te1
...     del te2
...     del te3
...     del te4
...     del tt1
...     del tt2
...     del tt3
...     del tt4
>>> except:
...     pass
```

The UWB measure from the WHERE1 project are handled in the class `UWBMeasure`

```python
>>> K=UWBMeasure(10)
```

The delays associated with the 4 links are obatained directly

```python
>>> K.de
array([ 46.19349565,  23.96692223,  71.6705632 ,  55.09300165])
```

As well as more information if needed.

```python
>>> K.info()
Date_Time : [u'31-Jul-2008 08:17:50']
Tx_height : [u'120cm']
Tx_position : [u'P010']
Tx :  [-25.3654  12.304    1.2   ]
------Tx1 ------
delays     (ns): 46.1934956511
range  (meters): 13.8580486953
visibility     : NLOS2
angular (degree)  : 2.80776139034
LQI Meth1 12.7383607831  (dB)
LQI Meth2 -1.2959656826  (dB)
------Tx2 ------
delays     (ns): 23.9669222261
range  (meters): 7.19007666782
visibility     : NLOS2
angular (degree)  : 3.55280136793
LQI Meth1 21.4526880977  (dB)
LQI Meth2 5.56578610451  (dB)
------Tx3 ------
delays     (ns): 71.6705632005
range  (meters): 21.5011689601
visibility     : NLOS2
angular (degree)  : 2.98073648115
LQI Meth1 11.5021659794  (dB)
LQI Meth2 4.23241262908  (dB)
------Tx4 ------
delays     (ns): 55.0930016539
range  (meters): 16.5279004962
visibility     : NLOS
angular (degree)  : 3.31710006758
LQI Meth1 28.2789637953  (dB)
LQI Meth2 4.89658815244  (dB)
```

```python
>>> ### Simulation section
... fig=plt.figure(figsize=(10,5))
>>> f,a = K.show(delay=K.de)
```

```python
>>> #toa1 = K.toa_new()
... #toa2 = K.toa_th(1e-1,1e-1)
```

```python
>>> K.tau_Emax()
array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
       [ 0.,  0.,  0., ...,  0.,  0.,  0.],
       [ 0.,  0.,  0., ...,  0.,  0.,  0.],
       [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
```

```python
>>> np.vstack((K.rx))
array([[  0.    ,   0.    ,   1.2   ],
       [-12.2724,   7.7632,   1.2   ],
       [-18.7747,  15.178 ,   1.2   ],
       [ -4.1418,   8.8603,   1.2   ],
       [ -9.0914,  15.1899,   1.2   ]])
```

The code below reads data from the M1-WHERE2 measurement campaign.

```python
>>> M=UWBMeasure(15)
```

```python
>>> for k in range(300):
...     try:
...         M  = UWBMeasure(k)
...     except:
...         break
...     tx = M.tx
...     D  = M.rx-tx[np.newaxis,:]
...     D2 = D*D
...     dist = np.sqrt(np.sum(D2,axis=1))[1:]
...     Emax = M.Emax()
...     Etot = M.Etot()[0]
...     try:
...         td1 = np.hstack((td1,dist[0]))
...         td2 = np.hstack((td2,dist[1]))
...         td3 = np.hstack((td3,dist[2]))
...         td4 = np.hstack((td4,dist[3]))
>>> 
...         te1 = np.hstack((te1,Emax[0]))
...         te2 = np.hstack((te2,Emax[1]))
...         te3 = np.hstack((te3,Emax[2]))
...         te4 = np.hstack((te4,Emax[3]))
>>> 
...         tt1 = np.hstack((tt1,Etot[0]))
...         tt2 = np.hstack((tt2,Etot[1]))
...         tt3 = np.hstack((tt3,Etot[2]))
...         tt4 = np.hstack((tt4,Etot[3]))
...         #tdist = np.hstack((tdist,dist))
...         #te = np.hstack((te,Emax))
...     except:
...         td1=np.array(dist[0])
...         td2=np.array(dist[1])
...         td3=np.array(dist[2])
...         td4=np.array(dist[3])
...         te1 =np.array(Emax[0])
...         te2 =np.array(Emax[1])
...         te3 =np.array(Emax[2])
...         te4 =np.array(Emax[3])
...         tt1 =np.array(Etot[0])
...         tt2 =np.array(Etot[1])
...         tt3 =np.array(Etot[2])
...         tt4 =np.array(Etot[3])
```

The IR-UWB applied waweform is available in the raw data structure and can be extracted as follow. This exracttion is important in order to proceeed to the ray tracing simulation with the same waveform as the one used in the measurement campaign.

```python
>>> td1
```

```python
>>> from pylayers.signal.bsignal import *
>>> s=M.RAW_DATA.tx[0]
>>> t=M.RAW_DATA.timetx[0]*1e9
>>> plt.plot(t,s)
>>> plt.xlabel('time (ns)')
>>> se=TUsignal(t,s)
```

```python
>>> te = t[1]-t[0]
>>> cs = np.cumsum(s*s)
>>> E = cs[-1]*te
>>> EdB = 10*np.log10(E*30)
>>> print EdB
>>> print E*30
>>> use =1/E
>>> print use
-10.2361907016
0.0947067492189
316.767286888
```

```python
>>> E2=se.Emax()
>>> print E2*30
>>> E2dB=10*np.log10(E2*30)
>>> print E2dB
```

```python
>>> se.plot(typ='v')
(<matplotlib.figure.Figure at 0x7f1da9c88410>,
 array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7f1da9da1c10>]], dtype=object))
```

```python
>>> td1
```

```python
>>> fig = plt.figure(figsize=(16,6))
>>> ax = fig.add_subplot(111)
>>> ax.semilogx(td1,te1+EdB,'.r',label='Rx1')
>>> ax.semilogx(td2,te2+EdB,'.b',label='Rx2')
>>> ax.semilogx(td3,te3+EdB,'.g',label='Rx3')
>>> ax.semilogx(td4,te4+EdB,'.c',label='Rx4')
>>> d = np.linspace(1,30,100)
>>> 
>>> LFS = -(32.4+20*np.log10(4)+20*np.log10(d))-4
>>> ax.semilogx(d,LFS)
>>> plt.legend()
>>> plt.grid()
```

```python
>>> fig = plt.figure(figsize=(16,6))
>>> ax = fig.add_subplot(111)
>>> ax.semilogx(td1,tt1+EdB,'.r',label='Rx1')
>>> ax.semilogx(td2,tt2+EdB,'.b',label='Rx2')
>>> ax.semilogx(td3,tt3+EdB,'.g',label='Rx3')
>>> ax.semilogx(td4,tt4+EdB,'.c',label='Rx4')
>>> d = np.linspace(1,30,100)
>>> 
>>> LFS = -(32.4+20*np.log10(4)+20*np.log10(d))-4
>>> ax.semilogx(d,LFS)
>>> plt.legend()
>>> plt.grid()
```

```python
>>> M.Etot()
```

On choisit un index pour Tx et Rx

```python
>>> #measure id
... tx_id = 100 #in M.valid_index
>>> rx_id = 2 #1,2,3,4
>>> M=UWBMeasure(tx_id)
>>> TX = M.tx
>>> RX =M.rx[rx_id]
```

```python
>>> TX
```

```python
>>> M.rx
```

```python
>>> L.showG('s',figsize=(8,4))
>>> plt.plot(TX[0],TX[1],'ob')
>>> plt.plot(RX[0],RX[1],'or')
>>> plt.title('TOF = '+ str(np.sqrt(np.sum((TX-RX)**2))/0.3))
```

```python
>>> TX
```

```python
>>> RX
```

```python
>>> #Lk = DLink()
... Lk = DLink(L=L,a=TX,b=RX,cutoff=4,verbose=False)
>>> Lk.Aa=Antenna('defant.vsh3')
>>> Lk.Ab=Antenna('defant.vsh3')
```

```python
>>> Lk.eval(force=['ray','Ct','H'],alg=5)
>>> #f,a = Lk.show(rays=True,labels=False)
```

```python
>>> #%timeit Lk.eval(force=True,alg=7,cutoff=3)
... #f,a = Lk.show(rays=True,labels=False)
```

```python
>>> Lk.R
```

```python
>>> #%timeit Lk.eval(force=True,alg=7,threshold=0.01)
... #f,a = Lk.show(rays=True,labels=False)
```

```python
>>> Lk.Si.keys()
```

```python
>>> U=Lk.R[4]['sig2d'][0]
```

```python
>>> print U.shape
```

```python
>>> s1 = U[:,:,0]
>>> print s1
```

```python
>>> from pylayers.antprop.signature import Signature
```

```python
>>> S=Signature(s1)
```

```python
>>> S
```

```python
>>> wav = Waveform(typ='W1compensate')
```

```python
>>> wav.show()
```

```python
>>> #ir = Lk.H.applywavB(wav.sfg)
```

```python
>>> Lk.H.isFriis
```

```python
>>> if Lk.H.isFriis:
...     ir = Lk.H.applywavB(wav.sf)
>>> else:
...     ir = Lk.H.applywavB(wav.sfg)
```

```python
>>> Lk.R.los
```

```python
>>> Lk.H.ak
```

```python
>>> Lk.H.taud
```

```python
>>> G=Lk.H.ift()
```

```python
>>> M.tdd.ch2.plot(typ='v')
>>> plt.xlim([10,130])
```

```python
>>> M.tx
```

```python
>>> M.rx
```

```python
>>> np.sqrt(np.sum((M.tx-M.rx[3,:])*(M.tx-M.rx[3,:]),axis=0))/0.3
```

```python
>>> Lk.H.ak
```

```python
>>> Lk.wav=wav
```

```python
>>> ir.plot(typ='v')
```

```python
>>> ir
```

```python
>>> fig = plt.figure(figsize=(10,7))
>>> ax1=fig.add_subplot(211)
>>> cmd='M.tdd.ch' + str(rx_id) + '.plot(typ=[\'l20\'],fig=fig,ax=ax1)'
>>> eval(cmd)
>>> ax2 = fig.add_subplot(212,sharex=ax1,sharey=ax1)
>>> #Lk.chanreal.plot(typ=['v'],fig=fig,ax=ax2)
... ir.plot(typ=['l20'],fig=fig,ax=ax2)
>>> plt.xlim(60,130)
```

```python

```
