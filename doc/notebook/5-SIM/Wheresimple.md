# WHERE1 UWB Measurement campaign M1

```python
>>> %matplotlib inline
>>> from pylayers.measures.mesuwb import *
```

```python
>>> from pylayers.measures.mesuwb import *
>>> from pylayers.gis.layout import *
>>> 
>>> from pylayers.simul.link import *
>>> from pylayers.signal.waveform import *
WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.
```

First of all we load the Layout of the environment. If the Layout associated graphs have already been built, one can load them with the `dumpr()` method.

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
...     del tt3### Simulation section
>>> 
...     del tt4
>>> except:
...     pass
```

```python
>>> K=UWBMeasure(2)
```

```python
>>> ### Simulation section
... fig=plt.figure(figsize=(15,8))
>>> K.show()
```

```python
>>> K.tau_Emax()
array([[ 72.61 ],
       [ 29.9  ],
       [ 87.835],
       [ 61.975]])
```

```python
>>> np.vstack((K.rx))
array([[  0.    ,   0.    ,   1.2   ],
       [-12.2724,   7.7632,   1.2   ],
       [-18.7747,  15.178 ,   1.2   ],
       [ -4.1418,   8.8603,   1.2   ],
       [ -9.0914,  15.1899,   1.2   ]])
```

The code below reads data from the M1-WHERE2 measurement campaign and store theme in arrays.

```python
>>> for k in range(300):
...     try:
...         M  = UWBMeasure(k)
...         tx = M.tx
...         D  = M.rx-tx[np.newaxis,:]
...         D2 = D*D
...         dist = np.sqrt(np.sum(D2,axis=1))[1:]
...         Emax = M.Emax()
...         Etot = M.Etot()[0]
...         try:
...             td1 = np.hstack((td1,dist[0]))
...             td2 = np.hstack((td2,dist[1]))
...             td3 = np.hstack((td3,dist[2]))
...             td4 = np.hstack((td4,dist[3]))
>>> 
...             te1 = np.hstack((te1,Emax[0]))
...             te2 = np.hstack((te2,Emax[1]))
...             te3 = np.hstack((te3,Emax[2]))
...             te4 = np.hstack((te4,Emax[3]))
>>> 
...             tt1 = np.hstack((tt1,Etot[0]))
...             tt2 = np.hstack((tt2,Etot[1]))
...             tt3 = np.hstack((tt3,Etot[2]))
...             tt4 = np.hstack((tt4,Etot[3]))
...             #tdist = np.hstack((tdist,dist))
...             #te = np.hstack((te,Emax))
...         except:
...             td1=np.array(dist[0])
...             td2=np.array(dist[1])
...             td3=np.array(dist[2])
...             td4=np.array(dist[3])
...             te1 =np.array(Emax[0])
...             te2 =np.array(Emax[1])
...             te3 =np.array(Emax[2])
...             te4 =np.array(Emax[3])
...             tt1 =np.array(Etot[0])
...             tt2 =np.array(Etot[1])
...             tt3 =np.array(Etot[2])
...             tt4 =np.array(Etot[3])
...     except:
...         pass
```

The IR-UWB applied waweform is available in the raw data and can be extracted as follow

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
```

```python
>>> E2=se.Emax()
>>> print E2*30
>>> E2dB=10*np.log10(E2*30)
>>> print E2dB
```

```python
>>> se.plot(typ='v')
```

```python
>>> fig = plt.figure(figsize=(18,6))
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
>>> plt.plot(te1,tt1,'.')
```

```python
>>> M.Etot()
```

```python
>>> #measure id
... tx_id = 100 #in M.valid_index
>>> rx_id = 3 #1,2,3,4
>>> M=UWBMeasure(tx_id)
>>> pylab
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
>>> L.showG('s',figsize=(16,8))
>>> plt.plot(TX[0],TX[1],'ob')
>>> plt.plot(RX[0],RX[1],'or')
>>> plt.title('TOF = '+ str(np.sqrt(np.sum((TX-RX)**2))/0.3))
```

```python
>>> Lk = DLink(L=L,a=TX,b=RX,cutoff=4,verbose=False)
>>> Lk.Aa=Antenna('defant.vsh3')
>>> Lk.Ab=Antenna('defant.vsh3')
```

```python
>>> Lk.eval(force=['ray','Ct','H'],alg=5)
>>> #f,a = Lk.show(rays=True,labels=False)
(array([  9.79138753e-05,   2.96104973e-04,   1.59308395e-04,
          2.27580461e-05,   1.34574073e-04,   4.82827749e-04,
          1.18095058e-03,   1.27977423e-04,   1.34909317e-04,
          5.48270185e-05,   9.81119862e-05,   3.61269565e-04,
          2.73004108e-04,   1.24319829e-04,   2.28615744e-04,
          2.76183582e-04,   2.32071461e-04,   2.86885128e-04,
          5.08065400e-04,   4.16943960e-05,   9.25585146e-05,
          3.67562947e-05,   3.26755640e-05,   7.18843013e-05,
          5.98546831e-05,   5.35179821e-05,   8.90165410e-05,
          3.64178931e-04,   3.64157003e-04,   1.49754933e-04,
          9.73558942e-05,   4.65239631e-05,   5.86669038e-05,
          2.05936803e-04,   9.48731849e-05,   2.36321979e-04,
          6.88923268e-04,   1.24565228e-04,   4.69513933e-05,
          2.31194688e-05,   2.94412022e-05,   2.94227854e-05,
          5.60825032e-05,   7.29054840e-05,   1.70818338e-04,
          1.70818449e-04,   5.34747295e-05,   4.74216436e-05,
          1.14993846e-04,   5.08387934e-05]),
 array([  67.6920456 ,   67.6920456 ,  116.7880677 ,  116.7880677 ,
          68.63195952,   89.0838292 ,   62.63975138,   69.09664151,
          69.67313591,   70.48894853,   75.87644272,   82.1700362 ,
          89.4423201 ,   89.8884232 ,   91.14669303,  102.17774607,
          63.14854276,   63.77882449,   64.67003199,   70.94146788,
          71.48668315,   71.48668315,   71.5030899 ,   76.29701541,
          76.81949336,   82.55855407,   83.04164527,   91.30130681,
          91.30130681,   91.49710187,   91.93323475,  102.49044732,
         102.87998732,  115.06282092,   65.1629729 ,   65.75514013,
          65.75514013,   65.7739541 ,   73.27135774,   73.27135774,
          78.46804802,   78.46804802,   84.56899461,   84.56899461,
          93.31516303,   93.31516303,  104.11672196,  104.11672196,
         115.34059458,  115.68687375]))
```

```python
>>> #%timeit Lk.eval(force=True,alg=7,cutoff=3)
... #f,a = Lk.show(rays=True,labels=False)
```

```python
>>> Lk.R
Rays3D
----------
8 / 4 : [0 1 2 3]
4 / 2 : [4 5]
5 / 10 : [ 6  7  8  9 10 11 12 13 14 15]
6 / 18 : [16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33]
7 / 16 : [34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49]
-----
ni : 310
nl : 670
```

```python
>>> #%timeit Lk.eval(force=True,alg=7,threshold=0.01)
... #f,a = Lk.show(rays=True,labels=False)
```

```python
>>> Lk.Si.keys()
[3, 4, 5, 6, 7, 8, 9, 10]
```

```python
>>> U=Lk.R[4]['sig2d'][0]
```

```python
>>> print U.shape
(2, 4, 2)
```

```python
>>> s1 = U[:,:,0]
>>> print s1
[[328 335  67  73]
 [  2   3   3   3]]
```

```python
>>> from pylayers.antprop.signature import Signature
```

```python
>>> S=Signature(s1)
```

```python
>>> S
[328 335  67  73]
[2 3 3 3]
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
True
```

```python
>>> if Lk.H.isFriis:
...     ir = Lk.H.applywavB(wav.sf)
>>> else:
...     ir = Lk.H.applywavB(wav.sfg)
```

```python
>>> Lk.R.los
False
```

```python
>>> Lk.H.ak
array([  9.79138753e-05,   2.96104973e-04,   1.59308395e-04,
         2.27580461e-05,   1.34574073e-04,   4.82827749e-04,
         1.18095058e-03,   1.27977423e-04,   1.34909317e-04,
         5.48270185e-05,   9.81119862e-05,   3.61269565e-04,
         2.73004108e-04,   1.24319829e-04,   2.28615744e-04,
         2.76183582e-04,   2.32071461e-04,   2.86885128e-04,
         5.08065400e-04,   4.16943960e-05,   9.25585146e-05,
         3.67562947e-05,   3.26755640e-05,   7.18843013e-05,
         5.98546831e-05,   5.35179821e-05,   8.90165410e-05,
         3.64178931e-04,   3.64157003e-04,   1.49754933e-04,
         9.73558942e-05,   4.65239631e-05,   5.86669038e-05,
         2.05936803e-04,   9.48731849e-05,   2.36321979e-04,
         6.88923268e-04,   1.24565228e-04,   4.69513933e-05,
         2.31194688e-05,   2.94412022e-05,   2.94227854e-05,
         5.60825032e-05,   7.29054840e-05,   1.70818338e-04,
         1.70818449e-04,   5.34747295e-05,   4.74216436e-05,
         1.14993846e-04,   5.08387934e-05])
```

```python
>>> Lk.H.taud
array([  67.6920456 ,   67.6920456 ,  116.7880677 ,  116.7880677 ,
         68.63195952,   89.0838292 ,   62.63975138,   69.09664151,
         69.67313591,   70.48894853,   75.87644272,   82.1700362 ,
         89.4423201 ,   89.8884232 ,   91.14669303,  102.17774607,
         63.14854276,   63.77882449,   64.67003199,   70.94146788,
         71.48668315,   71.48668315,   71.5030899 ,   76.29701541,
         76.81949336,   82.55855407,   83.04164527,   91.30130681,
         91.30130681,   91.49710187,   91.93323475,  102.49044732,
        102.87998732,  115.06282092,   65.1629729 ,   65.75514013,
         65.75514013,   65.7739541 ,   73.27135774,   73.27135774,
         78.46804802,   78.46804802,   84.56899461,   84.56899461,
         93.31516303,   93.31516303,  104.11672196,  104.11672196,
        115.34059458,  115.68687375])
```

```python
>>> G=Lk.H.ift()
```

```python
>>> M.tdd.ch3.plot(typ='v')
>>> plt.xlim([10,130])
(10, 130)
```

```python
>>> M.tx
array([-22.3797,  13.3897,   1.2   ])
```

```python
>>> M.rx
array([[  0.    ,   0.    ,   1.2   ],
       [-12.2724,   7.7632,   1.2   ],
       [-18.7747,  15.178 ,   1.2   ],
       [ -4.1418,   8.8603,   1.2   ],
       [ -9.0914,  15.1899,   1.2   ]])
```

```python
>>> np.sqrt(np.sum((M.tx-M.rx[3,:])*(M.tx-M.rx[3,:]),axis=0))/0.3
62.639751380717335
```

```python
>>> Lk.H.ak
array([  9.79138753e-05,   2.96104973e-04,   1.59308395e-04,
         2.27580461e-05,   1.34574073e-04,   4.82827749e-04,
         1.18095058e-03,   1.27977423e-04,   1.34909317e-04,
         5.48270185e-05,   9.81119862e-05,   3.61269565e-04,
         2.73004108e-04,   1.24319829e-04,   2.28615744e-04,
         2.76183582e-04,   2.32071461e-04,   2.86885128e-04,
         5.08065400e-04,   4.16943960e-05,   9.25585146e-05,
         3.67562947e-05,   3.26755640e-05,   7.18843013e-05,
         5.98546831e-05,   5.35179821e-05,   8.90165410e-05,
         3.64178931e-04,   3.64157003e-04,   1.49754933e-04,
         9.73558942e-05,   4.65239631e-05,   5.86669038e-05,
         2.05936803e-04,   9.48731849e-05,   2.36321979e-04,
         6.88923268e-04,   1.24565228e-04,   4.69513933e-05,
         2.31194688e-05,   2.94412022e-05,   2.94227854e-05,
         5.60825032e-05,   7.29054840e-05,   1.70818338e-04,
         1.70818449e-04,   5.34747295e-05,   4.74216436e-05,
         1.14993846e-04,   5.08387934e-05])
```

```python
>>> Lk.wav=wav
```

```python
>>> ir.plot(typ='v')
(<matplotlib.figure.Figure at 0x7f06ccd92350>,
 array([[<matplotlib.axes.AxesSubplot object at 0x7f071cdbc710>]], dtype=object))
```

```python
>>> ir
Usignal :  (9047,)  (9047,)
ax0 : 9047
```

```python
>>> fig = plt.figure(figsize=(10,7))
>>> ax1=fig.add_subplot(211)
>>> cmd='M.tdd.ch' + str(rx_id) + '.plot(typ=[\'v\'],fig=fig,ax=ax1)'
>>> eval(cmd)
>>> ax2 = fig.add_subplot(212,sharex=ax1,sharey=ax1)
>>> #Lk.chanreal.plot(typ=['v'],fig=fig,ax=ax2)
... ir.plot(typ=['v'],fig=fig,ax=ax2)
>>> plt.xlim(60,130)
(60, 130)
```
