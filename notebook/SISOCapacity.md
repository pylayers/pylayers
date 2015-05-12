# Mimo Capacity and Water filling

```python
>>> from pylayers.measures.mesmimo import *
>>> from pylayers.gis.layout import *
>>> from matplotlib.pylab import *
>>> %matplotlib inline
WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.
```

In a first step, we calculate the MIMO channel capacity,  assuming no knowledge of the channel state information (CSI) at the transmitter. For that purpose, we set as argument of the `capacity` method the scalar value $Pt$ of the total power  available over the $min(N_t,N_r)$ spatial channel and the whole frequency band $B$.

```python
>>> Hc =MIMO('S00_11_V.csv','/MIMO-FB/12-9-2013/S00/')
>>> rho,C1 = Hc.Bcapacity()
>>> print "Uniform capacity (Gb/s) : ",sum(C1)
Uniform capacity (Gb/s) :  19.5653452314
```

We define the constant

```python
>>> kB = 1.3806488e-23
>>> N0 = kB*273
>>> Pt = 1e-3
>>> B  = 400e6
>>> Pb = N0*B
```

## SISO Case

The Shannon capacity can be defined as

$$C = \int_{B} \log_2 (1 + SNR(f)) df $$

where   $SNR$ at the receiver is given by :

$$ SNR(f) = \frac{\rho|H(f)|^2 E(f)}{N_0} $$

where we assume that $\int_B |H(f)|^2 df =1$

with

$$\int_B E(f) df = P_t$$

is the energy is spread uniformly over the whole frequency band $E(f) = E$

$$P_t = E B$$

and then

$$ SNR(f) = \frac{\rho|H(f)|^2 P_t}{N_0 B} $$

$$C_u = \int_{B} \log_2 (1 + \frac{\rho P_t }{N_0 B}|H(f)|^2) df $$

####Water filling

It the channel is known at the transmitter we can try to maximimize the capacity in choosing the proper PSD of the transmitter. A Power spectral density being expressed in $W/Hz$ is fundamentally homogeneous to energy that is the reason of the notation $E(f)$

\begin{eqnarray} C_w = \max_{E(f)} \int_{B} \log_2 (1 + \frac{\rho|H(f)|^2 E(f)}{N_0}) df \\
                 \textrm{s.t} \int_B E(f) df = P_t
\end{eqnarray}

We express the Lagrangian as

$$ L = \int_{B} \log_2 (1 + \frac{\rho|H(f)|^2 E(f)}{N_0}) df - \lambda (\int_B E(f) df - P_t)$$

$$\frac{\partial L}{\partial E} =\frac{1}{\ln 2}\int_{B} \frac{\frac{\rho |H(f)|^2}{N_0}}{1 + \frac{\rho|H(f)|^2 E(f)}{N_0}} df  + \lambda B $$

Notice that at high SNR the integrand becomes independent of $|H(f)|^2$ which suggest that the optimal repartition of the power  will gradually become independent from the channel knowledge.

$$\frac{\partial L}{\partial E} = 0 \Rightarrow$$

$$\frac{1}{\ln 2}\int_{B} \frac{\frac{\rho |H(f)|^2}{N_0}}{1 + \frac{\rho|H(f)|^2 E(f)}{N_0}} df  - \lambda B = 0 $$

$$\frac{1}{\ln 2}\int_{B} \frac{\frac{\rho|H(f)|^2}{N_0}}{1 + \frac{\rho|H(f)|^2 E(f)}{N_0}} df  =  \lambda B $$

$$\frac{\rho|H(f)|^2}{N_0 + \rho|H(f)|^2 E(f)}=\lambda \ln2$$

$$  \rho|H(f)|^2 = \lambda \ln2 ( N_0 + \rho|H(f)|^2 E(f))$$

$$  E(f) = \frac{\frac{\rho|H(f)|^2}{\lambda \ln2} - N_0  } {\rho|H(f)|^2}$$

$$ E(f) = \frac{1}{\lambda \ln 2}  - \frac{N_0}{\rho|H(f)|^2}$$

$$ E(f,\mu) = (\mu  - \frac{N_0}{\rho|H(f)|^2})^+$$

$$ \tilde{E}(f,\mu) = P_t \frac{(\mu  - \frac{N_0}{\rho|H(f)|^2})^+}{\int_B (\mu  - \frac{N_0}{\rho|H(f)|^2})^+df} $$

$$ \tilde{E}(f,\mu) = \frac{P_t}{\alpha(\mu,\rho)} (\mu  - \frac{N_0}{\rho|H(f)|^2})^+$$

$$\frac{\rho|H(f)|^2}{N_0 + \rho|H(f)|^2 \tilde{E}(f,\mu,\rho)}$$

###  Practical example

+ We extract a SISO channel over the bandwidth 1.8 GHz, 2.2 GHz
+ We evaluate the channel capacity (without CSI) as a function of SNR for this channel
+ We evaluate the channel capacity (with CSI) as a function of SNR and $\mu$
+ Q : Is there a maximum of the capacity for a given combination of

```python
>>> Hc.Hcal.y.shape
(8, 4, 1601)
```

```python
>>> H2.shape
```

```python
>>> #Frequency base definitions
... fGHz = Hc.Hcal.x
>>> df = fGHz[1]-fGHz[0]
>>> for i in range(8):
...     for j in range(4):
...         # lets take a SISO channel from the MIMO matrix
...         H = Hc.Hcal.y[i,j,:]
...         #H = Hc.Hcal.y
...         #Hd  = np.conj(Hc.Hcal.y.swapaxes(0,1))
...         #HdH=np.einsum('ijk,jlk->ilk',Hd,H)
...         #T = np.real(trace(HdH))
...         #H2 =T / 32.
...         H2 = np.real(H*np.conj(H))
...         #H2 = H2-min(H2)+1e-7
...         #plt.plot(Hc.Hcal.x,H2)
...         plt.plot(fGHz,H2)
...         plt.xlabel('fGHz')
...         plt.ylabel(u'$|H(f)|^{2}$')
```

```python
>>> cu = np.log(1+(Pt*H2/N0))/np.log(2)
```

```python
>>> plt.plot(cu)
[<matplotlib.lines.Line2D at 0x7f9b54180b50>]
```

```python
>>> CU=np.sum(cu)*df
```

```python
>>> CU
17.703040830890171
```

## SNR at the receiver

```python
>>> SNR = 10*np.log10(Pt*H2/(N0*B))
>>> plt.plot(fGHz,SNR)
>>> plt.ylabel('SNR (dB)')
>>> plt.xlabel('f (GHz)')
```

```python
>>> # Channel normalization
... Hn = H2/(df*np.sum(H2))
>>> plt.plot(Hc.Hcal.x,Hn)
>>> plt.xlabel('fGHz')
>>> plt.ylabel(u'$|H(f)|^{2}$')
>>> plt.title('Normalized channel')
```

we set the range of the SNR adjustment parameter $\rho$ to address the region around $SNR=0dB$ (including both low and high SNR regions)

```python
>>> rho = np.logspace(-10,-8,100)
>>> SNRrho = 10*np.log10(Pt*Hn[None,:]*rho[:,None]/(N0*B))
```

```python
>>> i1 = 0
>>> i2 = 180
>>> i3 = 1200
>>> plt.semilogx(rho,SNRrho[:,i1],label=str(fGHz[i1])+' GHz')
>>> plt.semilogx(rho,SNRrho[:,i2],label=str(fGHz[i2])+' GHz')
>>> plt.semilogx(rho,SNRrho[:,i3],label=str(fGHz[i3])+' GHz')
>>> plt.ylabel('SNR (dB)',fontsize=18)
>>> plt.xlabel(r'$\rho$',fontsize=18)
>>> plt.legend(loc='best')
```

```python
>>> #rh0 = 1/np.sum(H2)
... #print rh0
```

$$\int_B E_{no CSI}(f) df = \frac{Pt B }{B} = P_t (Watt)$$

# Water Filling

```python
>>> # No CSIPt/B
... Encsi  = Pt/(B+df*1e9)
>>> #tightEncsi  = Pt/B
... # With CSI (water filling)
... mu = np.logspace(-11,-8,80)
>>> # mu,rho,f
... IC    = N0/(rho[None,:,None]*Hn[None,None,:])
>>> Ecsi  = np.maximum(mu[:,None,None]-IC,0)
>>> #normalization of the water filling energy distribution
... Ecsi_f=df*1e9*np.sum(Ecsi,axis=2)
>>> Ecsin = Pt*Ecsi/Ecsi_f[:,:,None]
>>> Ecsin[np.where(np.isnan(Ecsin))]=0
>>> assert abs(df*1e9*np.sum(Ecsin,axis=2)-Pt).all()<1e-8,'bad normalization'
```

We are aiming to compare the channel capacity depending whether or not we are exploiting a knowledge of the channel at the transmitter.

In the absence of CSI the total power is uniformily spread over the bandwidth $B$.

$$E_{no CSI}=\frac{Pt}{B} \; (J)$$

$$\int_B E_{no CSI}(f) df = \frac{Pt B }{B} = P_t (Watt)$$

When the CSI is available at the transmitter the water filling algorithm is used :

Let define $$\tilde{|H(f)|^2}=\frac{|H(f)|^2}{\int_B |H(f)|^2df}$$

$$E_{CSI}=\left(\mu - \frac{N_0}{\rho \tilde{|H(f)|^2}} \right)^{+}$$

We have now to normalize the power spectral density $E_{CSI}(f)$ such as the condition $$\int_B E_{CSI}(f) df = P_t $$
for all values of $\rho$ and $\mu$

$$\tilde{E}_{CSI}(f)=P_t \frac{\left(\mu - \frac{N_0}{\rho \tilde{|H(f)|^2}} \right)^{+}}{\int_B \left(\mu - \frac{N_0}{\rho \tilde{|H(f)|^2}} \right)^{+}df} (J)$$

```python
>>> Eu = np.ones((100,1601))*Encsi
>>> assert abs(sum(Eu,axis=1)*df*1e9-Pt<1e-8).all(), 'bad normalization'
```

```python
>>> i1 = 50
>>> plt.plot(fGHz,Eu[i1,:]*1e12,label='Uniform')
>>> plt.plot(fGHz,Ecsin[0,i1,:]*1e12,label='Water Filling')
>>> plt.title('Comparison of Energy density (CSI vs No CSI)'+r'$\rho$ ='+str(rho[i1]))
>>> plt.legend(loc='best')
>>> xlabel('frequency (GHz)')
>>> ylabel('Energy (pJ)')
```

```python
>>> #cu = np.log(1+(Pt*rho[:,None]*Hn[None,:])/(N0*B))/np.log(2)
... cu = np.log(1+(Eu*rho[:,None]*Hn[None,:])/N0)/np.log(2)
>>> cw = np.log(1+(Ecsin*rho[None,:,None]*Hn[None,None,:])/N0)/np.log(2)
```

```python
>>> CU=np.sum(cu,axis=1)*df
>>> CW=np.sum(cw,axis=2)*df
```

```python
>>> eta=CW/CU[None,:]
```

```python
>>> eta.max()
1.2276704892638046
```

```python
>>> plt.plot(CU)
>>> plt.plot(CW[4,:])
[<matplotlib.lines.Line2D at 0x7f9b52355f90>]
```

```python
>>> SNRu = (Pt*rho[:,None]*Hn[None,:])/(N0*B)
>>> SNRudB= 10*np.log10(np.sum(SNRu,axis=1)*df)
>>> ext = [SNRudB[0],SNRudB[-1],mu[0],mu[-1]]
>>> plt.imshow(eta,origin='lower',extent=ext)
>>> plt.colorbar()
>>> plt.axis('tight')
>>> plt.title('Capacity ratio')
>>> plt.xlabel('SNR(dB)')
>>> plt.ylabel(r'$\mu$')
```

```python
>>> i1 = 20
>>> plt.semilogx(mu*1e12,eta[:,i1])
>>> plt.xlabel('$\mu (pJ)$')
```

```python
>>> CU=np.sum(cu,axis=1)*df
```

```python
>>> plt.ion()
>>> CU=np.sum(cu,axis=1)*df
>>> plt.semilogx(rho,CU)
>>> CW=np.sum(cw,axis=2)*df
>>> plt.semilogx(rho,CW[10,:],'r')
[<matplotlib.lines.Line2D at 0x7f9b520ca550>]
```

CW.shape$$C_{BLAST} = \int_{B} \log_2 \textrm{det}(\mathbf{I}_t + \frac{Pt }{N_0 B N_t} \mathbf{H}^{\dagger}(f)\mathbf{H}(f))df \;\; \textrm{(bit/s)}$$

```python
>>> plt.plot(SNRudB,CU,label='No CSI')
>>> plt.plot(SNRudB,CW[10,:],label='With CSI (Water Filling) ')
>>> plt.xlabel('SNR (uniform)')
>>> plt.ylabel('Gbit/s')
>>> plt.legend(loc='best')
```

Some words about the units. 
Notice that $nJ$ is homogeneous to $mW/MHz$. For example for intentional UWB system, the emission mask is limited to -41.3 dBm/MHz between 3.1GHz and 10.6GHz. This corresponds to an energy level of $E=10^{-4.13}\times 10^{-9}=0.0741 pJ$. It is important to keep this values in mind when doing water filling because in some situation the total power limit is not the only constraint. There might also have limitation about the maximum level of energy at a given frequency from emission mask from the regulatory body.

$$C_u = \int_{B} \log_2 (1 + \frac{\rho P_t }{N_0 B}|H(f)|^2) df $$

$$C_{BLAST} = \int_{B} \log_2 \textrm{det}(\mathbf{I}_t + \frac{Pt }{N_0 B N_t} \mathbf{H}^{\dagger}(f)\mathbf{H}(f))df \;\; \textrm{(bit/s)}$$

```python

```
