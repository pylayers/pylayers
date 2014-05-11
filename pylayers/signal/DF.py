#!/usr/bin/python
# -*- coding: latin1 -*-
"""

DF Class
========

This class implements a LTI digital filter

.. autosummary::
    :toctree: generated/

    DF.__init__
    DF.filter
    DF.order
    DF.freqz
    DF.ellip_bp
    DF.remez
    DF.zplane
    DF.ir
"""
import numpy as np
from scipy import io
import scipy.signal as si
import matplotlib.pyplot as plt
#from EnergyDetector import *

class DF(object):
    """ Digital Filter Class

    Methods
    -------

    filter : filter a signal
    freqz  : Display transfer function
    ellip_bp : elliptic bandpath filter
    remez : FIR design

    """
    def __init__(self,b=np.array([1]),a=np.array([1])):
        self.b = b
        self.a = a

    @property
    def z(self):
        return(self._z)

    @property
    def p(self):
        return(self._p)

    @property
    def a(self):
        return(self._a)

    @property
    def b(self):
        return(self._b)

    @a.setter
    def a(self,a):
        self._a = a
        self._p = np.poly1d(a).r
        if len(a)==1:
            self.fir=True
        else:
            self.fir=False

    @b.setter
    def b(self,b):
        self._b = b
        self._z = np.poly1d(b).r

    @p.setter
    def p(self,p):
        self._p = p
        self._a = np.poly1d(p,r=True).c

    @z.setter
    def z(self,z):
        self._z = z
        self._b = np.poly1d(z,r=True).c

    def __mul__(self,df):
        """
        """

        b1 = self.b
        a1 = self.a
        b2 = df.b
        a2 = df.a

        pb1 = np.poly1d(b1)
        pa1 = np.poly1d(a1)
        pb2 = np.poly1d(b2)
        pa2 = np.poly1d(a2)

        rpb1 = pb1.r
        rpb2 = pb2.r
        rpa1 = pa1.r
        rpa2 = pa2.r

        F   = DF()
        F.p = np.hstack((rpa1,rpa2))
        F.z = np.hstack((rpb1,rpb2))

        return(F)


    def match(self):
        """ return a match filter
        """

        if self.fir:
            MF = DF()
            MF.b = np.conj(self.b[::-1])
            MF.a = self.a
            return(MF)

    def filter(self,x):
        """ filter x

        Parameters
        ----------

        x : np.array
            input signal

        Returns
        -------

        y : np.array
            output signal

        """

        y  = si.lfilter(self.b,self.a,x)
        return(y)

    def order(self):
        """ Returns filter order
        """
        return self.order

    def freqz(self):
        """ freqz : display filter transfer function

        """
        (w,h)    = si.freqz(self.b,self.a)
        self.w   = w
        self.h   = h
        plt.subplot(211)
        plt.plot(w/np.pi,20*np.log10(abs(h)+1e-15))
        plt.ylabel('dB')
        plt.title('Modulus')
        plt.grid()
        plt.subplot(212)
        plt.plot(w/np.pi,np.angle(h)*180./np.pi)
        plt.ylabel('deg')
        plt.xlabel('Relative frequency')
        plt.title('Phase')
        plt.grid()
        #show()

    def ellip_bp(self,wp,ws,gpass=0.5,gstop=20):
        """ Elliptic Bandpath filter

        Parameters
        ----------

        wp :
        ws :
        gpass
        gstop :

        See Also
        --------

        iirdesign

        """
        (b,a)   = si.iirdesign(wp,ws,gpass,gstop,analog=0,ftype='ellip',output='ba')
        self.b  = b
        self.a  = a

    def remez(self,numtaps=401,bands=(0,0.1,0.11,0.2,0.21,0.5),desired=(0.0001,1,0.0001)):
        """ FIR design Remez algorithm

        Parameters
        ----------

        numtaps : int
        bands  : tuple of N floats
        desired : tuple of N/2 floats

        Examples
        --------

        flt.remez(numtaps=401,bands=(0,0.1,0.11,0.2,0.21,0.5),desired=(0.0001,1,0.0001)):

        numtaps = 401
        bands   = (0,0.1,0.11,0.2,0.21,0.5)
        desired = (0.0001,1,0.0001))

        flt.remez( numtaps , bands , desired)

        """

        self.a = array([1])
        self.b = remez(numtaps, bands, desired, weight=None, Hz=1, type='bandpass', maxiter=25, grid_density=16)

    def zplane(self):
        """ Display filter in the complex plane

        Parameters
        ----------

        """
        A  = np.poly1d(self.a)
        B  = np.poly1d(self.b)
        ra = A.r
        rb = B.r

        t = np.arange(0,2*np.pi+0.1,0.1)
        plt.plot(np.cos(t),np.sin(t),'k')

        plt.plot(np.real(ra),np.imag(ra),'xr')
        plt.plot(np.real(rb),np.imag(rb),'ob')
        plt.axis('equal')
        plt.show()

    def ir(self,N,show=False):
        """ show impulse response
        """
        ip    = np.zeros(N)
        ip[0] = 1
        ir    = self.filter(ip)
        if show:
            plt.stem(np.arange(N),rip)
            plt.show()
        return(ir)


if __name__ == "__main__":
    fe  = 10
    fN  = fe/2.0
    wt  = 0.01
    ws = np.array([3.168,3.696])/fN
    wp = [ws[0]+wt,ws[1]-wt]
    flt = DF()
    gpass = 0.5
    gstop = 40
    flt.ellip_bp(wp,ws,gpass,gstop)
    flt.zplane()
    flt.freqz()
    x = np.random.randn(1000)
    y = flt.filter(x)
    fo = DF()
    y2 = fo.filter(x)

