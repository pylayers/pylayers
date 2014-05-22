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
from numpy import *
from scipy import io
from scipy.signal import *
from pylab import *
#from EnergyDetector import *

class DF(object):
    """ Digital Filter Class

    Methods
    -------

    filter : flilter a signal

    """
    def __init__(self,b=array([1]),a=array([1,-0.5])):
        self.b = b
        self.a = a

    def filter(self,x):
        y  = lfilter(self.b,self.a,x)
        return(y)

    def order(self):
        self.order = max(len(self.a),len(self.b))-1
        return self.order

    def freqz(self):
        """
        freqz : display filter transfer function
        """
        (w,h)    = freqz(self.b,self.a, whole = True)
        self.w   = w
        self.h   = h
        subplot(211)
        plot(w/(2*pi),20*log10(abs(h)+1e-15))
        ylabel('dB')
        title('Modulus')
        grid()
        subplot(212)
        plot(w/(2*pi),angle(h)*180./pi)
        ylabel('deg')
        xlabel('Relative frequency')
        title('Phase')
        grid()
        #show()

    def ellip_bp(self,wp,ws,gpass=0.5,gstop=20):
        """ Elliptic Bandpath filter

        wp :
        ws :
        gpass
        gstop :

        See Also
        --------

        iirdesign

        """
        (b,a)   = iirdesign(wp,ws,gpass,gstop,analog=0,ftype='ellip',output='ba')
        self.b  = b
        self.a  = a

    def remez(self,numtaps=401,bands=(0,0.1,0.11,0.2,0.21,0.5),desired=(0.0001,1,0.0001)):
        """
        FIR design Remez algorithm

         flt.remez(numtaps=401,bands=(0,0.1,0.11,0.2,0.21,0.5),desired=(0.0001,1,0.0001)):

         numtaps = 401
        bands   = (0,0.1,0.11,0.2,0.21,0.5)
        desired = (0.0001,1,0.0001))

        flt.remez( numtaps , bands , desired)

        """
        self.a = array(1)
        self.b = remez(numtaps, bands, desired, weight=None, Hz=1, type='bandpass', maxiter=25, grid_density=16)

    def zplane(self):
        """
        Display filter in the complex plane
        """
        A  = poly1d(self.a)
        B  = poly1d(self.b)
        ra = A.r
        rb = B.r

        t = arange(0,2*pi+0.1,0.1)
        plot(cos(t),sin(t),'k')

        plot(real(ra),imag(ra),'xr')
        plot(real(rb),imag(rb),'ob')
        axis('equal')
        show()

    def ir(self,N):
        ip    = zeros(N)
        ip[0] = 1
        rip   = self.filter(ip)
        stem(arange(N),rip)
        show() 


if __name__ == "__main__":
    fe  = 10
    fN  = fe/2.0
    wt  = 0.01
    ws = array([3.168,3.696])/fN
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

