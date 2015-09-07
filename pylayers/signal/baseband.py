#!/usr/bin/python
#-*- coding:Utf-8 -*-
r"""
.. currentmodule:: pylayers.signal.baseband
"""
from numpy import *
import numpy as np
import pylayers.signal.bsignal as bs


class Baseband(bs.TUsignal):
        """
        """
        def __init__(self,seqbin,Tns=1000):
            """
            seq :
            Tns : symbol duration (ns)
            """
            N = len(seq)
            x = np.linspace(0,(N-1)*Tns,N)
            bs.TUsignal.__init__(self,x,seqbin)

class OFDM(Baseband):
    """
    """
    def __init__(self,qamvec,**kwargs):
        """
        """
        defaults={'NFFT':128,
                  'NCP':12
                 }
        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        NFFT = kwargs['NFFT']
        NCP = kwargs['NCP']

        y_vec = N*fft.ifft(qamvec)  # factor N compensates for the internal scaling of function ifft()
        #  prefix cyclic extraction
        cyclic_prefix = y_vec[N-L:]

        # ajout du prfixe cyclique
        # pour former le symbole
        # complet
        self.y = np.concatenate( (cyclic_prefix, y_vec))

    def channel(self,h):
        """ convolution with the channel impulse response
        """
        self.r = convolve(self.y,h.y)


    def demodulate(self,sync):
        """
        """
        pass
