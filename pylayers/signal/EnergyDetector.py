from numpy import *
from scipy.signal import *
from pylab import *
from pylayers.signal.DF import *
from pylayers.util.project import *
import matplotlib.pyplot as plt


class EnergyDetector(PyLayers):
    """
    EnergyDetector Class

    This class implements an Energy Detector object


    """
    def __init__(self,filter1,beta,filter2,pfa):
        self.filter1  = filter1
        self.beta     = beta
        self.filter2  = filter2
        self.pfa      = pfa

    def apply(self,x):
        self.xf1 = self.filter1.filter(x)
        self.xq  = ((self.beta)**2)*(self.xf1*self.xf1)
        y        = self.filter2.filter(self.xq)
        return(y)

if __name__ == "__main__":
    """
    """
    #
    # Signal entree
    #

    Npoints    = 4000

    Var_w      = 1
    aff        = {'aff':0,'dsp':0,'densite':0,'seuil':0}

    t  = linspace(0,1,Npoints)
    wn = sqrt(Var_w)*randn(Npoints)

    #
    # Filter 1 : Input Band Pass Filter
    #
    w1       = 0.05
    w2       = 0.1
    wt       = 0.01
    wp       = [w1,w2]
    ws       = [w1-wt,w2+wt]
    gpass    = 0.5
    gstop    = 30
    type     = 'ellip'

    filter1  = DF()
    filter1.ellip_bp(wp,ws,gpass,gstop)

    #
    # Filter 2 : averaging filter
    #
    N       = 40
    b       = (1./N)*ones(N)
    a       = array([1])
    filter2 = DF(b,a)

    #
    #
    #
    beta     = 1
    pfa      = 0.001

    ED = EnergyDetector(filter1,beta,filter2,pfa)
    y  = ED.apply(wn)




