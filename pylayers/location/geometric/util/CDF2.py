# -*- coding:Utf-8 -*-
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pdb
import mplrc.ieee.transaction
from matplotlib import rcParams
rcParams['text.usetex'] = True
rcParams['text.latex.unicode'] = True


class CDF(object):

    def __init__(self, ld, filename):
        """
        cdf = CDF(ld)

        ld is a list of dictionnary

        d0 = ld[0]

        d0['bound']       : bornes en abscisses de la cdf 0
        d0['values']      : valeurs
        d0['xlabel']      :
        d0['ylabel']      :
        d0['legend']      : legend
        d0['title]        : title
        d0['filename]     : filename
        d0['linewidth']   : linewidth

        """
        self.ld = ld
        self.parmsh = {}
        self.parmsh['file'] = True
        self.filename = filename
        self.cdf = []

        for d in self.ld:
            bound = d['bound']
            values = d['values']
            Nv = len(values)
            cdf = np.array([])
            for k in bound:
                u = np.nonzero(values <= k)
                lu = len(u[0]) / (Nv * 1.0)
                cdf = np.hstack((cdf, lu))

            self.cdf.append(cdf)

    def show(self):
        """
        show()
        """
        f = plt.figure()
        leg = []
        c = []
        ax = f.add_subplot(111)
        for k in range(len(self.ld)):
            d = self.ld[k]
            bound = d['bound']
            marker = d['marker']
            markersize = d['markersize']
            markercolor = d['markercolor']
            markerfrequency = d['markerfrequency']
            linewidth = d['linewidth']
            line = d['line']
            color = d['color']
            legend = d['legend']
            cdf = self.cdf[k]
            c.append(
                ax.plot(bound, cdf, marker=marker, markevery=markerfrequency,
                        ms=markersize, mfc=markercolor, ls=line, c=color, linewidth=linewidth, label=legend))
        plt.xlabel(self.ld[0]['xlabel'])
        plt.ylabel(self.ld[0]['ylabel'])
        ax.legend(loc='best', scatterpoints=1, numpoints=1.)
        plt.grid()
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        plt.savefig('./cdf/' + self.filename + '/' + self.filename + '.pdf',
                    format='pdf', bbox_inches='tight', pad_inches=0)
        plt.savefig('./cdf/' + self.filename + '/' + self.filename + '.eps',
                    format='eps', bbox_inches='tight', pad_inches=0)


if __name__ == "__main__":
    d0 = {}
    d0['values'] = sp.randn(1000)
    d0['bound'] = np.arange(-10, 10, 0.1)
    d0['xlabel'] = 'xlabel'
    d0['ylabel'] = 'ylabel'
    d0['legend'] = 'legend '
    d0['title'] = 'title'
    d0['marker'] = 'r-'
    d0['linewidth'] = 3
    d0['filename'] = 'essai.png'
    d1 = {}
    d1['values'] = 4 * sp.randn(1000)
    d1['bound'] = np.arange(-10, 10, 0.1)
    d1['xlabel'] = 'xlabel'
    d1['ylabel'] = 'ylabel'
    d1['legend'] = 'legend '
    d1['title'] = 'title'
    d1['marker'] = 'bo'
    d1['linewidth'] = 3
    lv = [d0, d1]
    c = CDF(lv)
