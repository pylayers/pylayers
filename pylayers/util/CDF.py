# -*- coding:Utf-8 -*-
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pdb
#import mplrc.ieee.transaction
#
# mplrc is a python module which provides an easy way to change
# matplotlib's plotting configuration for specific publications.
# git clone https://github.com/arsenovic/mplrc.git
#
#
from matplotlib import rcParams


rcParams['text.usetex'] = True
rcParams['text.latex.unicode'] = True


class CDF(object):
    def __init__(self, ld, filename=''):
        """
        cdf = CDF(ld)
        
        Parameters
        ----------
        
        ld : list 
            list of dictionnary
        filename : string


        Notes
        -----

        d0 = ld[0]

        d0['bound']       : abscisse bounds of the cdf
        d0['values']      : valeurs
        d0['xlabel']      :
        d0['ylabel']      :
        d0['legend']      : legend
        d0['title]        : title
        d0['filename]     : filename
        d0['linewidth']   : linewidth

        """
        
        self.ld = ld
        self.filename = filename
        if self.filename == '':
            self.save=False
        else:
            self.save=True
            plt.rcParams['xtick.labelsize'] ='x-large'
            plt.rcParams['ytick.labelsize'] ='x-large'
            plt.rcParams['axes.labelsize']  ='large'
            plt.rcParams['font.weight']     ='normal'
            plt.rcParams['xtick.minor.size']=2
            plt.rcParams['legend.fontsize'] = 'xx-large'
            plt.rcParams['font.size']       =20
            plt.rcParams['grid.linewidth']  =3.5
            plt.rcParams['xtick.major.pad'] =20

       
        self.cdf = []
        for d in self.ld:
            if d.has_key('bound'):
                bound = d['bound']
            else:
                bound = np.arange(d['values'].min(),d['values'].max(),len(d['values']*0.1))
            values = d['values']
            Nv = len(values)
            cdf = np.array([])
            for k in bound:
                u = np.nonzero(values <= k)
                lu = len(u[0]) / (Nv * 1.0)
                cdf = np.hstack((cdf, lu))

            self.cdf.append(cdf)

    def show(self,**kwargs):
        """ show cdf
        """
        if 'fig' not in kwargs:
            f = plt.figure(**kwargs)
        else:
            f = kwargs['fig']

        if 'ax' not in kwargs:
            ax = f.add_subplot(111)
        else:
            ax == kwargs['ax']

        leg = []
        c = []

        for k in range(len(self.ld)):

            d = self.ld[k]
            if d.has_key('bound'):
                bound = d['bound']
            else:
                bound = np.arange(d['values'].min(),d['values'].max(),len(d['values']*0.1))
            if d.has_key('marker'):
                marker = d['marker']
            else:
                marker = ''
            if d.has_key('markersize'):
                markersize = d['markersize']
            else:
                markersize = 5
            if d.has_key('markercolor'):
                markercolor = d['markercolor']
            else:
                markercolor = 'k'
            if d.has_key('markerfrequency'):
                markerfrequency = d['markerfrequency']
            else:
                markerfrequency = 10
            if d.has_key('linewidth'):
                linewidth = d['linewidth']
            else:
                linewidth = 1
            if d.has_key('linestyle'):
                linestyle = d['linestyle']
            else:
                linestyle = '-'
            if d.has_key('color'):
                color = d['color']
            else:
                color ='k'
            if d.has_key('legend'):
                legend = d['legend']
            else:
                legend=''
            if k == 0:
                if d.has_key('x_label'):
                    xlabel=d['x_label']
                else:
                    xlabel=''
                if d.has_key('y_label'):
                    ylabel=d['y_label']
                else:
                    ylabel=''

#                       leg.append(legend)
            cdf = self.cdf[k]
            c.append(ax.plot(bound, cdf, marker=marker,
              markevery=markerfrequency, ms=markersize, mfc=markercolor,
                             ls=linestyle, c=color, linewidth=linewidth,
                             label=legend))
            plt.xlabel(xlabel)

        plt.ylabel(ylabel)
        ax.legend(loc='best', scatterpoints=1, numpoints=1.)
        plt.grid()
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        if self.save :
            plt.savefig(self.filename + '.pdf', format='pdf',
                        bbox_inches='tight', pad_inches=0)
            plt.savefig(self.filename + '.eps', format='eps',
                        bbox_inches='tight', pad_inches=0)


if __name__ == "__main__":
    d0 = {}
    d0['values'] = sp.randn(1000)
    d0['bound'] = np.arange(-10, 10, 0.1)
    d0['xlabel'] = 'xlabel'
    d0['ylabel'] = 'ylabel'
    d0['legend'] = 'legend '
    d0['markersize'] = 3
    d0['markercolor'] = 'red'
    d0['markerfrequency'] = 2
    d0['title'] = 'title'
    d0['color'] = 'black'
    d0['marker'] = 'o'
    d0['linestyle'] = '-'
    d0['linewidth'] = 3
    d0['filename'] = 'essai.png'
    d1 = {}
    d1['values'] = 4 * sp.randn(1000)
    d1['bound'] = np.arange(-10, 10, 0.1)
    d1['xlabel'] = 'xlabel'
    d1['ylabel'] = 'ylabel'
    d1['legend'] = 'legend '
    d1['markersize'] = 3
    d1['markercolor'] = 'blue'
    d1['linestyle'] = '-'
    d1['color'] = 'black'
    d1['markerfrequency'] = 2
    d1['title'] = 'title'
    d1['marker'] = 'o'
    d1['linewidth'] = 3
    lv = [d0, d1]
    c = CDF(lv, 'fig')
    c.show()
