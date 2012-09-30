#!usr/bin/python
# -*- coding: latin1 -*-
import os
import Image
import numpy as np
from pylayers.util import geomutil as geu
from pylayers.util import pyutil as pyu
import matplotlib.pyplot as plt
from pylayers.util.easygui import *


class SelectL(object):
    """ Associates a Layout and a figure

    o : toggle overlay
    p : create point
    x : save graph in .str2 file
    w : show all layers


    """
    def __init__(self, L, fig):
        """ SelectL is a class which associate a Layout and a figure

        Parameters
        ----------
        L   : Layout
        fig : figure

        """
        self.L = L
        self.fig = fig
        self.ax = self.fig.add_subplot(111)
        self.ax.axis(self.L.display['box'])
        plt.title('Init')
        self.text = self.ax.text(0.05, 0.95, 'Selected : none',
                                 transform=self.ax.transAxes, va='top')
        self.pt = []
        self.seg = []
        self.coseg = []
        self.pt1 = np.array([])
        self.pt2 = np.array([])
        self.selected_pt1 = -1
        self.selected_pt2 = -1
        self.selected_edge = -1
        self.current_layer = self.L.display['activelayer']
        self.npsel = 0
        self.nedge_sel = 0
        self.indp = 0
        self.state = 'Init'
        self.nsel = 0 
        self.update_state()

    def show(self, clear=False, dnodes=True, dedges=True,  font_size=10, title='Init'):
        """ show

        Parameters
        ----------
        clear     : boolean
        dnodes    : boolean
        dedges    : boolean
        dlabels   : boolean
        font_size : integer
        title     : string

        """
        #laxe  = self.ax.get_axes()
        #xmin,xmax=laxe.xaxis.get_view_interval()
        #ymin,ymax=laxe.yaxis.get_view_interval()
        #xmin,xmax,ymin,ymax = self.g.ax
        #ax = plt.axis('tight')
        axis = self.ax.axis()
        print('show : axis',axis) 
        #plt.axis(ax)
        self.L.display['clear'] = clear
        self.L.display['nodes'] = dnodes
        self.L.display['edges'] = dedges
        self.L.display['fontsize'] = font_size
        self.L.display['title'] = title
        self.L.display['ednodes'] = True
        self.ax = self.L.showGs(self.ax,axis=axis)

    def OnPress(self, event):
        """
        OnPress(event)
        """
        self.nsel = 0
        self.ptsel = np.array([])
        self.evt = event.key
        print "Evenement :", self.evt
        self.new_state()

    def OnClick(self, event):
        """
        OnClick(event)
        """
        self.nsel = 0
        self.ptsel = np.array([])
        xmin, xmax, ymin, ymax = self.ax.axis()
        print xmin,xmax,ymin,ymax
        dx = xmax - xmin
        dy = ymax - ymin
        dd = np.minimum(dx, dy)
        if event.button == 1 and event.inaxes:
            self.evt = 'lclic'
            x = event.xdata
            y = event.ydata
            self.ptsel = np.array((x, y))
            self.nsel = self.L.ispoint(self.ptsel, dd / 100)

        if event.button == 2 and event.inaxes:
            self.evt = 'cclic'
            x = event.xdata
            y = event.ydata
            self.ptsel = np.array((x, y))
            self.nsel = self.L.ispoint(self.ptsel, dd / 100)

        if event.button == 3 and event.inaxes:
            self.evt = 'rclic'
            x = event.xdata
            y = event.ydata
            self.ptsel = np.array((x, y))
            self.nsel = self.L.ispoint(self.ptsel, dd / 100)

        print "Selected point coord : ", self.ptsel
        print "Selected point number: ", self.nsel
        if self.nsel > 0:
            print "Selected Edge : ", self.nsel

        self.new_state()

    def update_state(self):
        """ update state
        """

        if self.state == 'Init':
            self.show(clear=True)
            self.selected_edge = 0
            self.selected_pt1 = 0
            self.selected_pt2 = 0
            try:
                del self.pt_previous
            except:
                pass
            self.ax.title.set_text('Init : '
                                   +self.L.display['activelayer'])
            try:
                self.p1[0].set_visible(False)
            except:
                pass
            try:
                self.p2[0].set_visible(False)
            except:
                pass
            try:
                self.segment[0].set_visible(False)
            except:
                pass
            #
            # If Layout has no point go to CP state
            #
            if self.L.Nn==0:
                self.state='CP'
                self.update_state()

        if self.state == 'SP1':
            self.show(clear=False)
            self.ax.title.set_text('Selected node : %d ' % (self.nsel))
            self.selected_pt1 = self.nsel
            self.pt1 = np.array(self.L.Gs.pos[self.nsel]).reshape(2, 1)
            self.pt_previous = self.pt1
            self.p1 = self.ax.plot([self.pt1[0]], [self.pt1[1]], 'o', visible=True)
            self.p1[0].set_color('yellow')
            self.p1[0].set_ms(10)
            self.p1[0].set_alpha(0.4)
            try:
                self.p2.set_visible(False)
            except:
                pass

        if self.state == 'SP2':
            self.p1[0].set_color('green')
            self.ax.title.set_text('Selected node : %d ' % (self.nsel))
            self.selected_pt2 = self.nsel
            self.pt2 = np.array(self.L.Gs.pos[self.nsel]).reshape(2, 1)
            self.pt_previous = self.pt2
            self.p2 = self.ax.plot([self.pt2[0]], [self.pt2[1]], 'o', visible=True)
            self.p2[0].set_color('green')
            self.p2[0].set_ms(10)
            self.p2[0].set_alpha(0.4)
            self.ax.title.set_text('SP2')

        if self.state == 'SS':
            try:
                self.p1[0].set_visible(False)
            except:
                pass
            try:
                self.p2[0].set_visible(False)
            except:
                pass
            self.selected_edge = self.nsel
            nse = self.nsel
            ta, he = self.L.Gs.neighbors(nse)
            pta = np.array(self.L.Gs.pos[ta])
            phe = np.array(self.L.Gs.pos[he])
            alpha = self.L.display['alpha']
            self.current_layer = self.L.Gs.node[nse]['name']
            self.L.display['activelayer'] = self.current_layer
            #self.seg       = linet(self.ax,pta,phe,alpha,'red',3.5)
            segdico = self.L.Gs.node[nse]
            self.show(clear=False)
            self.segment = self.ax.plot([pta[0],phe[0]],
                                        [pta[1],phe[1]],
                                        'r',linewidth=3, visible=True)
            if 'ss_name' in segdico:
                cosegname = segdico['ss_name']
                titre = 'Select Segment : %d (%d->%d) Layer : %s Coseg : %s ' % (nse, ta, he, self.current_layer, cosegname)
            else:
                titre = 'Select Segment : %d (%d->%d) Layer : %s' % (nse, ta, he, self.L.Gs.node[nse]['name'])

            self.ax.title.set_text(titre)
            self.L.show_nodes(ndlist=[nse], size=200, color='r', alpha=0.5)

        if self.state == 'SSS':
            nse = self.selected_edge
            segdico = self.L.Gs.node[nse]
            zmin    = segdico['ss_zmin']
            zmax    = segdico['ss_zmax']
            self.ax.title.set_text('SSS : '+self.L.Gs.node[nse]['name']+' ['+str(zmin)+','+str(zmax)+']')
            self.segment[0].set_color('blue')
        #
        # Create Point state
        #
        if self.state == 'CP':
            try:
                self.segment[0].set_visible(False)
            except:
                pass
            try:
                self.segment1[0].set_visible(False)
            except:
                pass
            try:
                self.segment2[0].set_visible(False)
            except:
                pass
            self.show(clear=False,title='CP lclic : free point, rclic same x, cclic same y')

        #
        # Create Point on Segment state
        #

        if self.state == 'CPS':
            self.selected_edge1 = self.nsel
            ta, he = self.L.Gs.neighbors(self.nsel)
            self.pta1 = np.array(self.L.Gs.pos[ta])
            self.phe1 = np.array(self.L.Gs.pos[he])
            self.current_layer = self.L.Gs.node[self.nsel]['name']
            self.L.display['activelayer'] = self.current_layer
            self.segment1 = self.ax.plot([self.pta1[0],self.phe1[0]],
                                        [self.pta1[1],self.phe1[1]],
                                        'g',linewidth=3, visible=True)
            try:
                self.segment2[0].set_visible(False)
            except:
                pass

        if self.state == 'CPSS':
            self.selected_edge2 = self.nsel
            ta, he = self.L.Gs.neighbors(self.nsel)
            self.pta2 = np.array(self.L.Gs.pos[ta])
            self.phe2 = np.array(self.L.Gs.pos[he])
            self.current_layer = self.L.Gs.node[self.nsel]['name']
            self.L.display['activelayer'] = self.current_layer
            self.segment2 = self.ax.plot([self.pta2[0],self.phe2[0]],
                                        [self.pta2[1],self.phe2[1]],
                                        'c',linewidth=3, visible=True)



        print self.state
        print self.nsel
        print self.selected_pt1
        print self.selected_pt2
        self.fig.canvas.draw()


    def new_state(self):
        """ layout editor state machine

        Parameters
        ----------

        'l'  : select activelayer
        'i'  : back to init state
        'e'  : edit segment
        'h'  : add subsegment
        'd'  : delete subsegment
        'r'  : refresh
        'o'  : toggle overlay
        'm'  : toggle mode (point or segment)  
        'z'  : change display parameters
        'q'  : quit interactive mode
        'x'  : save .str2 file
        'w'  : display all layers

        """
        sl = self.L.sl
        cold = pyu.coldict()
        print "In State ",self.state
        print "In Event ",self.evt
        #
        # Choose layers to visualized
        #
        if self.evt == "l":
            listchoices = self.L.name.keys()
            self.L.display['layers'] = multchoicebox('message',
                                                     'titre', listchoices)
            self.state = 'Init'
            self.update_state()
            return
        #
        # Increment layer
        #
        if self.evt=='=':
            N = len(self.L.display['layerset'])
            index = self.L.display['layerset'].index(self.L.display['activelayer'])
            self.L.display['activelayer'] = self.L.display['layerset'][(index+1) % N]
            self.current_layer = self.L.display['activelayer']
            self.update_state()
            return

        if self.evt=='$':
            N = len(self.L.display['layerset'])
            index = self.L.display['layerset'].index(self.L.display['activelayer'])
            self.L.display['activelayer'] = self.L.display['layerset'][(index-1) % N]
            self.current_layer = self.L.display['activelayer']
            self.update_state()
            return
        #
        # e : Edit
        #
        if self.evt == 'i':
            self.state = 'Init'
            self.update_state()
            return

        if self.evt == 'e':
            if (self.state == 'SS') | (self.state =='SSS'):
                self.L.edit_edge(self.selected_edge)
                self.state = 'Init'
                self.update_state()
                return
            if self.state == 'SP1':
                print "Write edit_node"
        #
        # h : add subsegment
        #
        if self.evt == 'h':
            if self.state == 'SS':
                self.L.add_subseg(self.selected_edge,self.current_layer)
                self.state = 'SSS'
                self.update_state()
                return
        #
        # d : delete 
        #
        if self.evt == 'd':
            if  self.state == 'SP1':
                self.state = 'Init'
                self.L.del_node(self.selected_pt1)
                self.update_state()
                return
            if self.state == 'SS':
                self.L.del_edge(self.selected_edge)
                self.state = 'Init'
                self.update_state()
                return
            if self.state == 'SSS':
                self.L.del_subseg(self.selected_edge)
                self.state = 'Init'
                self.update_state()
                return

        if self.evt == 'c':
            if self.state == 'Init':
                #pt = plt.ginput(4)
                #self.ax = plt.gca()
                x1 = self.ax.get_xbound()
                y1 = self.ax.get_ybound()
                #self.ax.autoscale(True)
                #axx = plt.axis()
                #ndlist,edlist = self.L.get_zone([x1[0],x1[1],y1[0],y1[1])
                print x1,y1
                #print axx
                #print ndlist
                #self.L.del_node(ndlist)
                #self.show(clear=True, title='Init')
                return

        #
        # r : Refresh
        #
        if self.evt == 'r':
            plt.axis('tight')
            self.show(clear=True)
            self.state = 'Init'
            self.update_state()
            return
        #
        # o : Toggle overlay
        #
        if self.evt == 'o':
            if self.L.display['overlay']:
                self.L.display['overlay'] = False
                self.update_state()
            else:
                self.L.display['overlay'] = True
                self.update_state()
            return

        #
        # m : Toggle mode edition Point | Segment 
        #
        if self.evt == 'm':
            if self.state == "Init":
                self.state = "CP"
            elif self.state == "CP":
                self.state = "Init"
            self.update_state()
            return

        if self.evt == 'z':
        # change display parameters
            self.L.displaygui()
            self.show(clear=True)
            return

        if self.evt == 'q':
        # quit interactive mode
            self.fig.canvas.mpl_disconnect(self.L.cid1)
            self.fig.canvas.mpl_disconnect(self.L.cid2)
            return
        if self.evt == 'x':
        # save structure
            racine, ext = os.path.splitext(self.L.filename)
            filename = racine + '.str2'
            self.L.savestr2(filename)
            print "structure saved in ", filename
            return

        if self.evt == 'n':
            self.L.display['ndlabel'] = not self.L.display['ndlabel']  
            self.L.display['edlabel'] = not self.L.display['edlabel']  
            self.show(clear=True, title=self.L.display['activelayer'])
            self.fig.canvas.draw()
            return

        if self.evt == 'w':
        # display all layer
            self.L.display['activelayer'] = self.L.name.keys()
            self.show(clear=True, title=self.L.display['activelayer'])
            return
        #
        # Left clic and selected node is a Point 
        #
        if (self.evt == 'lclic') & (self.nsel < 0):

            if self.state=='Init':
                self.state = 'SP1'
                self.update_state()
                return

            if self.state=='SP1':
                if self.nsel != self.selected_pt1:
                    self.state = 'SP2'
                    self.update_state()
                    return

                   

        #
        # Left clic and selected node is a segment 
        #

        if (self.evt == 'lclic') & (self.nsel > 0):
            if self.state=='Init':
                self.state = 'SS'
                self.update_state()
                return

            if self.state=='SS':
                self.nsel = self.selected_edge
                segdico = self.L.Gs.node[self.nsel]
                if 'ss_name' in segdico:
                    self.state = 'SSS'
                self.update_state()
                return
            
        #
        # Right clic and selected node is a point
        #

        if (self.evt == 'rclic') & (self.nsel < 0):
            if self.state=='SP1':
                if self.nsel==self.selected_pt1:
                    self.state = 'Init'
                    self.update_state()
                    return
        #
        # Right clic and selected node is a segment
        #

        if (self.evt == 'rclic') & (self.nsel > 0):
            if self.state=='SS':
                self.state = 'Init'
                self.update_state()
                return

            if self.state=='SSS':
                self.state = 'SS'
                self.update_state()
                return

            if self.state == 'CP':
            # create point on edge
                self.state = 'CPS'
                self.update_state()
                return
            
            if (self.state == 'CPS') & (self.nsel!= self.selected_edge1):
            # create point on edge
                self.state = 'CPSS'
                self.update_state()
                return
        #
        # Left clic
        #
        if (self.evt == 'lclic'):
            # add free node
            if self.state == 'CP':
                self.L.add_fnod(self.ptsel)
                self.pt_previous = self.ptsel
                self.update_state()
                return
           
            if self.state == 'SP2':
                ta = self.selected_pt1
                he = self.selected_pt2
                self.nsel  = self.L.add_edge(ta, he,name=self.current_layer)
                self.state = 'Init'
                self.update_state()
                return
            
            # create point on segment 
            if self.state == 'CPS':
                pt_new = geu.ptonseg(self.pta1, self.phe1, self.ptsel)
                pd1 = pt_new - self.pta1
                pd2 = self.phe1 - self.pta1
                alpha = np.sqrt(np.dot(pd1, pd1)) / np.sqrt(np.dot(pd2, pd2))
                if (pt_new != []):
                    # calculate alpha
                    self.L.add_none(self.selected_edge1, 1. - alpha)
                    self.current_layer = self.L.Gs.node[self.selected_edge1]['name']
                    self.state = 'Init'
                self.update_state()
                return

            
        #
        # Right Clic event  
        #
        if (self.evt == 'rclic'):
            if self.state == 'CP':
                try:
                    self.ptsel[0] = self.pt_previous[0]
                    self.L.add_fnod(self.ptsel)
                    self.pt_previous = self.ptsel
                    self.update_state()
                    return
                except:
                    pass
            
            if self.state=='SP2':
                if self.nsel == self.selected_pt1:
                    self.p1[0].set_visible(False)
                    self.p2[0].set_visible(False)
                    self.nsel = self.selected_pt2
                    self.state = 'SP1'
                    self.update_state()
                    return
                if self.nsel == self.selected_pt2:
                    self.p1[0].set_visible(False)
                    self.p2[0].set_visible(False)
                    self.nsel = self.selected_pt1
                    self.state = 'SP1'
                    self.update_state()
                    return
             
            if self.state == 'CPS':
                self.state = 'CP'
                self.update_state()
                return
            
            if self.state == 'CPSS':
                self.state = 'CPS'
                self.update_state()
                return
        #
        # Center Clic event 
        #
        if (self.evt == 'cclic'):
            if self.state == 'CP':
                try:
                    self.ptsel[1] = self.pt_previous[1]
                    self.L.add_fnod(self.ptsel)
                    self.pt_previous = self.ptsel
                    self.update_state()
                except:
                    pass
        #
        # Left clic and selected node is a point 
        #


