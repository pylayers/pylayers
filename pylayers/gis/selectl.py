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

    Methods
    -------

    show(clear=F,dnodes=T,dedges=T,dlabels=T,font_size=10):
    OnPress(event)
    OnClick(event)
    call_editor()
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
        self.text = self.ax.text(0.05, 0.95, 'selected : none',
                                 transform=self.ax.transAxes, va='top')
        self.pt = []
        self.seg = []
        self.coseg = []
        self.pt1 = np.array([])
        self.pt2 = np.array([])
        self.selected_pt1 = -1
        self.selected_pt2 = -1
        self.selected_edge = -1
        self.current_layer = "WALL"
        self.npsel = 0
        self.nedge_sel = 0
        self.indp = 0
        self.state = 'Init'
        self.ax.axis('tight')
        self.show()

    def show(self, clear=False, dnodes=True, dedges=True, dlabels=False, font_size=10, title='Init'):
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
        #xmin,xmax,ymin,ymax = self.L.ax
        self.L.display['clear'] = clear
        self.L.display['nodes'] = dnodes
        self.L.display['edges'] = dedges
        self.L.display['ndlabel'] = dlabels
        self.L.display['edlabel'] = dlabels
        self.L.display['fontsize'] = font_size
        self.L.display['title'] = title
        self.L.display['ednodes'] = True
        self.L.showGs()

    def OnResize(self, event):
        """ OnResize
        """
        print "width", event.width
        print "height", event.height

    def OnPress(self, event):
        """
        OnPress(event)
        """
        self.nsel = 0
        self.ptsel = np.array([])
        self.evt = event.key
        print "Evenement :", self.evt
        self.call_editor()

    def OnClick(self, event):
        """
        OnClick(event)
        """
        self.nsel = 0
        self.ptsel = np.array([])
        xmin, xmax, ymin, ymax = plt.axis()
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

        print
        "Selected point coord : ", self.ptsel
        print "Selected point number: ", self.nsel
        if self.nsel > 0:
            print "Selected Edge : ", self.nsel

        self.call_editor()

    def call_editor(self):
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
        'p'  : toggle to point creation mode
        'z'  : change display parameters
        'q'  : quit interactive mode
        'x'  : save .str2 file
        'w'  : display all layers

            <----p----> CP
            <---p-----
            ----ml---> SP1
            <---ml---- SP1
                         ----ml-----> SP2
                         <---ml------
                                        -----mr----> SS
            <--------------------lm-----------------
            <--------------------mc-----------------
            <--------------------rc-----------------   Create point on segment
            <--------------------e------------------   edit seg
            ----ml---> SS

        """
        sl = self.L.sl
        cold = pyu.coldict()
        #print "call :",self.L.ax
        #print "In State ",self.state
        #print "In Event ",self.evt
        #
        #
        #
        if self.evt == "l":
            listchoices = self.L.name.keys()
            self.L.display['activelayer'] = multchoicebox(
                'message', 'titre', listchoices)
            return()
        #
        # e : Edit
        #
        if self.evt == 'i':
            self.state = 'Init'
            print "Out State", self.state
            return

        if self.evt == 'e':
            if self.state == 'SS':
                self.L.edit_edge(self.selected_edge)
                self.state = 'Init'
                self.show(clear=True, title='Init')
                print "Out State", self.state
                return
        #
        # h : add subsegment
        #
        if self.evt == 'h':
            if self.state == 'SS':
                self.L.add_subseg(self.selected_edge)
                self.show(clear=True, title='Init')
                self.state = 'Init'
                print "Out State", self.state
                return
        #
        # d : delete subsegment
        #
        if self.evt == 'd':
            if self.state == 'SS':
                self.L.del_subseg(self.selected_edge)
                self.show(clear=True, title='Init')
                self.state = 'Init'
                print "Out State", self.state
                return

        if self.evt == 'c':
            if self.state == 'Init':
                #pt = plt.ginput(4)
                #self.ax = plt.gca()
                x1 = self.ax.get_xlim()
                y1 = self.ax.get_ylim()
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
            #plt.axis(self.L.ax)
            plt.axis('tight')
            self.show(clear=True)
            return
        #
        # o : Toggle overlay
        #
        if self.evt == 'o':
            if self.L.display['overlay']:
                self.L.display['overlay'] = False
            else:
                self.L.display['overlay'] = True
            self.show(clear=True)
            return

        #
        # p : Point
        #
        if self.evt == 'p':
            if self.state == "Init":
                self.state = "CP"
                #self.text.set_text('Create Point (Lclic : same x Rclic same y Cclic Free Point')
                self.show(clear=True, title='Create Point (lclic : same x rclic same y cclic Free Point')
            elif self.state == "CP":
                self.state = "Init"
                #self.text.set_text('Init State')
                self.show(clear=True, title='Init State')
            self.fig.canvas.draw()
            print "Out State", self.state
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
        if self.evt == 'w':
        # display all layer
            self.L.display['activelayer'] = self.L.name.keys()
            self.show(clear=True, title=self.L.display['activelayer'])
            return

        if self.state == 'Init':
            #
            # Initial State
            #        -------> SP1  (lclick)
            #        -------> SS   (rclick)
            #
            self.selected_edge = 0
            self.selected_pt1 = 0
            self.selected_pt2 = 0

            if (self.evt == 'lclic') & (self.nsel < 0):
            # select point 1
                self.state = 'SP1'
                self.selected_pt1 = self.nsel
                self.pt1 = np.array(
                    self.L.Gs.pos[self.nsel]).reshape(2, 1)
                self.pt_previous = self.pt1
                self.p1 = self.ax.plot([self.pt1[0]
                                        ], [self.pt1[1]], 'o', visible=True)
                self.p1[0].set_color('yellow')
                self.p1[0].set_ms(10)
                self.p1[0].set_alpha(0.4)
                title = 'SP1 Node : %d ' % (self.nsel)
                self.show(clear=False, title=title)
                self.L.show_nodes(ndlist=[self.nsel],
                                  alpha=0.5, color='y', size=200)
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return

            if (self.evt == 'lclic') & (self.nsel > 0):
            # select seg
                self.state = 'SS'
                self.selected_edge = self.nsel
                nse = self.nsel
                ta, he = self.L.Gs.neighbors(nse)
                pta = np.array(self.L.Gs.pos[ta])
                phe = np.array(self.L.Gs.pos[he])
                alpha = self.L.display['alpha']
                self.current_layer = self.L.Gs.node[nse]['name']
                #self.seg       = linet(self.ax,pta,phe,alpha,'red',3.5)
                segdico = self.L.Gs.node[nse]
                if 'ss_name' in segdico:
                    cosegname = segdico['ss_name']
                    titre = 'Select Segment : %d (%d->%d) Layer : %s Coseg : %s' % (nse, ta, he, self.current_layer, cosegname)
                else:
                    titre = 'Select Segment : %d (%d->%d) Layer : %s' % (nse, ta, he, self.L.Gs.node[nse]['name'])

                self.show(clear=True, title=titre)
                self.L.show_nodes(ndlist=[nse], size=200, color='r', alpha=0.5)
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            #self.ax.axis(self.L.ax)
            #self.fig.canvas.draw()

        if self.state == 'CP':
            #
            #  State : Create Point
            #
            #
            if (self.evt == 'lclic') & (self.nsel == 0):
            # add point same x
                self.ptsel[0] = self.pt_previous[0]
                self.L.add_fnod(self.ptsel)
                self.show(clear=False)
                self.pt_previous = self.ptsel
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            if (self.evt == 'rclic') & (self.nsel == 0):
            # add point same x
                self.ptsel[1] = self.pt_previous[1]
                self.L.add_fnod(self.ptsel)
                self.show(clear=False)
                self.pt_previous = self.ptsel
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            if (self.evt == 'cclic') & (self.nsel == 0):
            # add point
                self.L.add_fnod(self.ptsel)
                self.pt_previous = self.ptsel
                self.show(clear=False)
                self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return

        if (self.state == 'SP1'):
            # select point 1
            if (self.evt == 'lclic') & (self.nsel < 0):
                if self.nsel == self.selected_pt1:
                # unselect pt1
                    self.state = 'Init'
                    self.selected_pt1 = 0
                    self.p1[0].set_visible(False)
                    self.show(clear=True, title='Init State')
                    #self.ax.axis(self.L.ax)
                    self.fig.canvas.draw()
                    print "Out State", self.state
                    return
                else:
                # select pt2
                    self.state = 'SP2'
                    self.selected_pt2 = self.nsel
                    self.pt2 = np.array(self.L.Gs.pos[self.nsel]).reshape(2, 1)
                    self.p2 = self.ax.plot([self.pt2[0]],
                                           [self.pt2[1]], 'o', visible=True)
                    self.p2[0].set_color('green')
                    self.p2[0].set_ms(10)
                    self.p2[0].set_alpha(0.4)
                    self.p2[0].set_visible(True)
                    self.show(clear=True, title='SP2 : cclic -> delete point / rclic -> add segment Layer %s' % (self.current_layer))
                    self.L.show_nodes(ndlist=[self.selected_pt1, self.selected_pt2], alpha=0.5, color='g', size=200)
                    #self.ax.axis(self.L.ax)
                    self.fig.canvas.draw()
                    print "Out State", self.state
                    return
            if (self.evt == 'cclic') & (self.selected_pt1 == self.nsel):
            # delete selected point
                self.state = 'Init'
                self.selected_pt1 = 0
                self.L.del_node(self.nsel)
                self.p1[0].set_visible(False)
                self.show(clear=True, title='Init State')
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            if (self.evt == 'rclic'):
            # edit point
                #self.L.Gs.edit_node(self.selected_pt1)
                self.state = 'Init'
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return

            #self.ax.axis(self.L.ax)
            #self.fig.canvas.draw()

        if (self.state == 'SP2'):
            #print self.selected_pt1
            #print self.selected_pt2
            # print self.pt1
            #print self.pt2
            if (self.evt == 'lclic') & (self.nsel < 0):
                if self.nsel == self.selected_pt2:
                # unselect pt2
                    self.selected_pt2 = 0
                    self.p2[0].set_visible(False)
                    self.state = 'SP1'
                    #self.ax.axis(self.L.ax)
                    self.fig.canvas.draw()
                    print "Out State", self.state
                    return
                if self.nsel == self.selected_pt1:
                # unselect pt1
                # pt2 --> pt1
                    self.selected_pt1 = self.selected_pt2
                    self.p1[0].set_visible(False)
                    self.p1 = self.p2
                    self.p1[0].set_color('yellow')
                    self.p1[0].set_visible(True)
                    self.state = 'SP1'
                    #self.ax.axis(self.L.ax)
                    self.fig.canvas.draw()
                    print "Out State", self.state
                    return
            if (self.evt == 'lclic') & (self.nsel == 0):
            # unselect pt2
                self.selected_pt2 = 0
                self.p2[0].set_visible(False)
                self.state = 'SP1'
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            if (self.evt == 'cclic') & (self.selected_pt2 == self.nsel):
            # delete pt2
                self.selected_pt2 = 0
                self.L.del_node(self.nsel)
                self.state = 'SP1'
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            if (self.evt == 'rclic'):
            # add segment
                self.state = 'Init'
                n1 = self.selected_pt1
                n2 = self.selected_pt2
                nn1 = np.array(self.L.Gs.neighbors(n1))
                nn2 = np.array(self.L.Gs.neighbors(n2))
                inn1nn2 = np.intersect1d(nn1, nn2)
                # check if segment already exists
                if len(inn1nn2) > 0:
                    bool = True
                    nume = inn1nn2[0]
                else:
                    bool = False
                #bool,iseg=self.L.isseg(self.pt1,self.pt2)
                #print "Le segment existe ? ",bool
                #print "Son numero est : ",iseg
                if bool:
                # segment already exists
                # set selected segment in red
                    self.state = 'SS'
                    ta, he = self.L.Gs.neighbors(nume)
                    segdico = self.L.Gs.node[nume]
                    if 'ss_name' in segdico:
                        cosegname = segdico['ss_name']
                        titre = 'Select Segment : %d (%d->%d) Layer : %s Coseg : %s' % (nume, ta, he, self.current_layer, cosegname)
                    else:
                        titre = 'Select Segment : %d (%d->%d) Layer : %s' % (nume, ta, he, self.L.Gs.node[nume]['name'])
                    self.show(clear=True, title=titre)
                    self.selected_edge = nume
                    self.L.show_nodes(ndlist=[nume],
                                      size=200, color='r', alpha=0.5)
                else:
                # segment creation
                    #del(self.p1)
                    #del(self.p2)
                    self.p1[0].set_visible(False)
                    self.p2[0].set_visible(False)
                    ta = self.selected_pt1
                    he = self.selected_pt2
                    nume = self.L.add_edge(ta, he, name=self.current_layer)
                    self.show(clear=False)
                    self.state = 'SS'
                    titre = 'Select Segment : %d (%d->%d) Layer : %s' % (nume, ta, he, self.L.Gs.node[nume]['name'])
                    self.show(clear=True, title=titre)
                    self.selected_edge = nume
                    self.L.show_nodes(ndlist=[nume],
                                      size=200, color='r', alpha=0.5)
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return

            #self.ax.axis(self.L.ax)
            #self.fig.canvas.draw()

        if (self.state == 'SS'):
            #
            # State Select Segment (SS)
            #
            #  lclick on selected edge : add point on edge
            #  cclick outside selected edge : unselect
            #  rclick on selected edge : edit
            #
            ta, he = self.L.Gs.neighbors(self.selected_edge)
            pta = np.array(self.L.Gs.pos[ta])
            phe = np.array(self.L.Gs.pos[he])
            if (self.evt == 'lclic'):
            # Unselect edge
                self.state = 'Init'
                self.show(clear=True, title='Init')
                print "Out State", self.state
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                return
            if (self.evt == 'cclic'):
                if self.nsel == self.selected_edge:
                # delete edge
                    print "Delete Edge :", self.nsel
                    self.L.del_edge(self.selected_edge)
                    self.show(clear=True)
                    self.state = 'Init'
                    self.text.set_text('Init')
                    print "Out State", self.state
                    #self.ax.axis(self.L.ax)
                    self.fig.canvas.draw()
                    return

            if (self.evt == 'rclic'):
            # create point on edge
                pt_new = geu.ptonseg(pta, phe, self.ptsel)
                pd1 = pt_new - pta
                pd2 = phe - pta
                alpha = np.sqrt(np.dot(pd1, pd1)) / np.sqrt(np.dot(pd2, pd2))
                if (pt_new != []):
                    # calculer alpha
                    self.L.add_none(self.selected_edge, 1. - alpha)
                    self.current_layer = self.L.Gs.node[
                        self.selected_edge]['name']
                    #self.L.dels(self.selected_edge)
                    #self.L.adds([[ta],[self.L.nn-1]],self.current_layer)
                    #self.L.adds([[self.L.nn-1],[he]],self.current_layer)
                    #self.seg.set_visible['False']
                    self.text.set_text('selected none')
                    self.show(clear=True, title='Init')
                    self.state = 'Init'
                else:
                    self.state = 'SS'
                    self.show(clear=True, title="SS : Try again ")

                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return

            if (self.evt == 'alt'):
                print "Edit Coseg"
                print "Out State", self.state
                #self.ax.axis(self.L.ax)
                self.fig.canvas.draw()
                return

            return
