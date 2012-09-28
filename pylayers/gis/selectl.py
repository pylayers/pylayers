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
    def __init__(self, g, fig):
        """
        Select is a class which associate a Layout and a figure
        """
        self.g = g
        self.fig = fig
        self.spl = self.fig.add_subplot(111)
        #self.title = self.spl.title='Title is not set'
        self.text = self.spl.text(0.05, 0.95, 'selected : none',
                                  transform=self.spl.transAxes, va='top')
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

    def show(self, clear=False, dnodes=True, dedges=True, dlabels=False, font_size=10, title='Init'):
        """ show 

        Parameters
        ----------
        clear     : boolean
        dnodes    : boolean 
        dedges    : boolean 
        dlabels   : boolean 
        font_size :  
        """
        #laxe  = self.spl.get_axes()
        #xmin,xmax=laxe.xaxis.get_view_interval()
        #ymin,ymax=laxe.yaxis.get_view_interval()
        #xmin,xmax,ymin,ymax = self.g.ax
        ax = plt.axis()
        self.g.display['clear'] = clear
        self.g.display['nodes'] = dnodes
        self.g.display['edges'] = dedges
        self.g.display['ndlabel'] = dlabels
        self.g.display['edlabel'] = dlabels
        self.g.display['fontsize'] = font_size
        self.g.display['title'] = title
        self.g.display['ednodes'] = True
        self.g.showGs()
        #plt.axis((xmin,xmax,ymin,ymax))
        plt.axis(ax)

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
            self.nsel = self.g.ispoint(self.ptsel, dd / 100)

        if event.button == 2 and event.inaxes:
            self.evt = 'cclic'
            x = event.xdata
            y = event.ydata
            self.ptsel = np.array((x, y))
            self.nsel = self.g.ispoint(self.ptsel, dd / 100)

        if event.button == 3 and event.inaxes:
            self.evt = 'rclic'
            x = event.xdata
            y = event.ydata
            self.ptsel = np.array((x, y))
            self.nsel = self.g.ispoint(self.ptsel, dd / 100)

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

        Notes
        -----
            ----p----> CP
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
        sl = self.g.sl
        cold = pyu.coldict()
        #print "call :",self.g.ax
        #print "In State ",self.state
        #print "In Event ",self.evt
        #
        #
        #
        if self.evt == "l":
            listchoices = self.g.name.keys()
            self.g.display['activelayer'] = multchoicebox(
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
                self.g.edit_edge(self.selected_edge)
                self.state = 'Init'
                self.show(clear=True, title='Init')
                print "Out State", self.state
                return
        #
        # h : add subsegment
        #
        if self.evt == 'h':
            if self.state == 'SS':
                self.g.add_subseg(self.selected_edge)
                self.show(clear=True, title='Init')
                self.state = 'Init'
                print "Out State", self.state
                return
        #
        # d : delete subsegment
        #
        if self.evt == 'd':
            if self.state == 'SS':
                self.g.del_subseg(self.selected_edge)
                self.show(clear=True, title='Init')
                self.state = 'Init'
                print "Out State", self.state
                return

        if self.evt == 'c':
            if self.state == 'Init':
                ax = plt.axis()
                ndlist,edlist = self.g.get_zone(ax)
                self.g.del_node(ndlist)

        #
        # r : Refresh
        #
        if self.evt == 'r':
            #plt.axis(self.g.ax)
            plt.axis('tight')
            self.show(clear=True)
            return
        #
        # o : Toggle overlay
        #
        if self.evt == 'o':
            if self.g.display['overlay']:
                self.g.display['overlay'] = False
            else:
                self.g.display['overlay'] = True
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
            self.g.displaygui()
            self.show(clear=True)
            return

        if self.evt == 'q':
        # quit interactive mode
            self.fig.canvas.mpl_disconnect(self.g.cid1)
            self.fig.canvas.mpl_disconnect(self.g.cid2)
            return
        if self.evt == 'x':
        # save structure
            racine, ext = os.path.splitext(self.g.filename)
            filename = racine + '.str2'
            self.g.savestr2(filename)
            print "structure saved in ", filename
            return
        if self.evt == 'w':
        # display all layer
            self.g.display['activelayer'] = self.g.name.keys()
            self.show(clear=True, title=self.g.display['activelayer'])
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
                    self.g.Gs.pos[self.nsel]).reshape(2, 1)
                self.pt_previous = self.pt1
                self.p1 = self.spl.plot([self.pt1[0]
                                         ], [self.pt1[1]], 'o', visible=True)
                self.p1[0].set_color('yellow')
                self.p1[0].set_ms(10)
                self.p1[0].set_alpha(0.4)
                title = 'SP1 Node : %d ' % (self.nsel)
                self.show(clear=False, title=title)
                self.g.show_nodes(ndlist=[self.nsel],
                                  alpha=0.5, color='y', size=200)
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return

            if (self.evt == 'lclic') & (self.nsel > 0):
            # select seg
                self.state = 'SS'
                self.selected_edge = self.nsel
                nse = self.nsel
                ta, he = self.g.Gs.neighbors(nse)
                pta = np.array(self.g.Gs.pos[ta])
                phe = np.array(self.g.Gs.pos[he])
                alpha = self.g.display['alpha']
                self.current_layer = self.g.Gs.node[nse]['name']
                #self.seg       = linet(self.spl,pta,phe,alpha,'red',3.5)
                segdico = self.g.Gs.node[nse]
                if 'ss_name' in segdico:
                    cosegname = segdico['ss_name']
                    titre = 'Select Segment : %d (%d->%d) Layer : %s Coseg : %s' % (nse, ta, he, self.current_layer, cosegname)
                else:
                    titre = 'Select Segment : %d (%d->%d) Layer : %s' % (nse, ta, he, self.g.Gs.node[nse]['name'])

                self.show(clear=True, title=titre)
                self.g.show_nodes(ndlist=[nse], size=200, color='r', alpha=0.5)
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            #self.spl.axis(self.g.ax)
            #self.fig.canvas.draw()

        if self.state == 'CP':
            #
            #  State : Create Point
            #
            #
            if (self.evt == 'lclic') & (self.nsel == 0):
            # add point same x
                self.ptsel[0] = self.pt_previous[0]
                self.g.add_fnod(self.ptsel)
                self.show(clear=False)
                self.pt_previous = self.ptsel
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            if (self.evt == 'rclic') & (self.nsel == 0):
            # add point same x
                self.ptsel[1] = self.pt_previous[1]
                self.g.add_fnod(self.ptsel)
                self.show(clear=False)
                self.pt_previous = self.ptsel
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            if (self.evt == 'cclic') & (self.nsel == 0):
            # add point
                self.g.add_fnod(self.ptsel)
                self.pt_previous = self.ptsel
                self.show(clear=False)
                self.spl.axis(self.g.ax)
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
                    #self.spl.axis(self.g.ax)
                    self.fig.canvas.draw()
                    print "Out State", self.state
                    return
                else:
                # select pt2
                    self.state = 'SP2'
                    self.selected_pt2 = self.nsel
                    self.pt2 = np.array(self.g.Gs.pos[self.nsel]).reshape(2, 1)
                    self.p2 = self.spl.plot([self.pt2[0]],
                                            [self.pt2[1]], 'o', visible=True)
                    self.p2[0].set_color('green')
                    self.p2[0].set_ms(10)
                    self.p2[0].set_alpha(0.4)
                    self.p2[0].set_visible(True)
                    self.show(clear=True, title='SP2 : cclic -> delete point / rclic -> add segment Layer %s' % (self.current_layer))
                    self.g.show_nodes(ndlist=[self.selected_pt1, self.selected_pt2], alpha=0.5, color='g', size=200)
                    #self.spl.axis(self.g.ax)
                    self.fig.canvas.draw()
                    print "Out State", self.state
                    return
            if (self.evt == 'cclic') & (self.selected_pt1 == self.nsel):
            # delete selected point
                self.state = 'Init'
                self.selected_pt1 = 0
                self.g.del_node(self.nsel)
                self.p1[0].set_visible(False)
                self.show(clear=True, title='Init State')
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            if (self.evt == 'rclic'):
            # edit point
                #self.g.Gs.edit_node(self.selected_pt1)
                self.state = 'Init'
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return

            #self.spl.axis(self.g.ax)
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
                    #self.spl.axis(self.g.ax)
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
                    #self.spl.axis(self.g.ax)
                    self.fig.canvas.draw()
                    print "Out State", self.state
                    return
            if (self.evt == 'lclic') & (self.nsel == 0):
            # unselect pt2
                self.selected_pt2 = 0
                self.p2[0].set_visible(False)
                self.state = 'SP1'
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            if (self.evt == 'cclic') & (self.selected_pt2 == self.nsel):
            # delete pt2
                self.selected_pt2 = 0
                self.g.del_node(self.nsel)
                self.state = 'SP1'
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return
            if (self.evt == 'rclic'):
            # add segment
                self.state = 'Init'
                n1 = self.selected_pt1
                n2 = self.selected_pt2
                nn1 = np.array(self.g.Gs.neighbors(n1))
                nn2 = np.array(self.g.Gs.neighbors(n2))
                inn1nn2 = np.intersect1d(nn1, nn2)
                # check if segment already exists
                if len(inn1nn2) > 0:
                    bool = True
                    nume = inn1nn2[0]
                else:
                    bool = False
                #bool,iseg=self.g.isseg(self.pt1,self.pt2)
                #print "Le segment existe ? ",bool
                #print "Son numero est : ",iseg
                if bool:
                # segment already exists
                # set selected segment in red
                    self.state = 'SS'
                    ta, he = self.g.Gs.neighbors(nume)
                    segdico = self.g.Gs.node[nume]
                    if 'ss_name' in segdico:
                        cosegname = segdico['ss_name']
                        titre = 'Select Segment : %d (%d->%d) Layer : %s Coseg : %s' % (nume, ta, he, self.current_layer, cosegname)
                    else:
                        titre = 'Select Segment : %d (%d->%d) Layer : %s' % (nume, ta, he, self.g.Gs.node[nume]['name'])
                    self.show(clear=True, title=titre)
                    self.selected_edge = nume
                    self.g.show_nodes(ndlist=[nume],
                                      size=200, color='r', alpha=0.5)
                else:
                # segment creation
                    #del(self.p1)
                    #del(self.p2)
                    self.p1[0].set_visible(False)
                    self.p2[0].set_visible(False)
                    ta = self.selected_pt1
                    he = self.selected_pt2
                    nume = self.g.add_edge(ta, he, name=self.current_layer)
                    self.show(clear=False)
                    self.state = 'SS'
                    titre = 'Select Segment : %d (%d->%d) Layer : %s' % (nume, ta, he, self.g.Gs.node[nume]['name'])
                    self.show(clear=True, title=titre)
                    self.selected_edge = nume
                    self.g.show_nodes(ndlist=[nume],
                                      size=200, color='r', alpha=0.5)
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return

            #self.spl.axis(self.g.ax)
            #self.fig.canvas.draw()

        if (self.state == 'SS'):
            #
            # State Select Segment (SS)
            #
            #  lclick on selected edge : add point on edge
            #  cclick outside selected edge : unselect
            #  rclick on selected edge : edit
            #
            ta, he = self.g.Gs.neighbors(self.selected_edge)
            pta = np.array(self.g.Gs.pos[ta])
            phe = np.array(self.g.Gs.pos[he])
            if (self.evt == 'lclic'):
            # Unselect edge
                self.state = 'Init'
                self.show(clear=True, title='Init')
                print "Out State", self.state
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                return
            if (self.evt == 'cclic'):
                if self.nsel == self.selected_edge:
                # delete edge
                    print "Delete Edge :", self.nsel
                    self.g.del_edge(self.selected_edge)
                    self.show(clear=True)
                    self.state = 'Init'
                    self.text.set_text('Init')
                    print "Out State", self.state
                    #self.spl.axis(self.g.ax)
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
                    self.g.add_none(self.selected_edge, 1. - alpha)
                    self.current_layer = self.g.Gs.node[
                        self.selected_edge]['name']
                    #self.g.dels(self.selected_edge)
                    #self.g.adds([[ta],[self.g.nn-1]],self.current_layer)
                    #self.g.adds([[self.g.nn-1],[he]],self.current_layer)
                    #self.seg.set_visible['False']
                    self.text.set_text('selected none')
                    self.show(clear=True, title='Init')
                    self.state = 'Init'
                else:
                    self.state = 'SS'
                    self.show(clear=True, title="SS : Try again ")

                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                print "Out State", self.state
                return

            if (self.evt == 'alt'):
                print "Edit Coseg"
                print "Out State", self.state
                #self.spl.axis(self.g.ax)
                self.fig.canvas.draw()
                return

            return
