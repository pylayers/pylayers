#!usr/bin/python
# -*- coding: latin1 -*-
import os
import pdb
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
    z : change display parameter (GUI)


    """
    def __init__(self,L,fig,ax):
        """ SelectL is a class which associates a Layout and a figure

        Parameters
        ----------
        L   : Layout
        fig : figure

        """
        self.L = L
        ax.axis(self.L.display['box'])
        plt.title('Init')
        self.text = ax.text(0.05, 0.95, 'Selected : none',
                                 transform=ax.transAxes, va='top')
        self.set_origin = False
        self.set_x = False
        self.set_y = False
        self.pt = []
        self.seg = []
        self.coseg = []
        self.pt1 = np.array([])
        self.pt2 = np.array([])
        self.selected_pt1 = 0
        self.selected_pt2 = 0
        self.selected_edge1 = 0
        self.selected_edge2 = 0
        self.current_layer = self.L.display['activelayer']
        self.npsel = 0
        self.nedge_sel = 0
        self.indp = 0
        self.state = 'Init'
        self.nsel = 0 
        self.update_state()

    def show(self,fig,ax,clear=False, dnodes=True, dedges=True,  font_size=10, title='Init'):
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
        #laxe  = ax.get_axes()
        #xmin,xmax=laxe.xaxis.get_view_interval()
        #ymin,ymax=laxe.yaxis.get_view_interval()
        #xmin,xmax,ymin,ymax = self.g.ax
        #ax = plt.axis('tight')
        axis = ax.axis()
        #print('show : axis',axis) 
        #plt.axis(ax)
        self.L.display['clear'] = clear
        #self.L.display['nodes'] = dnodes
        #self.L.display['edges'] = dedges
        self.L.display['fontsize'] = font_size
        self.L.display['title'] = title
        #self.L.display['ednodes'] = False
        #self.L.display['nodes'] = False
        fig,ax = self.L.showGs(fig,ax,axis=axis)
        return(fig,ax)

    def OnPress(self,event):
        """
        OnPress(event)
        """
        fig = plt.gcf()
        ax  = plt.gca()
        self.nsel = 0
        self.ptsel = np.array([])
        self.evt = event.key
        #print "Evenement :", self.evt
        self.new_state()

    def OnClick(self, event):
        """
        OnClick(event)
        """
        fig = plt.gcf()
        ax  = plt.gca()
        self.nsel = 0
        self.ptsel = np.array([])
        xmin, xmax, ymin, ymax = ax.axis()
        #print xmin,xmax,ymin,ymax
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

        #print "Selected point coord : ", self.ptsel
        #print "Selected point number: ", self.nsel
        if self.nsel > 0:
            print "Selected Edge : ", self.nsel

        self.new_state()

    def update_state(self):
        """ update state
        """
        fig = plt.gcf()
        ax = plt.gca()
        if self.state == 'Init':
            fig,ax = self.show(fig,ax,clear=True)
            self.selected_edge1 = 0
            self.selected_pt1 = 0
            self.selected_pt2 = 0
            try:
                del self.pt_previous
            except:
                pass
            ax.title.set_text('Init : '
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
            fig,ax = self.show(fig,ax,clear=False)
            ax.title.set_text('Selected node : %d ' % (self.nsel))
            self.selected_pt1 = self.nsel
            self.pt1 = np.array(self.L.Gs.pos[self.nsel]).reshape(2, 1)
            self.pt_previous = self.pt1
            self.p1 = ax.plot([self.pt1[0]], [self.pt1[1]], 'o', visible=True)
            self.p1[0].set_color('yellow')
            self.p1[0].set_ms(10)
            self.p1[0].set_alpha(0.4)
            try:
                self.p2.set_visible(False)
            except:
                pass

        if self.state == 'SP2':
            self.p1[0].set_color('green')
            ax.title.set_text('Selected node : %d ' % (self.nsel))
            self.selected_pt2 = self.nsel
            self.pt2 = np.array(self.L.Gs.pos[self.nsel]).reshape(2, 1)
            self.pt_previous = self.pt2
            self.p2 = ax.plot([self.pt2[0]], [self.pt2[1]], 'o', visible=True)
            self.p2[0].set_color('green')
            self.p2[0].set_ms(10)
            self.p2[0].set_alpha(0.4)
            ax.title.set_text('SP2')

        if self.state == 'SS':
            try:
                self.p1[0].set_visible(False)
            except:
                pass
            try:
                self.p2[0].set_visible(False)
            except:
                pass
            self.selected_edge1 = self.nsel
            nse = self.nsel
            ta, he = self.L.Gs.neighbors(nse)
            pta = np.array(self.L.Gs.pos[ta])
            phe = np.array(self.L.Gs.pos[he])
            alpha = self.L.display['alpha']
            self.current_layer = self.L.Gs.node[nse]['name']
            self.L.display['activelayer'] = self.current_layer
            #self.seg       = linet(ax,pta,phe,alpha,'red',3.5)
            segdico = self.L.Gs.node[nse]
            fig,ax=self.show(fig,ax,clear=False)
            self.segment = ax.plot([pta[0],phe[0]],
                                        [pta[1],phe[1]],
                                        'r',linewidth=3, visible=True)
            if 'ss_name' in segdico:
                cosegname = segdico['ss_name']
                titre = 'Select Segment : %d (%d->%d) Layer : %s Coseg : %s ' % (nse, ta, he, self.current_layer, cosegname)
            else:
                titre = 'Select Segment : %d (%d->%d) Layer : %s' % (nse, ta, he, self.L.Gs.node[nse]['name'])

            ax.title.set_text(titre)
            self.L.show_nodes(ndlist=[nse], size=200, color='r', alpha=0.5)

        if self.state == 'SSS':
            nse = self.selected_edge1
            segdico = self.L.Gs.node[nse]
            zmin    = segdico['ss_zmin']
            zmax    = segdico['ss_zmax']
            ax.title.set_text('SSS : '+self.L.Gs.node[nse]['name']+' ['+str(zmin)+','+str(zmax)+']')
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
            fig,ax=self.show(fig,ax,clear=False,title='CP lclic : free point, rclic same x, cclic same y')

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
            self.segment1 = ax.plot([self.pta1[0],self.phe1[0]],
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
            self.segment2 = ax.plot([self.pta2[0],self.phe2[0]],
                                        [self.pta2[1],self.phe2[1]],
                                        'c',linewidth=3, visible=True)



        #print self.state
        #print self.nsel
        #print self.selected_pt1
        #print self.selected_pt2
        fig.canvas.draw()
        return(fig,ax)


    def new_state(self):
        """ layout editor state machine

        Parameters
        ----------

        'l'  : select activelayer
        'i'  : back to init state
        'e'  : edit segment
        't'  : translate  structure
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
        fig = plt.gcf()
        ax  = plt.gca()
        sl = self.L.sl
        cold = pyu.coldict()
        #print "In State ",self.state
        #print "In Event ",self.evt

        #
        # flip layout in y
        #
        if self.evt == 'v':
            for n in self.L.Gs.pos:
                self.L.Gs.pos[n]=(self.L.Gs.pos[n][0],-self.L.Gs.pos[n][1])
            self.update_state()
            return
        #
        # offset layout
        #
        if self.evt == 't':
            offx,offy = offsetbox()
            for n in self.L.Gs.pos:
                self.L.Gs.pos[n]=(self.L.Gs.pos[n][0]+offx,self.L.Gs.pos[n][1]+offy)
            self.update_state()
            return

        #
        # Choose layers to visualized
        #
        if self.evt == 'l':
            listchoices = self.L.name.keys()
            self.L.display['layers'] = multchoicebox('message',
                                                     'titre', listchoices)
            self.state = 'Init'
            self.update_state()
            return
        #
        # 'f' toggle points nodes display
        #
        if self.evt=='f':
            self.L.display['nodes'] = not self.L.display['nodes']
            print self.L.display['nodes']
            self.update_state()
            return

        #
        # 'g' toggle segment nodes dislay
        #
        if self.evt=='g':
            self.L.display['ednodes'] = not self.L.display['ednodes'] 
            print self.L.display['ednodes']
            self.update_state()
            return

        #
        # '=' Increment layer
        #
        if self.evt=='=':
            N = len(self.L.display['layerset'])
            index = self.L.display['layerset'].index(self.L.display['activelayer'])
            self.L.display['activelayer'] = self.L.display['layerset'][(index+1) % N]
            self.current_layer = self.L.display['activelayer']
            self.update_state()
            return

        #
        # '=' Decrement layer
        #
        if self.evt=='$':
            N = len(self.L.display['layerset'])
            index = self.L.display['layerset'].index(self.L.display['activelayer'])
            self.L.display['activelayer'] = self.L.display['layerset'][(index-1) % N]
            self.current_layer = self.L.display['activelayer']
            self.update_state()
            return
        #
        # 'i' : Back to init state 
        #
        if self.evt == 'i':
            self.state = 'Init'
            self.update_state()
            return

        #
        #  'e'
        #       if state == Init
        #           egalize points coordinates
        #
        #       if state == SS
        #           edit segment properties
        #
        if self.evt == 'e':
            if (self.state == 'Init'):
                #
                # averaging one point coordinate along the smallest dimension
                #
                x1 = ax.get_xbound()
                y1 = ax.get_ybound()
                # get node list and edge list
                ndlist, edlist = self.L.get_zone([x1[0],x1[1],y1[0],y1[1]])
                for k,nd in enumerate(ndlist):
                    try:
                        tp = np.vstack((tp,np.array(self.L.Gs.pos[nd])))
                    except:
                        tp = np.array(self.L.Gs.pos[nd])
                mtp = np.sum(tp,axis=0)/(k+1)
                stp = np.sqrt(np.sum((tp-mtp)*(tp-mtp),axis=0)/(k+1))
                # if the standard deviation is lower than 10cm
                # averaging coordinates along the shortest axis
                if min(stp) < 0.10:
                    ind = np.where(stp==min(stp))[0][0]
                    for nd in ndlist:
                        x = self.L.Gs.pos[nd][0]
                        y = self.L.Gs.pos[nd][1]
                        if ind ==0:
                            self.L.Gs.pos[nd]=(mtp[0],y)
                        if ind ==1:
                            self.L.Gs.pos[nd]=(x,mtp[1])
                    plt.axis('tight')
                    fig,ax = self.show(fig,ax,clear=True)
                    self.update_state()
                return()

            if (self.state == 'SS') | (self.state =='SSS'):
                self.L.edit_segment(self.selected_edge1)
                self.state = 'Init'
                self.update_state()
                return

            if self.state == 'SP1':
                self.L.edit_point(self.selected_pt1)
                self.state = 'Init'
                self.update_state()
                return
        #
        # "b" : enter a segment node value with keyboard
        #
        if self.evt == 'b':
            if self.state == 'Init':
                self.nsel = eval(raw_input("seg number :"))
                #self.L.edit_segment(nseg)
                self.state='SS'
                self.update_state()
                return

        #
        # j : vertical and horizontal scaling (Init)
        #
        if self.evt == 'j':
            if self.state == 'Init':
                vscale,hscale = offsetbox(text1='Enter scaling values',
                                          text2=('vscale','hscale'),
                                          default=('1.0','1.0')
                                          )
                for n in self.L.Gs.pos:
                    self.L.Gs.pos[n]=(self.L.Gs.pos[n][0],self.L.Gs.pos[n][1]*vscale)
                plt.axis('tight')
                fig,ax = self.show(fig,ax,clear=True)
                self.update_state()
                return

        # Init
        # h : horizontal scaling factor
        #    add subsegment (SS)
        #
        if self.evt == 'h':
#            if self.state == 'Init':
#                hscale = eval(raw_input("horizontal scaling factor : "))
#                for n in self.L.Gs.pos:
#                    self.L.Gs.pos[n]=(self.L.Gs.pos[n][0]*hscale,self.L.Gs.pos[n][1])
#                plt.axis('tight')
#                fig,ax = self.show(fig,ax,clear=True)
#                self.update_state()
#                return()

            if self.state == 'SS':
                self.L.add_subseg(self.selected_edge1,self.current_layer)
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
                self.L.del_segment(self.selected_edge1)
                self.state = 'Init'
                self.update_state()
                return
            if self.state == 'SSS':
                self.L.del_subseg(self.selected_edge1)
                self.state = 'Init'
                self.update_state()
                return
        #
        # 'Del' delete points in the current axis region
        #
        if self.evt == 'delete':
            if self.state=='Init':
                x1 = ax.get_xbound()
                y1 = ax.get_ybound()
                ndlist, edlist = self.L.get_zone([x1[0],x1[1],y1[0],y1[1]])
                #print x1,y1
                #print ndlist
                self.L.del_node(ndlist)
                self.update_state()
                return

        #
        # r : Refresh
        #
        if self.evt == 'r':
            plt.axis('tight')
            fig,ax = self.show(fig,ax,clear=True)
            self.state = 'Init'
            self.update_state()
            return

        #
        # o : Toggle overlay
        #
        if self.evt == 'o':
            if self.state <> 'CP':
                if self.L.display['overlay']:
                    self.L.display['overlay'] = False
                    self.update_state()
                else:
                    self.L.display['overlay'] = True
                    self.update_state()
                return
            else:
                self.set_origin = True


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
        #
        # 'z' : change display parameters
        #
        if self.evt == 'z':
            self.L.displaygui()
            fig,ax = self.show(fig=fig,ax=ax,clear=True)
            return
        #
        # 'q' : quit interactive mode
        #
        if self.evt == 'q':
            fig.canvas.mpl_disconnect(self.L.cid1)
            fig.canvas.mpl_disconnect(self.L.cid2)
            return

        #
        # 'x' save structure
        #
        if self.evt == 'x':
            racine, ext = os.path.splitext(self.L.filename)
            filename = racine + '.str2'
            fileini = racine + '.ini'
            self.L.savestr2(filename)
            self.L.saveini(fileini)
            print "structure saved in ", filename
            print "structure saved in ", fileini
            return
        #
        # 'n' : toggle node label display
        #
        if self.evt == 'n':
            self.L.display['ndlabel'] = not self.L.display['ndlabel']
            self.L.display['edlabel'] = not self.L.display['edlabel']
            fig,ax = self.show(fig=fig,ax=ax,clear=True, title=self.L.display['activelayer'])
            fig.canvas.draw()
            return
        #
        # "w" : display all layers
        #
        if self.evt == 'w':
        # display all layer
            self.L.display['activelayer'] = self.L.name.keys()
            fig,ax = self.show(fig=fig,ax=ax,clear=True, title=self.L.display['activelayer'])
            return fig,ax
        #
        # Left clic and selected node is a point
        #
        if (self.evt == 'lclic') & (self.nsel < 0):

        #
        # select point 1 : Init -> SP1
        #
            if self.state=='Init':
                self.state = 'SP1'
                self.update_state()
                return
        #
        # select point 2 : SP1 --> SP2
        #

            if self.state=='SP1':
                if self.nsel != self.selected_pt1:
                    self.state = 'SP2'
                    self.update_state()
                    return
        #
        # Create point on selected segment orthogonaly to segment starting in
        # selected point
        #
            if self.state=='SS':
                connect = self.L.Gs.node[self.selected_edge1]['connect']
                if (self.nsel != connect[0]) & (self.nsel != connect[1]): 
                   #self.L.add_nfpe(self.nsel,conn
                   pass

        #
        # Left clic and selected node is a segment
        #

        if (self.evt == 'lclic') & (self.nsel > 0):
            if self.state=='Init':
                self.state = 'SS'
                self.update_state()
                return

            if self.state=='SS':
                self.nsel = self.selected_edge1
                segdico = self.L.Gs.node[self.nsel]
                if 'ss_name' in segdico:
                    self.state = 'SSS'
                else:
                    self.state = 'CPS'
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
            # or set origin
            if self.state == 'CP':
                if self.set_origin:
                    offx = self.ptsel[0]
                    offy = self.ptsel[1]
                    print offx,offy
                    xmin,xmax,ymin,ymax = self.L.display['box']
                    self.L.display['box'] = [xmin-offx,xmax-offx,ymin-offy,ymax-offy]
                    self.set_origin=False
                    self.set_x=True
                    plt.axis('tight')
                    fig,ax = self.show(fig,ax,clear=True)
                    self.update_state()
                    return
                if self.set_x:
                    offx = self.ptsel[0]
                    val  = eval(enterbox('enter x value'))
                    ratio = val/offx
                    print ratio
                    xmin,xmax,ymin,ymax = self.L.display['box']
                    self.L.display['box'] = [ratio*xmin,ratio*xmax,ymin,ymax]
                    self.set_x=False
                    self.set_y=True
                    plt.axis('tight')
                    fig,ax = self.show(fig,ax,clear=True)
                    self.update_state()
                    return
                if self.set_y:
                    offx = self.ptsel[1]
                    val  = eval(enterbox('enter y value'))
                    ratio = val/offx
                    print ratio
                    xmin,xmax,ymin,ymax = self.L.display['box']
                    self.L.display['box'] = [xmin,xmax,ratio*ymin,ratio*ymax]
                    self.set_y=False
                    plt.axis('tight')
                    fig,ax = self.show(fig,ax,clear=True)
                    self.update_state()
                    return
                else:
                    self.L.add_fnod(tuple(self.ptsel))
                    self.pt_previous = self.ptsel
                    self.update_state()

                return
            if self.state == 'SP2':
                ta = self.selected_pt1
                he = self.selected_pt2
                self.nsel  = self.L.add_segment(ta, he,name=self.current_layer)
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
                    self.L.add_pons(self.selected_edge1, 1. - alpha)
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
                    self.L.add_fnod(tuple(self.ptsel))
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
        #
        # right click : back to CP from CPS
        #
            if self.state == 'CPS':
                self.state = 'CP'
                self.update_state()
                return
        #
        # right click : back to CPS from CPSS
        #
            if self.state == 'CPSS':
                self.state = 'CPS'
                self.update_state(fig,ax)
                return
        #
        # Center Clic event
        #
        if (self.evt == 'cclic'):
            if self.state == 'CP':
                try:
                    self.ptsel[1] = self.pt_previous[1]
                    self.L.add_fnod(tuple(self.ptsel))
                    self.pt_previous = self.ptsel
                    self.update_state()
                    return
                except:
                    pass
        #
        # Left clic and selected node is a point
        #


