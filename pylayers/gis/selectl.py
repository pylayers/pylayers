#!usr/bin/python
# -*- coding: latin1 -*-
import os
import pdb
import Image
import numpy as np
from pylayers.util import geomutil as geu
from pylayers.util import pyutil as pyu
import pylayers.util.plotutil as plu

import matplotlib.pyplot as plt
from pylayers.util.easygui import *
from matplotlib.widgets import RectangleSelector

import copy


class SelectL(object):
    """ Associates a Layout and a figure

    'l'  : select activelayer
    'i'  : back to init state
    'e'  : edit segment
    'CTLR + t'  : translate  structure
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
    def __init__(self,L,fig,ax):
        """ SelectL is a class which associates a Layout and a figure

        Parameters
        ----------

        L   : Layout
        fig : figure
        ax  : axes

        """
        self.fig=fig
        self.ax=ax

        self.L = L
        self.text = self.ax.text(0.05, 0.95, 'Selected : none',
                                 transform=self.ax.transAxes, va='top')
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
        self.evt=''
        self.statename={'Init':'Point/Segment Selection',
                'CP':'Create Point',
                'SP1':'Select Point 1',
                'SP2':'Select Point 2, Click Again for Creating Segment',
                'SS':'Select Segment',
                'SSS':'Select Sub Segment',
                'CPS':'Click again for Split Segment',
                'CPSS':'Create Point On Sub Segment',
                'SMP': 'Multiple Points Selection',
                'SMS': 'Multiple Segments Selection'
                }
        self.help={'':'',
                'Init':'Select Point or Segment/ F1: Multiple selection/ F2: Create Point/ CTRL+q: Quit',
                'CP':'Create Point/ +CTRL same x/ +SHIFT same y',
                'SP1':'Select Point/ Click another point to create segment',
                'SP2':'Click Again for Creating Segment',
                'SS':'e: edit segment properties, h: add a sub-segment',
                'SSS':'Select Sub Segment',
                'CPS':'Click again for Split Segment',
                'CPSS':'Create Point On Sub Segment',
                'SMP': 't: toggle point/segment, Shift + select : add selected points, CTRL + select : remove selected points',
                'SMS': 't: toggle point/segment, e: Edit Selected Segments Propeties'
                }
        self.nsel = 0
        self.ax.axis(self.L.display['box'])
        plt.title(self.statename[self.state])
        self.undoGs=[]
        self.bundo=False
        self.update_state()
        self.shift_is_held = False
        self.ctrl_is_held = False
        self.alt_is_held = False
        self.selectpt=[]
        self.selectseg=[]
        self.selected='pt'
        # save matplotlib config
        self.rcconf = {}
        self.rcconf['keymap.save']= plt.rcParams['keymap.save']
        plt.rcParams['keymap.save']=[]


        self.ddoc = {'l'  : 'select activelayer',
            'i'  :' back to init state',
            'j'  :' vertical and horizontal scaling',
            'e'  :' edit segment',
            'b'  :' edit segment keyboard',
            'CTRL + t'  :' translate  structure',
            'h'  :' add subsegment',
            'd'  :' delete selected object',
            'r'  :' refresh',
            'o'  :' toggle overlay (<> CP mode) set origin (CP mode) ',
            'm'  :' toggle mode (point or segment)',
            'n'  : 'toggle node label display ',
            'z'  : 'change display parameters',
            'x'  : 'save .str2 and .ini file',
            'w'  :' display all layers',
            'v'  :' flip layout w.r.t y axis',
            'f'  :' toggle points nodes display',
            'g'  :' toggle segments nodes display',
            '='  :' increment layer ',
            ','  : 'this help',
            'delete' :'delete selected',
            '$'  :' decrement layer '}

        

    def show(self,fig,ax,clear=False, dnodes=True, dedges=True,  font_size=14, title=''):
        """ show layout

        Parameters
        ----------
        clear     : boolean
        dnodes    : boolean
        dedges    : boolean
        dlabels   : boolean
        font_size : integer
        title     : string

        """
        if title=='':
            title = self.statename[self.state]
        axis = self.ax.axis()
        self.L.display['clear'] = clear
        self.L.display['fontsize'] = font_size
        self.L.display['title'] = title
        self.fig,self.ax = self.L.showGs(fig=self.fig,ax=self.ax,axis=axis,subsegnb=True)
        return(self.fig,self.ax)


    def plotselptseg(self,pt,color='y',ms=10,marker='o'):
        """ plot selected point or segments

        Parameters
        ----------

            pt : list
            list of points or segmetns to plot
        """
        if len(pt)>0:
            
            pts = np.array([self.L.Gs.pos[x] for x in pt])
            p1 = self.ax.plot(pts[:,0], pts[:,1],marker=marker, 
                                visible=True, 
                                color =color,
                                ms=10,
                                alpha=0.4)
            self.fig.canvas.draw()

        return self.fig,self.ax


    def OnPress(self,event,verbose=True):
        """ Keyboard event handler


        Parameters
        ----------

        event
        verbose

        """

        # fig = plt.gcf()
        # ax  = plt.gca()
        # selected

        self.nsel = 0
        self.ptsel = np.array([])
        self.evt = event.key
        if event.key == 'shift':
            self.shift_is_held = True
        if event.key == 'control':
            self.ctrl_is_held = True
        if event.key == 'alt':
            self.alt_is_held = True


        if verbose:
            try:
                print "Evenement :", self.evt,self.ddoc[self.evt]
            except:
                print self.evt+ 'N/A'
        self.new_state()


    def OnRelease(self, event):
        if event.key == 'shift':
           self.shift_is_held = False
        if event.key == 'control':
            self.ctrl_is_held = False
        if event.key == 'alt':
            self.alt_is_held = False

    def OnClick(self, event):
        """ handle OnClick event

        Parameters
        ----------

        event :

        See Also
        --------

        pylayers.gis.layout.Layout.ispoint

        """
        fig = self.fig#plt.gcf()
        ax  = self.ax#plt.gca()
        self.nsel = 0
        self.ptsel = np.array([])
        xmin, xmax, ymin, ymax = self.ax.axis()
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
            print "Selected segment : ", self.nsel

        self.new_state()




    def format_coord(self,x, y):
        col = int(x+0.5)
        row = int(y+0.5)
        string = 'x=%1.4f, y=%1.4f'%(x, y)
        try:
            string = string + ' ' + self.L.Gs.node[self.nsel]['name']
        except:
            pass
        try:
            string = string + ' with ' +str(len(self.L.Gs.node[self.nsel]['ss_name'])) + 'subseg(s)'
        except:
            pass
        string = string + ' ///' +self.help[self.state]
        return string
        
        # if col>=0 and col<numcols and row>=0 and row<numrows:
        #     z = X[row,col]
        #     return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
        # else:
        #     return 'x=%1.4f, y=%1.4f'%(x, y)

    def update_state(self):
        """ update state
        """
        # fig = plt.gcf()
        # ax = plt.gca()
        if not self.bundo:
            self.undoGs.append(self.L.Gs.copy())
            if len(self.undoGs) > 50:
                self.undoGs.pop(0)

        self.ax.format_coord=self.format_coord

        if self.state == 'Init':
            self.fig,self.ax = self.show(self.fig,self.ax,clear=True)
            self.ax.title.set_text(self.statename[self.state])
            self.selected_edge1 = 0
            self.selected_pt1 = 0
            self.selected_pt2 = 0
            self.selectpt=[]
            self.selectseg=[]
            try:
                del self.pt_previous
            except:
                pass
            try:
                self.selector.set_active(False)
                print 'inhib select'
            except:
                pass
            #ax.title.set_text(self.state)
            #ax.title.set_text('Init : '
            #                       +self.L.display['activelayer'])
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
            if self.L.Np==0:
                self.state='CP'
                self.update_state()

        if self.state == 'SP1':
            self.fig,self.ax = self.show(self.fig,self.ax,clear=False)
            self.ax.title.set_text(self.statename[self.state])
            print 'Selected node : '+str(self.nsel)
            #ax.title.set_text(self.nsel))
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
            self.ax.title.set_text(self.statename[self.state])
            #ax.title.set_text('Selected node : %d ' % (self.nsel))
            print 'Selected node : ' + str(self.nsel)
            self.selected_pt2 = self.nsel
            self.pt2 = np.array(self.L.Gs.pos[self.nsel]).reshape(2, 1)
            self.pt_previous = self.pt2
            self.p2 = self.ax.plot([self.pt2[0]], [self.pt2[1]], 'o', visible=True)
            self.p2[0].set_color('green')
            self.p2[0].set_ms(10)
            self.p2[0].set_alpha(0.4)
            #ax.title.set_text('SP2')

        if self.state == 'SS':
            self.ax.title.set_text(self.statename[self.state])
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
            self.fig,self.ax=self.show(self.fig,self.ax,clear=False)
            self.segment = self.ax.plot([pta[0],phe[0]],
                                        [pta[1],phe[1]],
                                        'r',linewidth=3, visible=True)
            if 'ss_name' in segdico:
                cosegname = segdico['ss_name']
                titre = 'Select Segment : %d (%d->%d) Layer : %s Coseg : %s ' % (nse, ta, he, self.current_layer, cosegname)
            else:
                titre = 'Select Segment : %d (%d->%d) Layer : %s' % (nse, ta, he, self.L.Gs.node[nse]['name'])
            print titre 
            #ax.title.set_text(titre)
            self.L.show_nodes(ndlist=[nse], size=200, color='r', alpha=0.5)

        if self.state == 'SSS':
            self.ax.title.set_text(self.statename[self.state])
            nse = self.selected_edge1
            segdico = self.L.Gs.node[nse]
            z  = segdico['ss_z']
            #ax.title.set_text('SSS : '+self.L.Gs.node[nse]['name']+' ['+str(z[0])+']')
            print self.L.Gs.node[nse]['name']+' ['+str(z[0])+']'
            self.segment[0].set_color('blue')
        #
        # Create Point state
        #
        if self.state == 'CP':
            self.ax.title.set_text(self.statename[self.state])
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
            print 'lclic : free point, +CTRL same x, +SHIFT: same y'
            self.fig,self.ax=self.show(self.fig,self.ax,clear=False) 
            self.L.g2npy()

        #
        # Create Point on Segment state
        #

        if self.state == 'CPS':
            self.ax.title.set_text(self.statename[self.state])
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
            self.ax.title.set_text(self.statename[self.state])
            self.selected_edge2 = self.nsel
            ta, he = self.L.Gs.neighbors(self.nsel)
            self.pta2 = np.array(self.L.Gs.pos[ta])
            self.phe2 = np.array(self.L.Gs.pos[he])
            self.current_layer = self.L.Gs.node[self.nsel]['name']
            self.L.display['activelayer'] = self.current_layer
            self.segment2 = self.ax.plot([self.pta2[0],self.phe2[0]],
                                        [self.pta2[1],self.phe2[1]],
                                        'c',linewidth=3, visible=True)


        if 'SM' in self.state:
            self.fig,self.ax = self.show(self.fig,self.ax,clear=True)
            self.ax.title.set_text(self.statename[self.state])
            self.selected_edge1 = 0
            self.selected_pt1 = 0
            self.selected_pt2 = 0
            try:
                del self.pt_previous
            except:
                pass
            # self.state='SM'
        #print self.state
        #print self.nsel
        #print self.selected_pt1
        #print self.selected_pt2
        self.fig.canvas.draw()
        return(self.fig,self.ax)


    def new_state(self):
        """ layout editor state machine

        Parameters
        ----------


        'l'  : select activelayer
        'i'  : back to init state
        'j'  : vertical and horizontal scaling
        'e'  : edit segment
        'b'  : edit segment keyboard
        'CTRL + t'  : translate  structure
        'h'  : add subsegment
        'd | Del'  : delete subsegment
        'r | F5'  : refresh
        'o'  : toggle overlay (<> CP mode)
               set origin (CP mode) 
        'm'  : toggle mode (point or segment)
        'n'  : toggle node label display 
        'z'  : change display parameters
        'CTRL+q'  : quit
        'x | CTRL+s'  : save .str2 and .ini file
        'w'  : display all layers
        'v'  : flip layout w.r.t y axis
        'f'  : toggle points nodes display
        'g'  : toggle segments nodes display
        '='  : increment layer 
        '$'  : decrement layer 
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
        if self.evt == ',':
            for k in self.ddoc.keys():
                print k,self.ddoc[k]

        if self.evt == 'v':
            for n in self.L.Gs.pos:
                self.L.Gs.pos[n]=(self.L.Gs.pos[n][0],-self.L.Gs.pos[n][1])
            self.update_state()
            return
        #
        # translation of layout (open a box)
        #
        # if self.evt == 't' :
        #     offx,offy = offsetbox()
        #     for n in self.L.Gs.pos:
        #         self.L.Gs.pos[n]=(self.L.Gs.pos[n][0]+offx,self.L.Gs.pos[n][1]+offy)
        #     self.update_state()
        #     return

        if self.evt=='escape':
            self.state='Init'
            self.update_state()
            self.fig.canvas.draw()
            return


        if self.evt=='ctrl+z':
            self.bundo=True
            print len(self.L.Gs)
            if len (self.undoGs) >2:
                oGs=self.undoGs.pop(-1)
                oGs=self.undoGs.pop(-1)
                self.L.Gs=oGs
                self.L.g2npy()
            self.update_state()
            self.bundo=False
            return

        if self.evt=='t':
            if 'SM' in self.state:
                self.update_state()
                # fig=plt.gcf()
                # ax=plt.gca()
                if self.selected == 'pt':
                    self.plotselptseg(self.selectseg,color='r')
                    PP=self.L.pt[:,self.L.tahe[:,self.L.tgs[self.selectseg]]]
                    if PP.shape[-1]!=0:
                        self.fig,self.ax=plu.displot(PP[:,0],PP[:,1],fig=self.fig,ax=self.ax,color='r',linewidth=3,alpha=0.4)
                        plt.draw()
                    self.selected='seg'
                    self.state='SMS'
                else: 
                    self.fig,self.ax= self.plotselptseg(self.selectpt)
                    self.selected='pt'
                    self.state='SMP'
                self.ax.title.set_text(self.statename[self.state])
                # self.update_state()

        if self.evt == '3':
            self.L._show3()
            return

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
            print self.current_layer
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
                x1 = self.ax.get_xbound()
                y1 = self.ax.get_ybound()
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
                    self.fig,self.ax = self.show(self.fig,self.ax,clear=True)
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
            if self.state == 'SMS':
                outdata=self.L.edit_segment(self.selectseg[0])
                [self.L.edit_segment(s,outdata=outdata,gui=False) for s in self.selectseg]
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
                self.fig,self.ax = self.show(self.fig,self.ax,clear=True)
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
                result = self.L.add_subseg(self.selected_edge1,self.current_layer)
                if result:
                    self.state = 'SSS'
                else :
                    self.state = 'Init'
                self.update_state()
                return
        #
        # d : delete
        #
        if self.evt == 'd' or self.evt =='delete':
            if  self.state == 'SP1':
                self.state = 'Init'
                self.L.del_points(self.selected_pt1)
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

            if self.state=='SMP':
                # get boundary of the region 
                if hasattr(self,'selectpt'):

                    ptlist = self.selectpt
                    self.selectpt=[]
                    self.selectseg=[]
                    self.L.del_points(ptlist)
                    self.state = 'Init'
                    self.update_state()
                    return
                else :
                    print 'no selected region'

            if self.state=='SMS':
                seglist = self.selectseg
                self.selectpt=[]
                self.selectseg=[]
                self.L.del_segment(seglist)
                self.state = 'Init'
                self.update_state()
                return
            else :
                print 'no selected region'
        #
        # r : Refresh
        #
        if self.evt == 'r' or self.evt == 'f5':
            #plt.axis('tight')
            plt.axis(self.L.display['box'])
            self.fig,self.ax = self.show(self.fig,self.ax,clear=True)
            self.state = 'Init'
            self.update_state()
            return

        #
        # o : Toggle overlay
        #
        if self.evt == 'o' and not self.ctrl_is_held:
            self.state='Init'
            self.update_state()
            if self.L.display['overlay']:
                self.L.display['overlay'] = False
                self.update_state()
            else:
                self.L.display['overlay'] = True
                self.update_state()
            return

        if self.evt == 'o' :
            self.set_origin = True


        if self.evt == 'f2':
            self.state = "CP"
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
        #
        # 'z' : change display parameters
        #
        if self.evt == 'z':
            self.L.displaygui()
            self.fig,self.ax = self.show(fig=self.fig,ax=self.ax,clear=True)
            return
        #
        # 'q' : quit interactive mode
        #
        # if self.evt == 'q':
        #     plt.rcParams.update(self.rcconf)
        #     fig.canvas.mpl_disconnect(self.L.cid1)
        #     fig.canvas.mpl_disconnect(self.L.cid2)
        #     return

        if self.evt == 'ctrl+q':
            plt.rcParams.update(self.rcconf)
            self.fig.canvas.mpl_disconnect(self.L.cid1)
            self.fig.canvas.mpl_disconnect(self.L.cid2)
            plt.close()
            return

        #
        # 'x' save structure
        #
        if self.evt == 'x' or self.evt =='ctrl+s':
            racine, ext = os.path.splitext(self.L.filename)
            filename = racine + '.str2'
            fileini = racine + '.ini'

            # Commented because ss_ce not updated 
            #self.L.savestr2(filename)

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
            print self.L.display['activelayer']
            self.fig,ax = self.show(fig=self.fig,ax=self.ax,clear=True)
            self.fig.canvas.draw()
            return
        #
        # "w" : display all layers
        #
        if self.evt == 'w':
        # display all layer
            self.L.display['activelayer'] = self.L.name.keys()
            print self.L.display['activelayer']
            self.fig,self.ax = self.show(fig=self.fig,ax=self.ax,clear=True)
            return self.fig,self.ax
        #
        # Left clic and selected node is a point
        #
        if (self.evt == 'lclic') & (self.nsel < 0):

        #
        # select point 1 : Init -> SP1
        #
            if self.state=='Init':
                # yellow point 
                self.state = 'SP1'
                self.update_state()
                return
        #
        # select point 2 : SP1 --> SP2
        #

            if self.state=='SP1':
                if self.nsel != self.selected_pt1:
                    # green point 
                    self.state = 'SP2'
                    self.update_state()
                    return
                else:
                    self.state = 'Init'
                    # yellow point 
                    self.update_state()
                    return
        #
        # Create point on selected segment orthogonaly to segment starting in
        # selected point
        # 
        # Not finished 
        #
            if self.state=='SS':
                # get the connection of the selected segment
                connect = self.L.Gs.node[self.selected_edge1]['connect']
                if (self.nsel != connect[0]) & (self.nsel != connect[1]): 
                   self.L.add_nfpe(self.nsel,self.nsel,self.selected_edge1,self.selected_edge2)
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
        if (self.evt == 'lclic') and not (self.shift_is_held or self.alt_is_held or self.ctrl_is_held ):
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
                    self.fig,self.ax = self.show(self.fig,self.ax,clear=True)
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
                    self.fig,self.ax = self.show(self.fig,self.ax,clear=True)
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
                    self.fig,self.ax = self.show(self.fig,self.ax,clear=True)
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

                segexist = self.L.isseg(ta,he)
                print segexist
                # if segment do not already exist, create it
                if not segexist: 
                    self.nsel  = self.L.add_segment(ta, he,name=self.current_layer)
                else:
                    print "segment ("+str(ta)+","+str(he)+") already exists"
                self.L.g2npy()
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
        if (self.evt == 'rclic') or (self.evt == 'lclic' and self.ctrl_is_held ):
            if self.state == 'CP':
                try:
                    self.ptsel[0] = self.pt_previous[0]
                    self.L.add_fnod(tuple(self.ptsel))
                    self.pt_previous = self.ptsel
                    self.update_state()
                    return
                except:
                    return

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
        # right click : back to SS from CPS
        #
            if self.state == 'CPS':
                self.state = 'SS'
                self.update_state()
                return
        #
        # right click : back to CPS from CPSS
        #
            if self.state == 'CPSS':
                self.state = 'CPS'
                self.update_state(self.fig,self.ax)
                return
        #
        # Center Clic event
        #
        if (self.evt == 'cclic') or (self.evt == 'lclic' and self.shift_is_held ):
            if self.state == 'CP':
                try:
                    self.ptsel[1] = self.pt_previous[1]
                    self.L.add_fnod(tuple(self.ptsel))
                    self.pt_previous = self.ptsel
                    self.update_state()
                    return
                except:
                    return
        #
        # Left clic and selected node is a point
        #


        def point_select_callback(eclick, erelease):
            'eclick and erelease are the press and release events'
            self.update_state()
            if not (self.shift_is_held or self.ctrl_is_held):
                self.selectpt=[]
                self.selectseg=[]
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata

            # print x1,x2,y1,y2
            if x1>x2:
                x1,x2=x2,x1
            if y1>y2:
                y1,y2=y2,y1
            
            
            # try:
            selectpt,selectseg = self.L.get_zone([x1,x2,y1,y2])

            if not self.ctrl_is_held:
                self.selectpt.extend(selectpt)
                self.selectseg.extend(selectseg)
                self.selectseg=filter(lambda x: self.L.Gs.node[x]['connect'][0] in self.selectpt
                                 and self.L.Gs.node[x]['connect'][1] in self.selectpt,
                                 self.selectseg)

                self.selectpt=np.unique(self.selectpt).tolist()
                self.selectseg=np.unique(self.selectseg).tolist()
            else: 
                [self.selectpt.pop(self.selectpt.index(x)) for x in selectpt if x in self.selectpt]
                [self.selectseg.pop(self.selectseg.index(x)) for x in selectseg if x in self.selectseg]
            # except:
            #     print 'empty selection'
            print self.selectpt,self.selectseg
            self.plotselptseg(self.selectpt)
            self.selected='pt'
            print self.state
            
                

        def toggle_selector(event):
            if toggle_selector.RS.active:
                toggle_selector.RS.set_active(False)
            if not toggle_selector.RS.active:
                toggle_selector.RS.set_active(True)
        
        if self.evt == 'f1':
            # avoid conflict between zoom and selection 
            # fm=plt.get_current_fig_manager()
            # if fm.toolbar._active == 'PAN':
            #     fm.toolbar.pan()
            # if fm.toolbar._active == 'ZOOM':
            #     fm.toolbar.zoom()

            self.state='SMP'
            toggle_selector.RS = RectangleSelector(self.ax, point_select_callback,
                                               drawtype='box', useblit=True,
                                               button=[1,3], # don't use middle button
                                               minspanx=5, minspany=5,
                                               spancoords='pixels')
            self.selector = toggle_selector.RS
            self.update_state()

        if self.evt == 'f9':
            print self.selectpt, self.selectseg
            #print self.selectsl
            # plt.connect('key_press_event', toggle_selector)
 