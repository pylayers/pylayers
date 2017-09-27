# -*- coding: utf-8 -*-
import sys, os, random
try:
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi import mlab
except:
    'Mayavi not installed'
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from pylayers.gis.layout import *
from pylayers.gui.editor_select import SelectL2
from pylayers.util.project import *
import pylayers.util.pyutil as pyu
import os
import sys



class SubSegWin(QDialog):    # any super class is okay
    def __init__(self,idx,subsegdata={},parent=None):
        super(SubSegWin, self).__init__(parent)
        #
        self.gparent=parent.parent
        self.parent=parent
        self.idx=idx
        # mulsti segment selection indicator
        self.mulseg=parent.mulseg
        # dictionnary to pass subseg data
        self.subsegdata=subsegdata
        # self.Nss=Nss
        # self.zmin=zmin
        # self.zmax=zmax
        # # self._autocalc_height_val()
        self._init_subseg_prop()
        self._init_layout()

    def _init_layout(self):

        vbox = QVBoxLayout()

        # # Indicate Ceil
        # hboxceil = QHBoxLayout()
        # ceillabel = QLabel('Ceil')
        # ceillabel.setStyleSheet("font: bold 14px;")
        # hboxceil.addWidget(ceillabel)
        # hboxceil.setAlignment(Qt.AlignCenter)

        # vbox.addLayout(hboxceil)

        # vbox.addWidget(self.Hline())
        vbox.addWidget(self.Hline())


            #slab
            # hboxtitle=QHBoxLayout()
            # hboxtitle.addWidget(QLabel('Sub-Segment'+str(ss+1)))
        hbox1 = QHBoxLayout()
        hbox1.addWidget(self.comboslab)

        # slab prop
        hboxl2 = QHBoxLayout()
        hbox2 = QHBoxLayout()
        label=['zmin','zmax','offset','']
        for iw,w in enumerate([ self.heightmin,self.heightmax,
                                self.offset]):
            hboxl2.addWidget(QLabel(label[iw]))
            hbox2.addWidget(w)

        # vbox.addLayout(hboxtitle)
        vbox.addLayout(hbox1)
        vbox.addLayout(hboxl2)
        vbox.addLayout(hbox2)

        # vbox.addWidget(self.Hline())
        # # Indicate Floor
        # hboxfloor = QHBoxLayout()
        # floorlabel = QLabel('Floor')
        # floorlabel.setStyleSheet("font: bold 14px;")
        # hboxfloor.addWidget(floorlabel)
        # hboxfloor.setAlignment(Qt.AlignCenter)
        # vbox.addLayout(hboxfloor)



        # validation
        # buttono=QPushButton("Cancel")
        buttonc=QPushButton("Delete")
        # buttono.clicked.connect(self.valide)
        buttonc.clicked.connect(self.cancel)

        hboxDial = QHBoxLayout()
        hboxDial.addWidget(buttonc)
        # hboxDial.addWidget(buttono)



        # create Layout
        vbox.addLayout(hboxDial)
        self.setLayout(vbox)

    def Hline(self):
        """ Create  horizontal line widget
        """
        hline = QFrame()
        hline.setFrameShape(QFrame.HLine)
        hline.setFrameShadow(QFrame.Sunken)
        return hline
    # def _autocalc_height_val(self):
    #     """ split height proportionnaly to the number of subsegs
    #         when new subseg
    #         TO BE DONE
    #     """
    #     self.lma=[]
    #     self.lmi=[]
    #     for s in range(self.Nss):
    #         self.lma.append(self.zmax)
    #         self.lmi.append(self.zmin)
    #     if self.Nss >1:


    #         for s in range(self.Nss):
    #             self.lma.append(((self.Nss-s)*(self.zmax))/self.Nss)
    #             self.lmi.append(((self.Nss-s)*self.zmin)/self.Nss)
    #         self.lma=self.lma[::-1]
    #         self.lmi=self.lmi[::-1]
    #     else:
    #         self.lma.append(self.zmax)
    #         self.lmi.append(self.zmin)
    #     print self.lma
    #     print self.lmi


    def _init_subseg_prop(self):
        
        #TODO
        # sort subsegments from floor to ceil
        # add connect to impose begin of previous segment
        # sort subseg by height

  
  



        self.comboslab = QComboBox()
        for s in self.gparent.L.sl.keys():
            self.comboslab.addItem(s)
        idx=self.comboslab.findText(self.subsegdata['ss_name'])
        self.comboslab.setCurrentIndex(idx)


        self.heightmin = QDoubleSpinBox()
        self.connect(self.heightmin, SIGNAL('valueChanged(double)'), self.color_restore)

        self.heightmin.setObjectName("zmin")
        self.heightmin.setSingleStep(0.01)
        self.heightmin.setRange(0.,self.gparent.L.maxheight)
        self.heightmin.setValue(self.subsegdata['ss_z'][0])

        self.heightmax = QDoubleSpinBox()
        self.connect(self.heightmax, SIGNAL('valueChanged(double)'), self.color_restore)
        self.heightmax.setSingleStep(0.01)
        self.heightmax.setObjectName("zmax")
        self.heightmax.setRange(0.,self.gparent.L.maxheight)
        self.heightmax.setValue(self.subsegdata['ss_z'][1])
        self.offset = QDoubleSpinBox()
        self.offset.setObjectName("offset")
        self.offset.setSingleStep(0.01)
        self.offset.setRange(-1.,1.)
        self.offset.setValue(self.subsegdata['ss_offset'])

        self.connect(self.heightmin, SIGNAL('valueChanged(double)'), self.parent.force_ss_minmax)
        self.connect(self.heightmax, SIGNAL('valueChanged(double)'), self.parent.force_ss_minmax)
        self.connect(self.heightmin, SIGNAL('valueChanged(double)'), self.parent.force_minmax_previous)
        self.connect(self.heightmax, SIGNAL('valueChanged(double)'), self.parent.force_minmax_previous)

    def color_restore(self, event):
        """ retore black color to SpinBox
        """
        sender = self.sender()
        sender.setStyleSheet("color: black")
        # val = self.lheightmax[-1].value()
        # print val


        # if not self.mulseg:
        #     self.gparent.L.edit_seg(self.gparent.selectl.nsel,self.subsegdata)
        # else:
        #     [self.gparent.L.edit_seg(s,self.subsegdata) for s in self.gparent.selectl.selectseg]
        # self.close()

    def cancel(self):
        XX = self.parent.dwidget.pop(self.idx)
        self.close()


class PropertiesWin(QDialog):    # any super class is okay

    def __init__(self,mulseg=False,parent=None):
        super(PropertiesWin, self).__init__(parent)
        # to imporve here. Probably something to inherit from parent App
        self.parent=parent
        # determine if multiple se gments are selected
        self.mulseg=mulseg
        self.widx = 0
        self._init_subsegs()
        self._init_layout()



    def _init_layout(self):
    # main button
        self.addButton = QPushButton('Add H seg')
        self.addButton.clicked.connect(self.addWidget)

        # scroll area widget contents - layout
        self.scrollLayout = QFormLayout()

        # scroll area widget contents
        self.scrollWidget = QWidget()
        self.scrollWidget.setLayout(self.scrollLayout)

        # scroll area
        self.scrollArea = QScrollArea()
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setWidget(self.scrollWidget)

        # main layout
        self.mainLayout = QVBoxLayout()

        # add all main to the main vLayout
        self.mainLayout.addWidget(self.addButton)


        hboxfloor = QHBoxLayout()
        ceillabel = QLabel('Ceil')
        ceillabel.setStyleSheet("font: bold 14px;")
        ceillabel.setAlignment(Qt.AlignCenter)
        self.mainLayout.addWidget(ceillabel)

        self.mainLayout.addWidget(self.scrollArea)

        hboxfloor = QHBoxLayout()
        floorlabel = QLabel('Floor')
        floorlabel.setStyleSheet("font: bold 14px;")
        floorlabel.setAlignment(Qt.AlignCenter)
        self.mainLayout.addWidget(floorlabel)



        buttono=QPushButton("OK")
        buttonc=QPushButton("Cancel")
        buttono.clicked.connect(self.valide)
        buttonc.clicked.connect(self.cancel)

        hboxDial = QHBoxLayout()
        hboxDial.addWidget(buttonc)
        hboxDial.addWidget(buttono)


        # create Layout
        vbox = QVBoxLayout()

        vbox.addLayout(self.mainLayout)
        vbox.addLayout(hboxDial)

        self.setLayout(vbox)


        self.Nss = len(self.subsegdata['ss_name'])

        # order subseg in order of zmin to zmax
        z = np.array(self.subsegdata['ss_z'])

        self.sszo = np.argsort(z[:,0])[::-1]

        self.dwidget = {}

        for ss in self.sszo:
            subsegdata = {'ss_name':self.subsegdata['ss_name'][ss],
                          'ss_z':self.subsegdata['ss_z'][ss],
                          'ss_offset':self.subsegdata['ss_offset'][ss],
                          'ss_ID':self.subsegdata['ss_ID'][ss]}
            self.addWidget(subsegdata)

    def addWidget(self,subsegdata={}):
        if subsegdata==False:
            subsegdata = {'ss_name':'AIR',
                         'ss_z':(0,3.0),
                         'ss_offset':0.,
                         'ss_ID':-1}
        self.dwidget.update({self.widx:SubSegWin(parent=self,idx = self.widx, subsegdata=subsegdata)})
        self.scrollLayout.addRow(self.dwidget[self.widx])
        self.widx = self.widx +1




    def force_ss_minmax(self,event):
        """ Force sub-segment zmin < zmax
        """
        sender = self.sender()
        wid = None
        for k,w in self.dwidget.items():
            if w.heightmin == sender:
                wid = self.dwidget[k]
                break
            elif w.heightmax == sender:
                wid = self.dwidget[k]
                break

        
        zmin =wid.heightmin.value()
        zmax =wid.heightmax.value()
        if zmin >= zmax:
            wid.heightmin.setValue(zmax)

    def force_minmax_previous(self,event):
        sender = self.sender()
        wid = None
        widkey = None
        for k,w in self.dwidget.items():
            if w.heightmin == sender:
                wid = self.dwidget[k]
                widkey = k
                break
            elif w.heightmax == sender:
                wid = self.dwidget[k]
                widkey = k
                break
        # here it is supposed that larger key number are at the lowest level
        lw = self.dwidget.keys()
        lw.sort()
        ulw = lw.index(k)


        # if sender is a zmin
        if 'min' in sender.objectName():
            # not last line <=> bottom of the segment
            if ulw != lw[-1]:
                #next widget
                nwid = self.dwidget[lw[ulw+1]]
                nwid.heightmax.setValue(wid.heightmin.value())


        elif 'max' in sender.objectName():
            # not 1st line <=> top of the segment
            if ulw != lw[0]:
                #previous widget
                pwid = self.dwidget[lw[ulw-1]]
                pwid.heightmin.setValue(wid.heightmax.value())



    def _init_subsegs(self):
        self.subsegdata={}
        self.subsegdata['ss_name']=[]
        self.subsegdata['ss_offset']=[]
        self.subsegdata['ss_z']=[]
        self.subsegdata['ss_ID']=[]
        sub = self.parent.L.Gs.node[self.parent.selectl.nsel]
        Nss = sub['iso']

        # if len(Nss) == 1 => this is the segment and not subseg
        if len(Nss) != 0:
            for uss,ss in enumerate(Nss):
                node = self.parent.L.Gs.node[ss]
                self.subsegdata['ss_name'].append(node['name'])
                self.subsegdata['ss_offset'].append(node['offset'])
                self.subsegdata['ss_z'].append(node['z'])
                self.subsegdata['ss_ID'].append(ss)
        else:
            self.subsegdata['ss_name'].append(sub['name'])
            self.subsegdata['ss_offset'].append(sub['offset'])
            self.subsegdata['ss_z'].append(sub['z'])
            self.subsegdata['ss_ID'].append(self.parent.selectl.nsel)

    def valide(self):
        """ ok click
        """


        node = self.parent.L.Gs.node[self.parent.selectl.nsel]
        n1,n2 = node['connect']
        iso = node['iso']

        newID = [self.dwidget[x].subsegdata['ss_ID'] for x in self.dwidget]


        # import ipdb
        # ipdb.set_trace()
        # no segments just quit discarding changes
        if len(self.dwidget) == 0:
            self.close()

        # diffseg = len(self.dwidget) - len(iso) 
        # # try to reaffect segment number to existing ones
        # if diffseg == 0:
        #     # same number of created seg than existing before edit
        #     segID = iso
        # elif diffseg > 0 :
        #     # segments have been added
        #     segID = iso + [-1]*diffseg
        # elif diffseg < 0:
        #     # segments have been removed
        #     pass
        for w in self.dwidget.values():
            # get nodes
            zmin = w.heightmin.value()
            zmax = w.heightmax.value()
            offset = w.offset.value()
            name = str(w.comboslab.currentText())
            self.parent.L.add_segment(n1,
                                        n2,
                                        num=w.subsegdata['ss_ID'],
                                        maxnum=-1,
                                        transition = False,
                                        name=name, 
                                        z=(zmin, zmax), 
                                        offset=offset,
                                        verbose=False)
            newID.append(w.subsegdata['ss_ID'])
        # determine if a segment has been removed, and then removeit from Gs
        # seg to be removed
        bsegtbr = [i not in newID for i in iso]
        utbr = np.where(bsegtbr)[0]
        if len(utbr) > 0 :
            segtbr = np.array(iso)[utbr].tolist()
            self.parent.L.del_segment(segtbr,g2npy=False)
        self.parent.L.g2npy()


        self.close()

    def cancel(self):
        """ cancel click
        """
        self.close()




class SaveQuitWin(QDialog):    # any super class is okay
    def __init__(self,exit=False,parent=None):
        super(SaveQuitWin, self).__init__(parent)
        self.setWindowTitle('Do you want to save?')
        self.parent=parent
        self.exit=exit

        if self.exit :
            buttonq=QPushButton("Quit Editor")
        else:
            buttonq=QPushButton("Close Layout")
        buttons=QPushButton("Save")
        buttonc=QPushButton("Cancel")
        buttonq.clicked.connect(self.quit)
        buttons.clicked.connect(self.save)
        buttonc.clicked.connect(self.cancel)


        hboxDial = QHBoxLayout()
        hboxDial.addWidget(buttonq)
        hboxDial.addWidget(buttons)
        hboxDial.addWidget(buttonc)


        # create Layout

        self.setLayout(hboxDial)
        print exit
    def quit(self):
        try:
            self.parent.fig.clear()
            self.parent.fig.canvas.draw()
            del self.parent.L
            del self.parent.main_frame

            del self.parent.main_frame
        except:
            pass
        if self.exit :
            self.close()
            self.parent.exitl()
        self.close()

    def save(self):
        self.parent.save()
        if self.exit:
            self.close()
            self.parent.exitl()
        self.close()

    def cancel(self):
        self.close()


class NewLayout(QDialog):    # any super class is okay
    def __init__(self,parent=None,overlay={}):
        super(NewLayout, self).__init__(parent)
        self.parent=parent
        self.doverlay=overlay
        if self.doverlay == {}:
            self.setWindowTitle('New Layout')
            self._init_choices()
            self._init_layoutwin()
        else : 
            self.new()





    def _init_choices(self):
        self.width = QSpinBox()
        self.width.setObjectName("width [m]")
        self.width.setRange(1, 10000)
        self.width.setValue(10)

        self.height = QSpinBox()
        self.height.setObjectName("height [m]")
        self.height.setRange(1, 10000)
        self.height.setValue(10)

    def _init_layoutwin(self):

        vbox = QVBoxLayout()

        # Indicate Ceil
        hboxlabel = QHBoxLayout()
        height = QLabel('Height')
        width = QLabel('Width')
        # ceillabel.setStyleSheet("font: bold 14px;")
        hboxlabel.addWidget(height)
        # hboxlabel.setAlignment(Qt.AlignCenter)
        hboxlabel.addWidget(width)
        # hboxlabel.setAlignment(Qt.AlignCenter)

        vbox.addLayout(hboxlabel)

        hboxlw = QHBoxLayout()
        hboxlw.addWidget(self.height)
        hboxlw.addWidget(self.width)


        vbox.addLayout(hboxlw)

        # validation
        buttonn=QPushButton("New")
        buttonc=QPushButton("Cancel")
        buttonn.setAutoDefault(True)
        buttonn.setDefault(True)
        buttonn.clicked.connect(self.new)
        buttonc.clicked.connect(self.cancel)


        hboxDial = QHBoxLayout()
        hboxDial.addWidget(buttonc)
        hboxDial.addWidget(buttonn)

        vbox.addLayout(hboxDial)

        # create Layout

        self.setLayout(vbox)


    def new(self):
        self.parent.L=Layout()

        if self.doverlay.has_key('overlay_file'):
            self.parent.L.display['overlay_file']=self.doverlay['overlay_file']
            flip =''
            if self.doverlay['flipv']:
                flip = flip + 'v'
            if self.doverlay['fliph']:
                flip = flip + 'h'
            axis = self.doverlay['ax'].axis()
            ax = self.doverlay['ratiox']
            ay = self.doverlay['ratioy']
            dx = self.doverlay['origin'][0]*ax
            dy = self.doverlay['origin'][1]*ay
            axis = (-dx,(axis[1]*ax)-dx,-dy,(axis[3]*ay)-dy)
            self.parent.L.display['overlay_axis'] = axis
            self.parent.L.boundary(xlim=self.parent.L.display['overlay_axis'])
            self.parent.L.display['overlay_flip']=flip
            self.parent.L.display['overlay']=True

        else : 
            self.parent.L.display['overlay_file']=''
            self.parent.L.display['overlay_flip']=''
            lim = (0., self.width.value(), 0.,self.height.value())
            self.parent.L.boundary(xlim=lim)
            self.parent.L.display['overlay']=False

        self.parent.filename=''
        self.parent.create_main_frame()
        self.parent.on_draw()
        self.parent.setWindowTitle(self.parent.L._filename + '- Pylayers : Stand Alone Editor (Alpha)')
        self.parent.resize(self.parent.fig.canvas.width(),self.parent.fig.canvas.height())
        self.close()



    def cancel(self):
        self.close()

class Overset(QMainWindow):
    def __init__(self,parent=None):
        super(Overset, self).__init__(parent)
        self.setWindowTitle('Set Overlay')
        self.x0=np.array([0,0])
        self.x0selected=False
        self.xselected=False
        self.xselected=False
        self.click=False
        self.parent=parent
        self.toolbar()
        self.openoverlay()
        self.showfig()
        self.flipv=False
        self.fliph=False

    def openoverlay(self):
        filename = QFileDialog.getOpenFileName(self,'Open Layout Overlay',pyu.getlong('',pstruc['DIRIMAGE']),'(*.png);;(*.jpg);;(*.jpeg)')

        if filename != '':
            self._fileoverlay = pyu.getshort(str(filename))
            print 'overlay loaded'



    def toolbar(self):
        ###############################
        ### Toolbar
        ###############################
        # origin
        self.xx0 = QDoubleSpinBox()
        self.xx0.setObjectName("x0 [m]")
        self.xx0.setSingleStep(0.01)
        self.xx0.setRange(-1000., 1000.)
        self.xx0.setValue(0)
        self.connect(self.xx0, SIGNAL('valueChanged(double)'), self.refresh)


        self.yy0 = QDoubleSpinBox()
        self.yy0.setObjectName("y0 [m]")
        self.yy0.setSingleStep(0.01)
        self.yy0.setRange(-1000., 1000.)
        self.yy0.setValue(0)
        self.connect(self.yy0, SIGNAL('valueChanged(double)'), self.refresh)


        self.da = QDoubleSpinBox()
        self.da.setObjectName("d_a [m]")
        self.da.setSingleStep(0.01)
        self.da.setRange(1, 10000.)
        self.da.setValue(10)


        self.db = QDoubleSpinBox()
        self.db.setObjectName("d_b [m]")
        self.db.setSingleStep(0.01)
        self.db.setRange(1, 10000.)
        self.db.setValue(10)
        self.connect(self.da, SIGNAL('valueChanged(double)'), self.refresh)

        vbox = QVBoxLayout()

        # Indicate Ceil
        hboxlabel = QHBoxLayout()
        x0 = QLabel('origin x')
        y0 = QLabel('origin y')
        da = QLabel('x')
        db = QLabel('y')
        self.connect(self.db, SIGNAL('valueChanged(double)'), self.refresh)



        self.toolbar0 = QToolBar(self)
        self.toolbar0.addWidget(x0)
        self.toolbar0.addWidget(self.xx0)
        self.toolbar0.addWidget(y0)
        self.toolbar0.addWidget(self.yy0)
        self.toolbara = QToolBar(self)

        self.toolbara.addWidget(da)
        self.toolbara.addWidget(self.da)
        self.toolbarb = QToolBar(self)

        self.toolbarb.addWidget(db)
        self.toolbarb.addWidget(self.db)
        
        # validation

        self.toolbarval = QToolBar(self)

        # Indicate Ceil
        buttonn=QPushButton("Start Editing")
        buttonc=QPushButton("Cancel")
        buttonn.setAutoDefault(True)
        buttonn.setDefault(True)
        buttonn.clicked.connect(self.new)
        buttonc.clicked.connect(self.cancel)

        self.toolbarval.addWidget(buttonn)
        self.toolbarval.addWidget(buttonc)


    def showfig(self):
        self.main_frame = QWidget()
        self.figure = plt.figure()
        self.ax = self.figure.add_subplot(111)
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(self.main_frame)

        self.kpress = self.figure.canvas.mpl_connect('key_press_event', self.flip)

        self.mpress = self.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.mrelea = self.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.kmove = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)


        self.canvas.setFocusPolicy( Qt.ClickFocus )
        self.canvas.setFocus()


        self.navtoolbar = NavigationToolbar(self.canvas, self.main_frame)
        self.image = Image.open(os.path.join(basename,pstruc['DIRIMAGE'],self._fileoverlay))
        self.overlay = self.ax.imshow(self.image,origin='lower')
        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.navtoolbar)
        layout.addWidget(self.canvas)

        # layout.addWidget(self.toolbar)
        self.addToolBar(Qt.ToolBarArea(Qt.BottomToolBarArea), self.toolbar0)
        self.addToolBarBreak()
        self.addToolBar(Qt.ToolBarArea(Qt.BottomToolBarArea), self.toolbara)
        self.addToolBarBreak()
        self.addToolBar(Qt.ToolBarArea(Qt.BottomToolBarArea), self.toolbarb)
        self.addToolBarBreak()
        self.addToolBar(Qt.ToolBarArea(Qt.BottomToolBarArea), self.toolbarval)
        # layout.addWidget(self.button)
        self.main_frame.setLayout(layout)
        self.setCentralWidget(self.main_frame)

        self.axis = self.ax.axis()
        self.x = self.axis[1]*0.1
        self.y = self.axis[3]*0.1


        self.refresh()
        self.canvas.draw()


    def flip(self,event):
        if (event.key == 'v') or (event.key == 'V'):
            self.image = self.image.transpose(Image.FLIP_LEFT_RIGHT)
            self.flipv = ~self.flipv
        if (event.key == 'h') or (event.key == 'H'):
            self.image = self.image.transpose(Image.FLIP_TOP_BOTTOM )
            self.fliph = ~self.fliph
        if (event.key == 'v') or (event.key == 'V') or (event.key == 'h') or (event.key == 'H'):
            self.overlay.remove()
            self.overlay = self.ax.imshow(self.image,origin='lower')
            self.refresh()


    def on_press(self,event):
        mouse  = np.array([event.xdata,event.ydata])

        dx = abs(self.x0[0]-self.x)
        dy = abs(self.x0[1]-self.y)
        if event.button == 1:
            self.click=True
            if np.sqrt((mouse[0]-self.x0[0])**2+(mouse[1]-self.x0[1])**2)<10:
                self.x0selected=True
            elif np.sqrt((mouse[0]-self.x)**2+(mouse[1]-self.x0[1])**2)<0.3*dx:
                self.xselected=True
            elif np.sqrt((mouse[0]-self.x0[0])**2+(mouse[1]-self.y)**2)<0.3*dy:
                self.yselected=True


    def on_motion(self,event):
        if self.click:
            if self.x0selected:
                self.xx0.setValue(event.xdata)
                self.yy0.setValue(event.ydata)
            if self.xselected:
                self.x = event.xdata
            if self.yselected:
                self.y = event.ydata
        self.refresh()


    def on_release(self,event):
        self.click=False
        self.x0selected=False
        self.xselected=False
        self.yselected=False
        self.refresh()


    def refresh(self):

        ptrm = ['p0','p00','p1','p2','a0','a1','b0','b1']
        for pt in ptrm:
            try:
                eval('self.'+pt+'.remove()')
            except: 
                pass
        
        self.x0  = np.array([self.xx0.value(),self.yy0.value()])

        self.p00 = self.ax.scatter(self.x0[0],self.x0[1],marker='o',s=300,linewidths=8, facecolors='None', edgecolors='r')
        self.p0 = self.ax.scatter(self.x0[0],self.x0[1],marker='+',c='k',s=100,linewidths=5)
        self.p1 = self.ax.scatter(self.x,self.x0[1],marker='x',c='r',s=100,linewidths=5)
        self.p2 = self.ax.scatter(self.x0[0],self.y,marker='x',c='b',s=100,linewidths=5)

        mx = (self.x0[0] + self.x)/2.
        my = (self.x0[1] + self.y)/2.
        #da
        self.a0 = self.ax.annotate('', xy=(self.x0[0], self.x0[1]), xycoords='data',xytext=(self.x, self.x0[1]), textcoords='data',arrowprops={'arrowstyle': '<->'})
        # self.a1 = self.ax.annotate('a'+str(self.da.value()) + ' m', xy=(mx, self.x0[1]-100), xycoords='data',xytext=(5, 0), textcoords='offset points')
        self.a1 = self.ax.text(mx, self.x0[1]-50, 'x='+str(self.da.value()) + ' m', fontsize=15)
        #db
        self.b0 =self.ax.annotate('', xy=(self.x0[0], self.x0[1]), xycoords='data',xytext=(self.x0[0], self.y), textcoords='data',arrowprops={'arrowstyle': '<->'})
        # self.b1 =self.ax.annotate('b'+str(self.db.value()) + ' m', xy=(self.x0[0]-100,my), xycoords='data',xytext=(5, 0), textcoords='offset points')
        self.b1 = self.ax.text(self.x0[0]-50,my, 'y='+str(self.db.value()) + ' m', fontsize=15,rotation=90)
        self.canvas.draw()


    def compute_values(self):
        self.ratiox = self.da.value()/abs((self.x0[0]-self.x))
        self.ratioy = self.db.value()/abs((self.x0[1]-self.y))


    def cancel(self):
        self.close()

    def new(self):
        self.figure.canvas.mpl_disconnect(self.kpress)
        self.figure.canvas.mpl_disconnect(self.mpress)
        self.figure.canvas.mpl_disconnect(self.mrelea)
        self.figure.canvas.mpl_disconnect(self.kmove)

        
        self.compute_values()
        doverlay={'flipv':self.flipv,
                  'fliph':self.fliph,
                  'overlay_file':self._fileoverlay,
                  'origin':self.x0,
                  'ratiox':self.ratiox,
                  'ratioy':self.ratioy,
                  'ax':self.ax}
        self.parent.newlayout=NewLayout(parent=self.parent,overlay=doverlay)
        self.parent.newlayout.show()
        self.close()


class GridSet(QDialog):    # any super class is okay
    def __init__(self,parent=None):
        super(GridSet, self).__init__(parent)
        self.setWindowTitle('Set Grid')
        self.parent=parent
        self._init_choices()
        self._init_layoutwin()



    def _init_choices(self):
        self.xspacing = QDoubleSpinBox()
        self.xspacing.setObjectName("x spacing [m]")
        self.xspacing.setRange(0, 100)
        self.xspacing.setValue(self.parent.selectl.stepgridx)

        self.yspacing = QDoubleSpinBox()
        self.yspacing.setObjectName("y spacing [m]")
        self.yspacing.setRange(0, 100)
        self.yspacing.setValue(self.parent.selectl.stepgridy)

    def _init_layoutwin(self):

        vbox = QVBoxLayout()

        # Indicate Ceil
        hboxlabel = QHBoxLayout()
        yspacing = QLabel('y spacing')
        xspacing = QLabel('x spacing')
        # ceillabel.setStyleSheet("font: bold 14px;"
        hboxlabel.addWidget(xspacing)
        hboxlabel.addWidget(yspacing)
        # hboxlabel.setAlignment(Qt.AlignCenter)
        # hboxlabel.setAlignment(Qt.AlignCenter)

        vbox.addLayout(hboxlabel)

        hboxlw = QHBoxLayout()
        hboxlw.addWidget(self.xspacing)
        hboxlw.addWidget(self.yspacing)


        vbox.addLayout(hboxlw)

        # validation
        buttonn=QPushButton("Ok")
        buttonc=QPushButton("Cancel")
        buttonn.clicked.connect(self.ok)
        buttonc.clicked.connect(self.cancel)


        hboxDial = QHBoxLayout()
        hboxDial.addWidget(buttonc)
        hboxDial.addWidget(buttonn)

        vbox.addLayout(hboxDial)

        # create Layout

        self.setLayout(vbox)


    def ok(self):
        print self.xspacing.value(),self.yspacing.value()
        self.parent.selectl.stepgridx = self.xspacing.value()
        self.parent.selectl.stepgridy = self.yspacing.value()
        self.parent.selectl.gridOn=True
        if not self.parent.selectl.gridOn:
            self.parent.selectl.setgrid()
            self.parent.fig.canvas.draw()
        self.close()


    def cancel(self):
        self.close()



class AppForm(QMainWindow):
    def __init__(self, parent=None):
        super(AppForm,self).__init__()
        # QMainWindow.__init__(self, parent)
        self.setWindowTitle('Pylayers : Stand Alone Editor (Alpha)')
        self.filename=''

        self.create_menu()
        self.create_status_bar()
        self.shortcuts()
        if 'darwin' in sys.platform:
            self.create_toolbar()

        self.show3On = False
        # self.__loaddebug()



    def __loaddebug(self):
        self.L = Layout('defstr.lay')
        self.filename = self.L._filename
        self.create_main_frame()
        self.on_draw()
        self.setWindowTitle(self.L._filename + '- Pylayers : Stand Alone Editor (Alpha)')
        self.resize(self.fig.canvas.width(),self.fig.canvas.height())
        print 'loaded'


    def new(self):
        self.newlayout=NewLayout(parent=self)
        self.newlayout.show()


    def newover(self):
        self.overset = Overset(parent=self)
        self.overset.show()

    def open(self):
        filename = QFileDialog.getOpenFileName(self,'Open Pylayers Layout File',pyu.getlong('',pstruc['DIRLAY']),'(*.lay);;(*.ini)')

        if filename != '':
            _filename = pyu.getshort(str(filename))
            self.L = Layout(_filename)
            self.filename = self.L._filename
            self.create_main_frame()
            self.on_draw()
            self.setWindowTitle(self.L._filename + '- Pylayers : Stand Alone Editor (Alpha)')
            self.resize(self.fig.canvas.width(),self.fig.canvas.height())
            print 'loaded'

        # self.setgrid()
    def save(self,force=False):

        if self.filename == '' or force:
            filename = QFileDialog.getSaveFileName(self, 'Save Layout', pyu.getlong('',pstruc['DIRINI']),'(*.ini);;(*.osm)')
            try:
                _filename= pyu.getshort(str(filename))
            except:
                pass
        else :
            _filename=self.L._filename
        try:
            oldCursor = QCursor()
            QApplication.setOverrideCursor(QCursor(Qt.BusyCursor))
            self.L.saveini(_filename)
            self.L.saveosm(_filename.split('.')[0] + '.osm')
            # self.L = Layout(_filename)
            self.filename=self.L._filename
            self.setWindowTitle(self.L._filename + '- Pylayers : Stand Alone Editor (Alpha)')
            QApplication.setOverrideCursor(oldCursor)

            print 'saved'
        except:
            pass

    def chgover(self):
        """ change overlay file
        """
        filename = QFileDialog.getOpenFileName(self,'Open Pylayers Layout File',pyu.getlong('',pstruc['DIRIMAGE']),'(*.png);;(*.jpg)')

        if filename != '':
            _filename= pyu.getshort(str(filename))
            self.L.display['fileoverlay']=_filename
            print 'loaded overlay'

    def toggleoverl(self):
        """ toggle overlay
        """
        self.selectl.toggleoverlay()
        self.fig.canvas.draw()

    def closel(self,exit=False):
        dial_res=''
        self.sq = SaveQuitWin(parent=self,exit=exit)
        self.sq.show()

    def exitl(self):
        try:
            plt.rcParams.update(self.selectl.rcconf)
            self.selectl.fig.canvas.mpl_disconnect(self.cid1)
            self.selectl.fig.canvas.mpl_disconnect(self.cid2)
            self.selectl.fig.canvas.mpl_disconnect(self.cid3)
            self.selectl.fig.canvas.mpl_disconnect(self.cid4)
            self.selectl.fig.canvas.mpl_disconnect(self.cid5)
            self.selectl.fig.canvas.mpl_disconnect(self.cid6)
        except:
            pass        

        QApplication.quit()

    def edit_properties(self):
        """ edit wall properties
        """


        if (self.selectl.state == 'SS') and (self.selectl.nsel > 0):
            self.prop = PropertiesWin(parent=self,mulseg=False)
            self.prop.show()
        elif (self.selectl.state == 'SMS') and (self.selectl.selectseg!=[]):
            self.prop = PropertiesWin(parent=self,mulseg=True)
            self.prop.show()
        elif (self.selectl.state == 'SMP') and (self.selectl.selectseg!=[]):
            self.selectl.toggle()
            self.prop = PropertiesWin(parent=self,mulseg=True)
            self.prop.show()

        # self.on_draw()

    def editgrid(self):
        grid = GridSet(parent=self)
        grid.show()

    def togglegrid(self):
        self.selectl.togglegrid()

    def snapongrid(self):
        self.selectl.toggglesnapgrid()

    def toggleshow3(self):
        if not self.show3On:
            self.show3On = True
            self.show3()
        elif self.show3On:
            mlab.close()
            self.show3On = False

    def show3(self):
        if self.show3On:
            mlab.clf()
            self.L._show3()

    def updatelayerselector(self):
        slname={}
        slname['name']=str(self.layerselector.currentText())
        if self.selectl.state == 'Init' or self.selectl.state == 'SS':
            if self.selectl.nsel > 0:
                if (self.selectl.state == 'SS'):
                    self.L.edit_seg(self.selectl.nsel,slname)
        elif self.selectl.state == 'SMS' or self.selectl.state == 'SMP':
            [self.L.edit_seg(sl,slname) for sl in self.selectl.selectseg]



    def selectnodes(self):
        ''' select mode, managed by selectl
            here only cursor management
        '''
        QApplication.setOverrideCursor(QCursor(Qt.ArrowCursor))
        self.selectl.escape()
        string = self.selectl.help[self.selectl.state]
        self.statusBar().showMessage(string)


    def drawseg(self):
        ''' drawseg, managed by selectl
            here only cursor management
        '''

        QApplication.setOverrideCursor(QCursor(Qt.CrossCursor))
        self.L.display['activelayer']=str(self.layerselector.currentText())
        self.selectl.current_layer=self.L.display['activelayer']
        self.selectl.modeCP()
        string = self.selectl.help[self.selectl.state]
        self.statusBar().showMessage(string)




    def on_about(self):
        msg = """ This is the PyLayers' Stand-Alone Layout Editor (Alpha)

         This tool allows to edit/modyfy a building floor plan and/or constitutive materials.
         Once saved, the layout s ready to be used with PyLayers simunlation tools.



         Shortcuts:
         ------------

         F1 : Select mode
         F2 : New segment with current active Layer
         F3 : Edit segment properties

         g : toggle grid
         ctrl+g : choose grid properties

         CTRL + o : Open Layout
         CTRL + s : Save Layout
         CTRL + q : Quit Editor
         escape : back to a stable state

         More hints about editing can be found in the status bar.



         Thank you for using Pylayers and this tool

         The Pylayers' Dev Team
         www.pylayers.org

        """
        QMessageBox.about(self, "Pylayers' Stand-Alone Layout Editor (Alpha)", msg.strip())





    def on_draw(self):
        """ Redraws the figure
        """
        # str = unicode(self.textbox.text())
        # self.data = map(int, str.split())

        # x = range(len(self.data))

        # clear the axes and redraw the plot anew
        #
        # self.axes.clear()

        # self.axes.grid(self.grid_cb.isChecked())
        self.L.display['nodes']=True
        self.L.display['ednodes']=True
        self.L.display['subseg']=False
        self.L.display['isonb']=True
        self.L.display['ticksoff']=False


        self.fig,self.axes = self.selectl.show(self.fig,self.axes,clear=True)
        # self.axes.text(10,10,str(self.properties.currentText()))

        # self.L.showGs(fig=self.fig,ax=self.axes)
        # self.axes.bar(
        #     left=x,
        #     height=self.data,
        #     width=self.slider.value() / 100.0,
        #     align='center',
        #     alpha=0.44,
        #     picker=5)

        self.fig.canvas.draw()

    def on_release(self,event):
        string=''
        try:
            string = string + ' ' + self.L.Gs.node[self.selectl.nsel]['name']
        except:
            pass
        try:
            string = string + ' with ' +str(len(self.L.Gs.node[self.selectl.nsel]['ss_name'])) + 'subseg(s)'
        except:
            pass
        try:
            n1,n2 = self.L.Gs[self.selectl.nsel].keys()
            pn1 = np.array(self.L.Gs.pos[n1])
            pn2 = np.array(self.L.Gs.pos[n2])
            l="%.2f"%np.sqrt(np.sum((pn1-pn2)**2))

            string = string + '     length= ' + l + 'm     '
        except:
            pass
        string = string +'\t'+self.selectl.help[self.selectl.state]
        self.statusBar().showMessage(string)

        if self.selectl.nsel > 0 :
            try:
                idx=self.layerselector.findText(self.L.Gs.node[self.selectl.nsel]['name'])
                self.layerselector.setCurrentIndex(idx)
            except:
                pass

        if self.show3On:
            self.show3()


    def create_main_frame(self):

        self.main_frame = QWidget()
        self.create_toolbar()
        self.addToolBar(Qt.ToolBarArea(Qt.TopToolBarArea), self.toolbar)


        # Create the mpl Figure and FigCanvas objects.
        # 5x4 inches, 100 dots-per-inch
        #

        self.dpi = 100
        self.fig = Figure((20.0, 30.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)

        # Since we have only one plot, we can use add_axes
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)
        try:
            self.axes.axis(self.L.display['overlay_axis'])
        except:
            self.axes.axis(self.L.display['box'])

        # Bind the 'pick' event for clicking on one of the bars
        #
        # self.canvas.mpl_connect('pick_event', self.on_pick)
        self.selectl = SelectL2(self.L,fig=self.fig,ax=self.axes)

        self.cid1 = self.canvas.mpl_connect('button_press_event',
                                           self.selectl.OnClick)
        self.cid2 = self.canvas.mpl_connect('button_release_event',
                                           self.selectl.OnClickRelease)
        self.cid3 = self.canvas.mpl_connect('motion_notify_event',
                                           self.selectl.OnMotion)
        self.cid4 = self.canvas.mpl_connect('key_press_event',
                                           self.selectl.OnPress)
        self.cid5 = self.canvas.mpl_connect('key_release_event',
                                           self.selectl.OnRelease)
        self.cid6 = self.canvas.mpl_connect('button_release_event',
                                           self.on_release)
        self.canvas.setFocusPolicy( Qt.ClickFocus )
        self.canvas.setFocus()



        #Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

        vbox = QVBoxLayout()
        # vbox.addLayout(layerbox)
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)

        self.main_frame.setLayout(vbox)

        self.setCentralWidget(self.main_frame)



    def create_status_bar(self):
        self.status_text = QLabel("Open a Layout")
        self.statusBar().addWidget(self.status_text, 1)


    def shortcuts(self):
        esc = QShortcut(self)
        esc.setKey("escape")
        self.connect(esc, SIGNAL("activated()"), self.selectnodes)

    def refresh(self):
        f5 = QShortcut(self)
        f5.setKey("F5")
        self.connect(f5, SIGNAL("activated()"), self.selectl.refresh)


    def create_menu(self):
        self.file_menu = self.menuBar().addMenu("&File")
        self.edit_menu = self.menuBar().addMenu("&Edit")
        self.view_menu = self.menuBar().addMenu("&View")

        self.help_menu = self.menuBar().addMenu("&Help")
        # load_file_action = self.create_action("&Save plot",
        #     shortcut="Ctrl+S", slot=self.save_plot,
        #     tip="Save the plot")
        new_action = self.create_action("&New Layout", slot=self.new,
        shortcut="Ctrl+n", tip="new layout")
        new_overlay = self.create_action("&New Overlay", slot=self.newover,
        shortcut="Ctrl+Shift+n", tip="new layout from overlay")
        open_action = self.create_action("&Open", slot=self.open,
        shortcut="Ctrl+o", tip="Open Layout")
        save_action = self.create_action("&Save", slot=self.save,
        shortcut="Ctrl+s", tip="Save Layout")
        saveas_action = self.create_action("&Save as...", slot=lambda x=True:self.save(x),
        shortcut="Ctrl+Shift+s", tip="Save as")
        # open_action = self.create_action("&Open", slot=self.open,
        # shortcut="Ctrl+o", tip="Open Layout")
        close_action = self.create_action("&Close", shortcut='Ctrl+w', slot=self.closel, tip="Close Layout")
        quit_action = self.create_action("&Quit", slot=lambda x=True:self.closel(x),
            shortcut="Ctrl+Q", tip="Close the application")

        select_action = self.create_action("&Select Nodes", slot=self.selectnodes,
            shortcut="F1", tip="Select Nodes")
        draw_action = self.create_action("&Draw Segments", slot=self.drawseg,
            shortcut="F2", tip="Draw segements")


        refresh = self.create_action("&Refresh", slot=self.on_draw,
            shortcut="F10", tip="Refresh the application")
        properties= self.create_action("&Properties", slot=self.edit_properties,
            shortcut="F3", tip="Edit Wall properties")
        # show3= self.create_action("&Properties", slot=self.edit_properties,
        #     shortcut="F9", tip="3D show")

        about_action = self.create_action("&About",
            shortcut='F12', slot=self.on_about,
            tip='about')

        gridset_action = self.create_action("&Grid",
            shortcut='', slot=self.editgrid,
            tip='Set Grid',)
        snapongrid_action = self.create_action("&Snap On Grid",
            shortcut='s', slot=self.snapongrid,
            tip='Snap on Grid',checkable=True)

        gridtg_action = self.create_action("&Toggle Grid",
            shortcut='g', slot=self.togglegrid,
            tip='toggle Grid',checkable=True)

        view3D_action = self.create_action("&3D View",
            shortcut='3', slot=self.toggleshow3,
            tip='Display 3D view',checkable=True)

        chgoverlay_action = self.create_action("&Choose overlay", slot=self.chgover,
        shortcut="", tip="Choose ovelay")
        toggleover_action = self.create_action("&Toggle overlay", slot=self.toggleoverl,
        shortcut="", tip="Toggle ovelay display")

        self.add_actions(self.file_menu,
            ( new_action,new_overlay,open_action,None,save_action,saveas_action,None,close_action,quit_action,))

        self.add_actions(self.edit_menu,
            ( select_action,draw_action,properties,None,gridset_action,snapongrid_action,gridtg_action,None,refresh))

        self.add_actions(self.view_menu, (view3D_action,chgoverlay_action,toggleover_action))


        self.add_actions(self.help_menu, (about_action,))





    def create_toolbar(self):
        self.toolbar = QToolBar(self)
        ###############################
        ### Toolbar
        ###############################
        # get icons path
        iconpath = os.path.join(pylayersdir,'pylayers','gui','ico')
        # exit
        exitAction = QAction(QIcon(os.path.join(iconpath,'gnome_application_exit.png')), 'Quit', self)
        # exitAction.triggered.connect(lambda x=True:self.closel(x))
        self.toolbar.addAction(exitAction)

        #new
        newAction = QAction(QIcon(os.path.join(iconpath,'gnome_document_new.png')), 'new', self)
        newAction.triggered.connect(self.new)
        self.toolbar.addAction(newAction)

        #open
        openAction = QAction(QIcon(os.path.join(iconpath,'gnome_folder_open.png')), 'Open', self)
        openAction.triggered.connect(self.open)
        self.toolbar.addAction(openAction)

        #save
        saveAction = QAction(QIcon(os.path.join(iconpath,'gnome_document_save.png')), 'Save', self)
        saveAction.triggered.connect(self.save)
        self.toolbar.addAction(saveAction)


        self.toolbar.addSeparator()

        #select
        selectAction = QAction(QIcon(os.path.join(iconpath,'select.png')), 'Select', self)
        selectAction.triggered.connect(self.selectnodes)
        self.toolbar.addAction(selectAction)

        #draw
        drawAction = QAction(QIcon(os.path.join(iconpath,'gnome_list_add.png')), 'Draw Segments', self)
        drawAction.triggered.connect(self.drawseg)
        self.toolbar.addAction(drawAction)

        #edit
        editAction = QAction(QIcon(os.path.join(iconpath,'gnome_accessories_text_editor.png')), 'Edit Segments', self)
        editAction.triggered.connect(self.edit_properties)
        self.toolbar.addAction(editAction)

        self.toolbar.addSeparator()

        # self.addAction()
        #editgrid
        editgridAction = QAction(QIcon(os.path.join(iconpath,'editgrid.png')), 'Edit Grid', self)
        editgridAction.triggered.connect(self.editgrid)
        self.toolbar.addAction(editgridAction)

        #grid
        gridAction = QAction(QIcon(os.path.join(iconpath,'grid.png')), 'Toggle Grid', self)
        gridAction.triggered.connect(self.togglegrid)
        gridAction.setCheckable(True)
        self.toolbar.addAction(gridAction)

        #snapgrid
        snapgridAction = QAction(QIcon(os.path.join(iconpath,'grid_snap.png')), 'Snap On Grid', self)
        snapgridAction.triggered.connect(self.snapongrid)
        snapgridAction.setCheckable(True)
        self.toolbar.addAction(snapgridAction)


        self.toolbar.addSeparator()
        #show3D
        show3Action = QAction(QIcon(os.path.join(iconpath,'sugar_cube.png')), '3D View', self)
        show3Action.triggered.connect(self.toggleshow3)
        show3Action.setCheckable(True)
        self.toolbar.addAction(show3Action)


        self.toolbar.addSeparator()


        # Active layer Menu in toolbar
        layerbox = QHBoxLayout()

        layerlabel = QLabel('Active Layer')
        layerlabel.setStyleSheet("font: 16px;")
        layerlabel.setAlignment(Qt.AlignCenter)

        self.toolbar.addWidget(layerlabel)

        try:
            self.layerselector=QComboBox()
            for s in self.L.sl.keys():
                self.layerselector.addItem(s)
            self.toolbar.addWidget(self.layerselector)
        except:
            pass
        self.layerselector.activated.connect(self.updatelayerselector)


    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None,
                        icon=None, tip=None, checkable=False,
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)

        return action


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    # form.setGeometry(100,100,300,300)
    form.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
