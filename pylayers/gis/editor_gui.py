# -*- coding: utf-8 -*-
import sys, os, random
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from pylayers.gis.layout import *
from pylayers.gis.editor_select import SelectL2



class SubSegWin(QDialog):    # any super class is okay
    def __init__(self,Nss=1,zmin=0.,zmax=3.0,subsegdata={},parent=None):
        super(SubSegWin, self).__init__(parent)
        #
        self.gparent=parent.parent
        # mulsti segment selection indicator
        self.mulseg=parent.mulseg
        # dictionnary to pass subseg data
        self.subsegdata=subsegdata
        self.Nss=Nss
        self.zmin=zmin
        self.zmax=zmax
        self._autocalc_height_val()
        self._init_subseg_prop()
        self._init_layout()

    def _init_layout(self):

        vbox = QVBoxLayout()
        for ss in range(self.Nss):
            # slab
            hboxtitle=QHBoxLayout()
            hboxtitle.addWidget(QLabel('Sub-Segment'+str(ss+1)))
            hbox1 = QHBoxLayout() 
            hbox1.addWidget(self.lcomboslab[ss])

            # slab prop
            hboxl2 = QHBoxLayout() 
            hbox2 = QHBoxLayout() 
            label=['zmin','zmax','offset','']
            for iw,w in enumerate([ self.lheightmin[ss],self.lheightmax[ss],
                                    self.loffset[ss]]):
                hboxl2.addWidget(QLabel(label[iw]))
                hbox2.addWidget(w)
                
            vbox.addLayout(hboxtitle)
            vbox.addLayout(hbox1)
            vbox.addLayout(hboxl2)
            vbox.addLayout(hbox2)

        # validation
        buttono=QPushButton("OK")
        buttonc=QPushButton("Cancel")
        buttono.clicked.connect(self.valide)
        buttonc.clicked.connect(self.cancel)

        hboxDial = QHBoxLayout() 
        hboxDial.addWidget(buttonc)
        hboxDial.addWidget(buttono)


        # create Layout
        vbox.addLayout(hboxDial)
        self.setLayout(vbox)

    def _autocalc_height_val(self):
        """ split height proportionnaly to the number of subsegs
            TO BE DONE
        """
        self.lma=[]
        self.lmi=[]
        for s in range(self.Nss):
            self.lma.append(self.zmax)
            self.lmi.append(self.zmin)
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
        self.lheightmin=[]
        self.lheightmax=[]
        self.loffset=[]

        self.lcomboslab = [] 

        for ss in range(self.Nss):

            self.lcomboslab.append(QComboBox())
            for s in self.gparent.L.sl.keys():
                self.lcomboslab[ss].addItem(s)
            idx=self.lcomboslab[ss].findText(self.subsegdata['ss_name'][ss])
            self.lcomboslab[ss].setCurrentIndex(idx)


            self.lheightmin.append(QDoubleSpinBox())
            self.lheightmin[-1].setObjectName("zmin")
            self.lheightmin[-1].setSingleStep(0.01)
            self.lheightmin[-1].setRange(0.,self.gparent.L.maxheight)
            self.lheightmin[-1].setValue(self.subsegdata['ss_z'][ss][0])

            self.lheightmax.append(QDoubleSpinBox())
            self.lheightmax[-1].setSingleStep(0.01)
            self.lheightmax[-1].setObjectName("zmax")
            self.lheightmax[-1].setRange(0.,self.gparent.L.maxheight)
            self.lheightmax[-1].setValue(self.subsegdata['ss_z'][ss][1])

            self.loffset.append(QDoubleSpinBox())
            self.loffset[-1].setObjectName("offset")
            self.loffset[-1].setSingleStep(0.01)
            self.loffset[-1].setRange(-1.,1.)
            self.loffset[-1].setValue(self.subsegdata['ss_offset'][ss])


    def valide(self):
        # self.subsegdata={}
        # self.subsegdata['ss_name']=[]
        # self.subsegdata['ss_z']=[]
        # self.subsegdata['ss_offset']=[]
        
        # for ss in range(self.Nss):
        #     z = (self.lheightmin[ss].value(),self.lheightmax[ss].value())
        #     self.subsegdata['ss_name'].append(str(self.lcomboslab[ss].currentText()))
        #     self.subsegdata['ss_z'].append(z)
        #     self.subsegdata['ss_offset'].append(self.loffset[ss].value())
   
        # if not self.mulseg:
        #     self.gparent.L.edit_seg(self.gparent.selectl.nsel,self.subsegdata)
        # else:
        #     [self.gparent.L.edit_seg(s,self.subsegdata) for s in self.gparent.selectl.selectseg]
        self.close()

    def cancel(self):
        self.close()

class PropertiesWin(QDialog):    # any super class is okay
    def __init__(self,mulseg=False,parent=None):
        super(PropertiesWin, self).__init__(parent)
        # to imporve here. Probably something to inherit from parent App
        self.parent=parent
        # determine if multiple segments are selected
        self.mulseg=mulseg

        # combo box 
        self._init_slab_prop()
        self._init_subsegs()
        self._init_layout()
        print self.parent.selectl.selectseg


        # self.button.clicked.connect(self.create_child)
    


    def _init_subsegs(self):
        self.subsegdata={}
        self.subsegdata['ss_name']=[]
        self.subsegdata['ss_offset']=[]
        self.subsegdata['ss_z']=[]

        if self.parent.selectl.nsel in self.parent.L.lsss :
            sub = self.parent.L.Gs.node[self.parent.selectl.nsel]
            Nss = len(sub['ss_name']) 
        else : Nss=0

        for ss in range(Nss):
            self.subsegdata['ss_name'].append(sub['ss_name'][ss])
            self.subsegdata['ss_offset'].append(sub['ss_offset'][ss])
            self.subsegdata['ss_z'].append(sub['ss_z'][ss])

    def _init_slab_prop(self):

        # if mulseg, default value are default
        if self.mulseg:
            self.segdata={}
            self.segdata['name']=self.parent.L.sl.keys()[0]
            self.segdata['offset']=0.0
            self.segdata['z']=(0.0,self.parent.L.maxheight)
            self.segdata['transition']=False
            
        # else : default are read from self.Gs.node[node]
        else : 
            self.segdata = self.parent.L.Gs.node[self.parent.selectl.nsel]


        self.comboslab = QComboBox()
        for s in self.parent.L.sl.keys():
            self.comboslab.addItem(s)
        idx=self.comboslab.findText(self.segdata['name'])
        self.comboslab.setCurrentIndex(idx)



        self.heightmin = QDoubleSpinBox()       
        self.heightmin.setObjectName("zmin")
        self.heightmin.setSingleStep(0.01)
        self.heightmin.setRange(0.,800.)
        self.heightmin.setValue(self.segdata['z'][0])


        self.heightmax = QDoubleSpinBox()
        self.heightmax.setSingleStep(0.01)
        self.heightmax.setObjectName("zmax")
        self.heightmax.setRange(0.,800.)
        self.heightmax.setValue(self.segdata['z'][1])

        self.offset =  QDoubleSpinBox()
        self.offset.setObjectName("offset")
        self.offset.setSingleStep(0.01)
        self.offset.setRange(-1.,1.)
        self.offset.setValue(self.segdata['offset'])

        self.transition = QPushButton("Transition")
        self.transition.setCheckable(True)
        self.transition.setDefault(self.segdata['transition'])
        # self.transition.setText("0")

        self.nbsubseg = QSpinBox()
        self.nbsubseg.setObjectName("nbsubseg")
        self.nbsubseg.setRange(0,20.)
        try:
            self.nbsubseg.setValue(len(self.segdata['ss_name']))
        except:
            self.nbsubseg.setValue(0.)

        self.editssbutton = QPushButton("Edit Sub-Segments")
        self.editssbutton.clicked.connect(self.editsubseg)

        

        self.heightmin.setMinimumWidth(5)
        self.heightmax.setMinimumWidth(5)
        self.transition.setMinimumWidth(5)
        

        self.heightmin.setMaximumWidth(70)
        self.heightmax.setMaximumWidth(70)
        self.transition.setMaximumWidth(70)
        self.nbsubseg.setMaximumWidth(50)
        self.editssbutton.setMaximumWidth(120)

    def _init_layout(self):

        # slab
        hbox1 = QHBoxLayout() 
        hbox1.addWidget(self.comboslab)

        # slab prop
        hbox2 = QHBoxLayout() 
        hboxl2 = QHBoxLayout() 
        label=['zmin','zmax','offset','']
        for iw,w in enumerate([ self.heightmin,self.heightmax,self.offset,self.transition]):
            hboxl2.addWidget(QLabel(label[iw]))
            hbox2.addWidget(w)
            # hbox2.setAlignment(w, Qt.AlignVCenter)

        # subseg prop
        hbox3 = QHBoxLayout() 
        hboxl3 = QHBoxLayout() 

        hbox3.addWidget(self.nbsubseg)
        hbox3.addWidget(self.editssbutton)
        hboxl3.addWidget(QLabel('Number of \nSub-segments'))





        # validation
        buttono=QPushButton("OK")
        buttonc=QPushButton("Cancel")
        buttono.clicked.connect(self.valide)
        buttonc.clicked.connect(self.cancel)

        hboxDial = QHBoxLayout() 
        hboxDial.addWidget(buttonc)
        hboxDial.addWidget(buttono)


        # create Layout
        vbox = QVBoxLayout()
        vbox.addLayout(hbox1)
        vbox.addLayout(hboxl2)
        vbox.addLayout(hbox2)
        vbox.addLayout(hboxl3)
        vbox.addLayout(hbox3)
        vbox.addLayout(hboxDial)

        self.setLayout(vbox)

    def editsubseg(self):
        """ open a edit subseg window
        """

        Nss=self.nbsubseg.value()
        if Nss>0:
            zmin=self.heightmin.value()
            zmax=self.heightmax.value()
            if Nss>len(self.subsegdata['ss_name']):
                for i in range(Nss-len(self.subsegdata['ss_name'])):
                    self.subsegdata['ss_name'].append(self.parent.L.sl.keys()[0])
                    self.subsegdata['ss_offset'].append(0.)
                    self.subsegdata['ss_z'].append((zmin,zmax))

            self.subseg=SubSegWin(parent=self,subsegdata=self.subsegdata,Nss=Nss,zmin=zmin,zmax=zmax)

            self.subseg.show()

    def valide(self):
        """ ok click
        """
        z = (self.heightmin.value(),self.heightmax.value())
        if self.transition.isChecked():
            trans = True
        else :
            trans = False
        self.segdata = {'name':str(self.comboslab.currentText()),
                'z':z,
                'transition':trans,
                'offset':self.offset.value()
                }

        Nss = self.nbsubseg.value()

        if Nss>len(self.subsegdata['ss_name']):
            zmin=self.heightmin.value()
            zmax=self.heightmax.value()
            for i in range(Nss-len(self.subsegdata['ss_name'])):
                self.subsegdata['ss_name'].append(self.parent.L.sl.keys()[0])
                self.subsegdata['ss_offset'].append(0.)
                self.subsegdata['ss_z'].append((zmin,zmax))

        self.segdata.update({'ss_name':self.subsegdata['ss_name'][:Nss],
                             'ss_offset':self.subsegdata['ss_offset'][:Nss],
                             'ss_z':self.subsegdata['ss_z'][:Nss]})
        

        if not self.mulseg:
            self.parent.L.edit_seg(self.parent.selectl.nsel,self.segdata)
            self.parent.selectl.modeIni()
        else:
            [self.parent.L.edit_seg(s,self.segdata) for s in self.parent.selectl.selectseg]
            self.parent.selectl.modeSMS()
            self.parent.selectl.multsel()
        self.close()

    def cancel(self):
        """ cancel click
        """
        self.close()

class AppForm(QMainWindow):
    def __init__(self,L, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Pylayers : Stand Alone Editor (Beta)')
        self.L=L

        self.create_main_frame()
        self.create_menu()

        # self.create_status_bar()
        # self.textbox.setText('1 2 3 4')
        self.on_draw()
        # self.shortcuts()

    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)


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


    def on_about(self):
        msg = """ A demo of using PyQt with matplotlib:
        
         * Use the matplotlib navigation bar
         * Add values to the text box and press Enter (or click "Draw")
         * Show or hide the grid
         * Drag the slider to modify the width of the bars
         * Save the plot to a file using the File menu
         * Click on a bar to receive an informative message
        """
        QMessageBox.about(self, "About the demo", msg.strip())
        # self.edit_properties()

    # def on_pick(self, event):
    #     # The event received here is of the type
    #     # matplotlib.backend_bases.PickEvent
    #     #
    #     # It carries lots of information, of which we're using
    #     # only a small amount here.
    #     # 
    #     box_points = event.artist.get_bbox().get_points()
    #     msg = "You've clicked on a bar with coords:\n %s" % box_points
        
    #     QMessageBox.information(self, "Click!", msg)
    



    def on_draw(self):
        """ Redraws the figure
        """
        # str = unicode(self.textbox.text())
        # self.data = map(int, str.split())
        
        # x = range(len(self.data))

        # clear the axes and redraw the plot anew
        #
        self.axes.clear()        
        # self.axes.grid(self.grid_cb.isChecked())
        self.L.display['nodes']=True
        self.L.display['ednodes']=True
        self.L.display['subsegnb']=True

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
        
        self.canvas.draw()

    def create_main_frame(self):
        self.main_frame = QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #

        self.dpi = 100
        self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)

        
        # Bind the 'pick' event for clicking on one of the bars
        #
        # self.canvas.mpl_connect('pick_event', self.on_pick)
        self.selectl = SelectL2(self.L,fig=self.fig,ax=self.axes)

        # self.cid1 = self.canvas.mpl_connect('button_press_event',
        #                                    self.selectl.OnClick)
        # self.cid2 = self.canvas.mpl_connect('key_press_event',
        #                                    self.selectl.OnPress)
        # self.cid3 = self.canvas.mpl_connect('key_release_event',
        #                                    self.selectl.OnRelease)
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
        self.canvas.setFocusPolicy( Qt.ClickFocus )
        self.canvas.setFocus()

        
        #Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        # # Other GUI controls



        # # 
        # self.properties = QComboBox()
        # self.properties.addItem("Ubuntu")
        # self.properties.addItem("Mandriva")
        # self.properties.addItem("Fedora")
        # self.properties.addItem("Red Hat")
        # self.properties.addItem("Gentoo")

        # self.properties.setMinimumWidth(200)
        # self.connect(self.properties, SIGNAL('editingFinished ()'), self.on_draw)




        # self.textbox = QLineEdit()
        # self.textbox.setMinimumWidth(200)
        # self.connect(self.textbox, SIGNAL('editingFinished ()'), self.on_draw)
        
        # self.draw_button = QPushButton("&Draw")
        # self.connect(self.draw_button, SIGNAL('clicked()'), self.on_draw)
        
        # self.grid_cb = QCheckBox("Show &Grid")
        # self.grid_cb.setChecked(False)
        # self.connect(self.grid_cb, SIGNAL('stateChanged(int)'), self.on_draw)
        
        # slider_label = QLabel('Bar width (%):')
        # self.slider = QSlider(Qt.Horizontal)
        # self.slider.setRange(1, 100)
        # self.slider.setValue(20)
        # self.slider.setTracking(True)
        # self.slider.setTickPosition(QSlider.TicksBothSides)
        # self.connect(self.slider, SIGNAL('valueChanged(int)'), self.on_draw)
        
        # #
        # # Layout with box sizers
        # # 
        # hbox = QHBoxLayout()
        # hbox.addWidget(self.properties)
        # hbox.setAlignment(self.properties, Qt.AlignVCenter)
        # for w in [  self.textbox, self.draw_button, self.grid_cb,
        #             slider_label, self.slider]:
        #     hbox.addWidget(w)
        #     hbox.setAlignment(w, Qt.AlignVCenter)
        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        # vbox.addLayout(hbox)

        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("This is a demo")
        self.statusBar().addWidget(self.status_text, 1)

        
    # def shortcuts(self):
    #     shortcut = QShortcut(self)
    #     shortcut.setKey("Ctrl+D")
    #     self.connect(shortcut, SIGNAL("activated()"), self.edit_properties )    # QtGui.QShortcut(QtGui.QKeySequence("Ctrl+Return"), self.myWidget, self.doSomething)

    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        self.edit_menu = self.menuBar().addMenu("&Edit")
        self.help_menu = self.menuBar().addMenu("&Help")
        # load_file_action = self.create_action("&Save plot",
        #     shortcut="Ctrl+S", slot=self.save_plot, 
        #     tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")

        refresh = self.create_action("&Refresh", slot=self.on_draw, 
            shortcut="F10", tip="Refresh the application")
        properties= self.create_action("&Properties", slot=self.edit_properties, 
            shortcut="F9", tip="Edit Wall properties")
        # show3= self.create_action("&Properties", slot=self.edit_properties, 
        #     shortcut="F9", tip="3D show")

        about_action = self.create_action("&About", 
            shortcut='F12', slot=self.on_about, 
            tip='about')
        

        self.add_actions(self.file_menu, 
            ( quit_action,))
        
        self.add_actions(self.edit_menu, 
            ( refresh,properties))

        self.add_actions(self.help_menu, (about_action,))

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


def main(L):
    app = QApplication(sys.argv)
    form = AppForm(L=L)
    form.show()
    app.exec_()


if __name__ == "__main__":
    L=Layout('TA-Office.ini')
    main(L)
