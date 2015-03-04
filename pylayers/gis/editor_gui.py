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
from pylayers.util.project import *
import pylayers.util.pyutil as pyu


class SubSegWin(QDialog):    # any super class is okay
    def __init__(self,Nss=1,zmin=0.,zmax=3.0,subsegdata={},parent=None):
        super(SubSegWin, self).__init__(parent)
        #
        self.gparent=parent.parent
        self.parent=parent
        # mulsti segment selection indicator
        self.mulseg=parent.mulseg
        # dictionnary to pass subseg data
        self.subsegdata=parent.subsegdata
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
            when new subseg
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
        self.parent.subsegdata={}
        self.parent.subsegdata['ss_name']=[]
        self.parent.subsegdata['ss_z']=[]
        self.parent.subsegdata['ss_offset']=[]
        # # check z
        # zz=[]
        # for ss in range(self.Nss):
        #     zz.append((self.lheightmin[ss].value(),self.lheightmax[ss].value()))
        # zz=np.array(zz)
        # szz = np.sort(zz,axis=0)
        # # position overla
        # uo=np.where(szz[:-1,1]>szz[1:,0])
        
        if len(uo) == 0: 
            for ss in range(self.Nss):
                z = (self.lheightmin[ss].value(),self.lheightmax[ss].value())
                self.parent.subsegdata['ss_name'].append(str(self.lcomboslab[ss].currentText()))
                self.parent.subsegdata['ss_z'].append(z)
                self.parent.subsegdata['ss_offset'].append(self.loffset[ss].value())

            # if not self.mulseg:
            #     self.gparent.L.edit_seg(self.gparent.selectl.nsel,self.subsegdata)
            # else:
            #     [self.gparent.L.edit_seg(s,self.subsegdata) for s in self.gparent.selectl.selectseg]
            self.close()
        else :



    def cancel(self):
        self.close()


class PropertiesWin(QDialog):    # any super class is okay
    def __init__(self,mulseg=False,parent=None):
        super(PropertiesWin, self).__init__(parent)
        # to imporve here. Probably something to inherit from parent App
        self.parent=parent
        # determine if multiple se gments are selected
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
        else : 
            Nss=0

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
        print 'SSDATA :',self.subsegdata
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

        self.subsegdata.update({'ss_name':self.subsegdata['ss_name'][:Nss],
                              'ss_offset':self.subsegdata['ss_offset'][:Nss],
                              'ss_z':self.subsegdata['ss_z'][:Nss]})

        if not self.mulseg:
            self.parent.L.edit_seg(self.parent.selectl.nsel,self.segdata)
            self.parent.L.update_sseg(self.parent.selectl.nsel,self.subsegdata)
            self.parent.selectl.modeIni()
        else:
            [self.parent.L.edit_seg(s,self.segdata) for s in self.parent.selectl.selectseg]
            [self.parent.L.update_sseg(s,self.subsegdata) for s in self.parent.selectl.selectseg]
            self.parent.selectl.modeSMS()
            self.parent.selectl.multsel()
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
        self.parent.fig.clear()
        self.parent.fig.canvas.draw()
        del self.parent.L
        del self.parent.main_frame
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


class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Pylayers : Stand Alone Editor (Beta)')
        self.filename=''
        self.create_menu()
        self.create_status_bar()
        # self.shortcuts()


    def new(self):
        self.closel()
        self.L=Layout('void.ini')
        self.filename=''
        self.create_main_frame()
        self.on_draw()

    def open(self):
        filename = QFileDialog.getOpenFileName(self,'Open Pylayers Layout File',pyu.getlong('',pstruc['DIRINI']),'(*.ini);;(*.osm)')
        try:
            _filename= pyu.getshort(str(filename))
            self.L=Layout(_filename)
            self.filename=self.L.filename
            self.create_main_frame()
            self.on_draw()
            self.setWindowTitle(self.L.filename + '- Pylayers : Stand Alone Editor (Beta)')
            self.resize(self.fig.canvas.width(),self.fig.canvas.height())
            print 'loaded'
        except:
            pass

    def save(self,force=False):

        if self.filename == '' or force:
            filename = QFileDialog.getSaveFileName(self, 'Save Layout', pyu.getlong('',pstruc['DIRINI']),'*.ini')
            try:
                _filename= pyu.getshort(str(filename))
            except:
                pass
        else :
            _filename=self.L.filename

        try:
            self.L.saveini(_filename)
            self.L.saveosm(_filename.split('.')[0] + '.osm')
            self.L=Layout(_filename)
            self.filename=self.L.filename
            self.setWindowTitle(self.L.filename + '- Pylayers : Stand Alone Editor (Beta)')
            print 'saved'
        except:
            pass

    def closel(self,exit=False):
        dial_res=''
        self.sq = SaveQuitWin(parent=self,exit=exit)
        self.sq.show()

    def exitl(self):
        plt.rcParams.update(self.selectl.rcconf)
        self.selectl.fig.canvas.mpl_disconnect(self.cid1)
        self.selectl.fig.canvas.mpl_disconnect(self.cid2)
        self.selectl.fig.canvas.mpl_disconnect(self.cid3)
        self.selectl.fig.canvas.mpl_disconnect(self.cid4)
        self.selectl.fig.canvas.mpl_disconnect(self.cid5)
        self.selectl.fig.canvas.mpl_disconnect(self.cid6)
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


    def on_about(self):
        msg = """ This is the PyLayers' Stand-Alone Layout Editor (BETA)
        
         This tool allow to personalize  your own building layout in terms of
         floor plan and constitutive materials.
         Once edited and saved, the layout can be used into the PyLayers Ray tracing tool.



         Shortcuts:
         ----------
         F1 : Select mode
         F2 : New point mode
         F3 : Edit segment properties

         CTRL + o : Open Layout
         CTRL + s : Save Layout
         CTRL + q : Quit Editor
         escape : back to a stable state
 
         More hints about editing can be found in the status bar.



         Thank you for using Pylayers and this tool

         The Pylayers' Dev Team
         www.pylayers.org

        """
        QMessageBox.about(self, "Pylayers' Stand-Alone Layout Editor (BETA)", msg.strip())

    



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

    def on_release(self,event):
        string=''
        try:
            string = string + ' ' + self.L.Gs.node[self.selectl.nsel]['name']
        except:
            pass
        try:
            string = string + ' with ' +str(len(self.L.Gs.node[self.nsel]['ss_name'])) + 'subseg(s)'
        except:
            pass
        string = string +'\t'+self.selectl.help[self.selectl.state]
        self.statusBar().showMessage(string)
    def create_main_frame(self):
        self.main_frame = QWidget()
        
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
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        # vbox.addLayout(hbox)

        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("Open a Layout")
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
        new_action = self.create_action("&New Layout", slot=self.new, 
        shortcut="Ctrl+n", tip="new layout")
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

        refresh = self.create_action("&Refresh", slot=self.on_draw, 
            shortcut="F10", tip="Refresh the application")
        properties= self.create_action("&Properties", slot=self.edit_properties, 
            shortcut="F3", tip="Edit Wall properties")
        # show3= self.create_action("&Properties", slot=self.edit_properties, 
        #     shortcut="F9", tip="3D show")

        about_action = self.create_action("&About", 
            shortcut='F12', slot=self.on_about, 
            tip='about')
        

        self.add_actions(self.file_menu, 
            ( new_action,open_action,None,save_action,saveas_action,None,close_action,quit_action,))
        
        self.add_actions(self.edit_menu, 
            ( properties,None,refresh))

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


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    # form.setGeometry(100,100,300,300)
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()