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


class PropertiesWin(QDialog):    # any super class is okay
    def __init__(self, slabDB,parent=None):
        super(PropertiesWin, self).__init__(parent)
        self.slabDB=slabDB
        # combo box 
        self._init_combo()
        self._init_slab_prop()
        

        self._init_layout()

        # self.button.clicked.connect(self.create_child)
    
    def _init_combo(self):
        self.combo = QComboBox()
        for s in self.slabDB.keys():
            self.combo.addItem(s)

    def _init_slab_prop(self):
        self.heightmin = QLineEdit()       
        self.heightmin.setObjectName("zmin")
        self.heightmin.setText("0.0")

        self.heightmax = QLineEdit()
        self.heightmax.setObjectName("zmax")
        self.heightmax.setText("3.0")

        self.transition = QPushButton("Transition")
        self.transition.setCheckable(True)
        # self.transition.setText("0")

        self.heightmin.setMinimumWidth(10)
        self.heightmax.setMinimumWidth(10)
        self.transition.setMinimumWidth(10)

    def _init_layout(self):

        # slab
        hbox1 = QHBoxLayout() 
        hbox1.addWidget(self.combo)

        # slab prop
        hbox2 = QHBoxLayout() 
        for w in [ self.heightmin,self.heightmax,self.transition]:
            hbox2.addWidget(w)
            hbox2.setAlignment(w, Qt.AlignVCenter)

        # validation

        buttono=QPushButton("OK")
        buttonc=QPushButton("Cancel")
        buttono.clicked.connect(self.valide)
        buttonc.clicked.connect(self.cancel)

        hbox3 = QHBoxLayout() 
        hbox3.addWidget(buttonc)
        hbox3.addWidget(buttono)

        vbox = QVBoxLayout()
        vbox.addLayout(hbox1)
        vbox.addLayout(hbox2)
        vbox.addLayout(hbox3)

        self.setLayout(vbox)


    def valide(self):
        """ ok click
        """
        zmin = self.heightmin.text()
        zmax = self.heightmax.text()
        if self.transition.isChecked():
            trans = True
        else :
            trans = False
        print zmin,zmax,trans

    def cancel(self):
        """ cancel click
        """
        self.close()


    def spinbox(self):
        """ plus minus button to add subseg
        """
        pass


    def create_child(self):
        # here put the code that creates the new window and shows it.
        child = PropertiesWin(self)
        child.show()

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
        self.prop = PropertiesWin(slabDB=self.L.sl)
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
        print self.selectl.nsel,self.selectl.state 
        if (self.selectl.state == 'SS') and (self.selectl.nsel > 0):
            self.prop.show()
        elif (self.selectl.state == 'SMS') and (self.selectl.selectseg!=[]):
            self.prop.show()
    def on_about(self):
        msg = """ A demo of using PyQt with matplotlib:
        
         * Use the matplotlib navigation bar
         * Add values to the text box and press Enter (or click "Draw")
         * Show or hide the grid
         * Drag the slider to modify the width of the bars
         * Save the plot to a file using the File menu
         * Click on a bar to receive an informative message
        """
        # QMessageBox.about(self, "About the demo", msg.strip())
        self.edit_properties()

    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        # 
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        
        QMessageBox.information(self, "Click!", msg)
    



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

        self.cid1 = self.canvas.mpl_connect('button_press_event',
                                           self.selectl.OnClick)
        self.cid2 = self.canvas.mpl_connect('key_press_event',
                                           self.selectl.OnPress)
        self.cid3 = self.canvas.mpl_connect('key_release_event',
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
        
        # load_file_action = self.create_action("&Save plot",
        #     shortcut="Ctrl+S", slot=self.save_plot, 
        #     tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        refresh = self.create_action("&Refresh", slot=self.on_draw, 
            shortcut="F10", tip="Refresh the application")
        # save_layout = self.create_action("&MS", slot=self.selectl, 
        #      shortcut="F1", tip="Multiple Selection")
        properties= self.create_action("&Properties", slot=self.edit_properties, 
            shortcut="F9", tip="Edit Wall properties")
        self.add_actions(self.file_menu, 
            ( quit_action,refresh,properties))
        
        self.help_menu = self.menuBar().addMenu("&Help")

        about_action = self.create_action("&About", 
            shortcut='F12', slot=self.on_about, 
            tip='About the demo')
        
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
