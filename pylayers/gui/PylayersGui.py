# -*- coding: utf-8 -*-

from pylayers.simul.link import *
import pylayers.util.pyutil as pyu
import pylayers.signal.standard as std
from pylayers.util.project import *
import json


# TEST
import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT
from matplotlib.figure import Figure
from pyface.qt import QtGui,QtCore
from traitsui.qt4.editor import Editor
from traitsui.qt4.basic_editor_factory import BasicEditorFactory



# console ipython 
from IPython import embed_kernel





from traits.api import HasTraits, Button,Range,Enum, Instance, \
        on_trait_change,property_depends_on,Float,Str,Int,Bool,List
from traitsui.api import View, Item,HSplit,VSplit, RangeEditor, \
                        EnumEditor,Group,spring,HGroup,VGroup,Handler, \
                        InstanceEditor


from traitsui.menu import Action, ActionGroup, Menu, MenuBar, ToolBar


from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor, \
                MlabSceneModel
from tvtk.pyface.api import Scene




try:
    get_ipython
except NameError:
    banner=exit_msg=''
else:
    banner = '*** Nested interpreter ***'
    exit_msg = '*** Back in main IPython ***'

# First import the embed function
from IPython.frontend.terminal.embed import InteractiveShellEmbed

## INIT DLink object 

DL=DLink()

filename=pyu.getlong('wstd.json',pstruc['DIRSIMUL'])
fp = open(filename)
stds = json.load(fp)
av_wstds = ['None']+ stds.keys()



dchann = {w:[str(i) for i in std.Wstandard(w).chan.keys()] for w in av_wstds if w !='None'}
dchann.update({'None':['None']})




from qtconsole.rich_ipython_widget import RichJupyterWidget
from qtconsole.inprocess import QtInProcessKernelManager
from IPython.lib import guisupport

class QIPythonWidget(RichJupyterWidget):
    """ Convenience class for a live IPython console widget. We can replace the standard banner using the customBanner argument"""
    def __init__(self,customBanner=None,*args,**kwargs):
        if not customBanner is None: self.banner=customBanner
        super(QIPythonWidget, self).__init__(*args,**kwargs)
        self.kernel_manager = kernel_manager = QtInProcessKernelManager()

        kernel_manager.start_kernel()

        kernel_manager.kernel.gui = 'qt4'
        self.kernel_client = kernel_client = self._kernel_manager.client()
        kernel_client.start_channels()

        def stop():
            kernel_client.stop_channels()
            kernel_manager.shutdown_kernel()
            guisupport.get_app_qt4().exit()            
        self.exit_requested.connect(stop)

    def pushVariables(self,variableDict):
        """ Given a dictionary containing name / value pairs, push those variables to the IPython console widget """
        self.kernel_manager.kernel.shell.push(variableDict)
    def clearTerminal(self):
        """ Clears the terminal """
        self._control.clear()    
    def printText(self,text):
        """ Prints some plain text to the console """
        self._append_plain_text(text)        
    def executeCommand(self,command):
        """ Execute a command in the frame of the console widget """
        self._execute(command,False)


class JupyterWidget(QtGui.QWidget):
    """ Main GUI Widget including a button and IPython Console widget 
    inside vertical layout 
    """
    def __init__(self, parent=None):
        super(JupyterWidget, self).__init__(parent)
        layout = QtGui.QVBoxLayout(self)
        ipyConsole = QIPythonWidget()
        layout.addWidget(ipyConsole)
        # ipyConsole.pushVariables({'DL':DL})
        allvar = globals()
        allvar.update(locals())
        ipyConsole.pushVariables(allvar)

class _MPLFigureEditor(Editor):

    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # matplotlib commands to create a canvas 
        frame = QtGui.QWidget()
        mpl_canvas = FigureCanvas(self.value)
        mpl_toolbar = NavigationToolbar2QT(parent=frame,canvas = mpl_canvas)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(mpl_canvas)
        vbox.addWidget(mpl_toolbar)
        frame.setLayout(vbox)
        mpl_canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        mpl_canvas.setFocus()

        return frame#mpl_canvas

class MPLFigureEditor(BasicEditorFactory):

    klass = _MPLFigureEditor

class WstdHandler(Handler):

    channels = List(Str)

    def object_Wstd_Enum_changed(self, info):
        """
        This method listens for a change in the *state* attribute of the
        object (Address) being viewed.

        When this listener method is called, *info.object* is a reference to
        the viewed object (Address).

        """
        # Change the list of available cities

        self.channels = dchann[info.object.Wstd_Enum]

        # As default value, use the first city in the list:
        info.object.chann = self.channels[0]
        # info.object.DL.fGHz = 



class PylayersGUI(HasTraits):

    # slider/dropdown widgets etc

    # Layout
    laynames = [''] + np.sort(os.listdir(basename +'/struc/lay/')).tolist()#['','DLR.lay','defstr.lay','TC2_METIS.lay']#,
    Lay_Enum = Enum(laynames)


    ## Antenna file :
    av_ant = ['Omni','Gauss','aperture']
    antext= ['vsh3','sh3']
    for fname in os.listdir(basename +'/ant'):
        if fname.split('.')[-1] in antext:
            av_ant.append(fname)

    # Init Positions
    xmin = DL.L.ax[0]
    xmax = DL.L.ax[1]
    ymin = DL.L.ax[2]
    ymax = DL.L.ax[3]
    zmin = 0.
    zmax = DL.L.maxheight-0.1
    # Antenna 

    ## position a
    aX    = Range(low='xmin',high='xmax',value= float(xmin+xmax/2.))
    aY    = Range(low='ymin',high='ymax',value= float(ymin+ymax/2.))
    aZ    = Range(low='zmin',high='zmax',value= float(zmin+zmax/2.))

    ## rotation a
    agamma    = Range(float(-3.14), float(3.14), 0., )#mode='spinner')
    abeta    = Range(float(-3.14), float(3.14), 0., )#mode='spinner')
    aalpha    = Range(float(-3.14), float(3.14), 0., )#mode='spinner')

    ## file a:

    a_ant = Enum(av_ant)


    # Antenna B
    ## position b
    bX    = Range(low='xmin',high='xmax',value= float(xmin+xmax/2.))
    bY    = Range(low='ymin',high='ymax',value= float(ymin+ymax/2.))
    bZ    = Range(low='zmin',high='zmax',value= float(zmin+zmax/2.))

    ## rotation b

    bgamma    = Range(float(-3.14), float(3.14), 0., )#mode='spinner')
    bbeta    = Range(float(-3.14), float(3.14), 0., )#mode='spinner')
    balpha    = Range(float(-3.14), float(3.14), 0., )#mode='spinner')

    ## file b:
    b_ant = Enum(av_ant)


    # frequency

    fmmin = 0.
    fmmax = 80.

    

    fmin=Range(low = 'fmmin', high = 'fmmax',value = float(DL.Aa.fGHz[0]) )
    fmax=Range(low = 'fmmin', high = 'fmmax',value = float(DL.Aa.fGHz[-1]) )
    fstep=Range(low = 0.001,high = 100, value = 0.5)

    # advanced




    # init interface
    scene = Instance(MlabSceneModel, ())

    plot = Instance(PipelineBase)

    # @on_trait_change('scene.activated')
    # def init_plot(self):
    #     DL._show3()
    # When the scene is activated, or when the parameters are changed, we
    # update the plot.


    # def _open_changed(self):
    #     """ Handles the user clicking the 'Open...' button.
    #     """
    #     path = pyu.getlong('',pstruc['DIRSTR'])
    #     file_name = open_file(file_name= path ,extensions = FileInfo())
    #     if file_name != '':
    #         self.file_name = file_name


    @on_trait_change('Lay_Enum')
    def update_L(self):
        if self.Lay_Enum != ' ':
            mlab.clf()
            DL.L=Layout(self.Lay_Enum,bgraphs=True)

            self.xmin=DL.L.ax[0]
            self.xmax=DL.L.ax[1]
            self.ymin=DL.L.ax[2]
            self.ymax=DL.L.ax[3]
            self.zmin=0.
            self.zmax=DL.L.maxheight-0.1

            self.aX,self.aY,self.aZ=DL.a
            self.bX,self.bY,self.bZ=DL.b

            DL.a= np.array([self.aX,self.aY,self.aZ])
            DL.b= np.array([self.bX,self.bY,self.bZ])
            self.cutoff = DL.cutoff
            if not hasattr(DL,'_maya_fig'):
                DL._show3()

    @on_trait_change('aX,aY,aZ')
    def update_a(self):
        """ update position ant a
        """
        self.clear_fig()
        DL.a= np.array([self.aX,self.aY,self.aZ])
        self.cutoff = DL.cutoff


    @on_trait_change('bX,bY,bZ')
    def update_b(self):
        """ update position ant b
        """
        self.clear_fig()
        DL.b= np.array([self.bX,self.bY,self.bZ])
        self.cutoff = DL.cutoff

    @on_trait_change('aalpha,abeta,agamma')
    def update_Ta(self):
        """ update rot ant a
        """
        T = geu.MEulerAngle(self.aalpha,beta=self.abeta,gamma=self.agamma)
        DL.Ta=T
        self.clear_fig()
        # if DL.dexist['Ct']['exist']:
        #     DL.C.locbas(Tt=DL.Ta, Tr=DL.Tb)
        #     #T channel
        #     DL.H = DL.C.prop2tran(a=DL.Aa,b=DL.Ab,Friis=True)
        #     self.plt_all()

    @on_trait_change('balpha,bbeta,bgamma')
    def update_Tb(self):
        """ update rot ant b
        """
        T = geu.MEulerAngle(self.balpha,beta=self.bbeta,gamma=self.bgamma)
        DL.Tb=T
        self.clear_fig()


    @on_trait_change('a_ant,fmin,fmax,fstep')
    def update_Aa(self):
        DL.Aa=Antenna(self.a_ant)
        self.clear_fig()
        # if DL.Aa.fromfile:
        #     self.fmin=DL.Aa.fGHz[0]
        #     self.fmax=DL.Aa.fGHz[-1]
        #     self.fstep=min(1,DL.Aa.fGHz[1]-DL.Aa.fGHz[0])

    @on_trait_change('b_ant,fmin,fmax,fstep')
    def update_Ab(self):
        DL.Ab=Antenna(self.b_ant)
        self.clear_fig()

        # if DL.Ab.fromfile:
        #     self.fmin=DL.Ab.fGHz[0]
        #     self.fmax=DL.Ab.fGHz[-1]
        #     self.fstep=min(1,DL.Ab.fGHz[1]-DL.Ab.fGHz[0])

    @on_trait_change('fmin,fmax,fstep,chann')
    def update_fGHz(self):
        if self.Wstd_Enum != 'None':
            W=std.Wstandard(self.Wstd_Enum)
            # DL.fGHz = W.chan[eval(self.chann)].fghz
            Wchan = W.chan[eval(self.chann)]
            fcGHz = Wchan['fcGHz']
            BWGHz = Wchan['BMHz']
            GMHz = Wchan['GMHz']
            fGHz = Wchan.fghz

            DL.fGHz = np.array([fcGHz])

            self.BWGHz = BWGHz
            self.fmin = float(fGHz[0])
            self.fmax = float(fGHz[-1])
            self.fstep = float(fGHz[1]-fGHz[0])

        else:
            if self.fmin < self.fmax:
                DL.fGHz = np.arange(self.fmin,
                                    self.fmax,
                                    self.fstep
                                    )
            elif self.fmin == self.fmax:
                DL.fGHz=np.array([self.fmin])
            self.BWGHz = 5


    @on_trait_change('Beval')
    def DLeval(self):
    #     clear_output()
        DL.eval(verbose=False,
                force=self.force,
                cutoff=self.cutoff,
                threshold=self.threshold,
                diffraction = self.diffraction,
                nD=self.nD,
                nT=self.nT,
                nR=self.nR,
                applywav = self.applywav)
        DL._update_show3(delrays=True)
        ER = np.squeeze(DL.H.energy())
        DL.R._show3(ER=ER)
        self.plt_all()



    def plt_all(self):
        self.plt_cir()
        self.plt_doa()
        self.plt_dod()
        self.plt_dspread()
        self.plt_aspread()

    def plt_cir(self):
        self.figcir.clf()
        ax = self.figcir.add_subplot(111)
        DL.plt_cir(fig=self.figcir,ax=ax,BWGHz=self.BWGHz,Nf = 5000 )
        # ir = DL.H.getcir(BWGHz=5,Nf=1000)
        # ir.plot(fig=self.figcir,ax=ax)
        # ax.plot(DL.H.taud,20*np.log10(DL.H.y[:,0,0,0]),'or')
        self.figcir.canvas.draw()
        # DL.plt_doadod(d='doa')
        # DL.H.plot(fig=self.figcir,ax=ax)
        # self.figcir.canvas.draw()

    def plt_doa(self):
        self.figdoa.clf()
        ax = self.figdoa.add_subplot(111,polar=True)
        # DL.L.showG('s',ax=ax,fig=self.figure)
        # DL.H.plotd(d='doa',polar=True,fig=self.figdoa,ax=ax)
        DL.plt_doa(polar=True,fig=self.figdoa,ax=ax)
        self.figdoa.canvas.draw()

    def plt_dod(self):
        self.figdod.clf()
        ax = self.figdod.add_subplot(111,polar=True)
        DL.plt_dod(polar=True,fig=self.figdod,ax=ax)
        # DL.L.showG('s',ax=ax,fig=self.figure)
        # DL.H.plotd(d='dod',polar=True,fig=self.figdod,ax=ax)
        self.figdod.canvas.draw()

    def plt_dspread(self):
        self.figds.clf()
        ax = self.figds.add_subplot(111)
        DL.plt_dspread(fig=self.figds,ax=ax)
        self.figds.canvas.draw()

    def plt_aspread(self):
        self.figas.clf()
        ax = self.figas.add_subplot(111)
        DL.plt_aspread(fig=self.figas,ax=ax)
        self.figas.canvas.draw()


    def clear_fig(self,lf=['cir','doa','dod','as','ds']):
        for f in lf:
            eval('self.fig'+f+'.clf()')
            eval('self.fig'+f+'.canvas.draw()')


    #####
    ##### RENDERING 3D MAYAVI
    #####

    render3d = Item('scene', editor=SceneEditor(scene_class=Scene),
                     height=500, width=1500, show_label=False)



    # ###
    # ### Matplotlib figure
    # ###

    # figure = Instance(Figure(figsize=(8,20)), ())


    #####
    ##### Layout SELECTION
    #####


    # Layout
    GLay =  Group(Item('Lay_Enum',
                            style='simple',
                            label='file'),
                 show_labels=False,
                 label='Layout')



    #####
    ##### WIRELESS STANDARD
    #####


    # wireless standard
    Wstd_Enum = Enum('None', av_wstds)

    chann = Str

    # chann = Enum(av_chann)

    GWstd_None = Group(Item('fmin',label='fGHz min',style='text'),
                       Item('fmax',label='fGHz max',style='text'),
                       Item('fstep',label='fGHz step',style='text'),
                       label = 'Frequency',
                       show_border= True,
                       enabled_when = 'Wstd_Enum == \'None\''
                       )



    GWstd_std = Group(Item(name ='chann',editor=EnumEditor(name='handler.channels')
                           )  ,
                       label = 'channel',
                       show_border= True,
                       enabled_when = 'Wstd_Enum != \'None\''
                       )



    GWstd = Group(
                  Group(Item (name = 'Wstd_Enum',
                        label = 'Wireless Standard')),
                  GWstd_None,
                  GWstd_std,
                  label='Wireless Standard',
                  show_labels=True,
                  show_border=False)


    #####
    ##### ANTENNA
    ##### 

    xmin=Float
    xmax = Float
    ymin=Float
    ymax = Float
    zmin=Float
    zmax = Float


    # Ant A file

    Iax = Item('aX',
                 editor=RangeEditor(low_name='xmin',
                                    high_name='xmax',
                                    format='%.1f',
                                    label_width=28,
                                    mode='auto'),
                 label='x'
                 )
    Iay = Item('aY',
                 editor=RangeEditor(low_name='ymin',
                                    high_name='ymax',
                                    format='%.1f',
                                    label_width=28,
                                    mode='auto'),
                 label='y'
                 )
    Iaz = Item('aZ',
                 editor=RangeEditor(low_name='zmin',
                                    high_name='zmax',
                                    format='%.1f',
                                    label_width=28,
                                    mode='auto'),
                 label='z'
                 )
    GPos_a = VGroup(
               Iax,
               Iay,
               Iaz,
               id = 'a',
               label = 'Position',
               show_border=True,
               show_labels=True,
               layout='split'
               )


    Ifile_a = Item('a_ant',label='file')


    GRot_a = VGroup(
               Item('agamma',label='x-roll'),
               Item('abeta',label='y-roll'),
               Item('aalpha',label='z-roll'),
               id = 'Ta',
               label = 'Rotation',
               show_border=True,
               layout='split'
               )


    G_a = Group(Ifile_a,
                       GPos_a,
                       GRot_a,
                label='Antenna a',
                show_border=False
                )


    #### ANtenna B

    # Ant B positions

    Ibx = Item('bX',
                 editor=RangeEditor(low_name='xmin',
                                    high_name='xmax',
                                    format='%.1f',
                                    label_width=28,
                                    mode='auto'),
                 label='x'
                 )
    Iby = Item('bY',
                 editor=RangeEditor(low_name='ymin',
                                    high_name='ymax',
                                    format='%.1f',
                                    label_width=28,
                                    mode='auto'),
                 label='y'
                 )
    Ibz = Item('bZ',
                 editor=RangeEditor(low_name='zmin',
                                    high_name='zmax',
                                    format='%.1f',
                                    label_width=28,
                                    mode='auto'),
                 label='z'
                 )

    GPos_b = Group(
               Ibx,
               Iby,
               Ibz,
               id = 'b',
               label = 'Position',
               show_border=True,
               layout='split'
               )

    # Ant B file

    Ifile_b = Item('b_ant',label='file')


    GRot_b = Group(
               Item('bgamma',label='x-roll'),
               Item('bbeta',label='y-roll'),
               Item('balpha',label='z-roll'),
               id = 'Tb',
               label = 'Rotation',
               show_border=True,
               layout='split'
               )

    G_b = Group(Ifile_b,
                       GPos_b,
                       GRot_b,
                label='Antenna b',
                show_border=False,
                )




    #### 
    #### advanced CONFIRGURATION
    #### 


    force =Bool
    diffraction = Bool
    applywav = Bool
    applywav = Bool
    low_cutoff = 1
    high_cutoff = 30
    cutoff = Range(low='low_cutoff',high='high_cutoff',value=DL.cutoff)
    threshold = Range(0,1.,0.8)
    nD=2
    nR=10
    nT=10

    G_advanced = Group(VGroup(
                      Item('force',
                            label='force',
                            resizable=False,
                            style='simple'),
                       Item('cutoff',
                            label='cutoff',
                            editor=RangeEditor(low_name='low_cutoff',
                                    high_name='high_cutoff',
                                    label_width=28,
                                    mode='auto'),
                            width=0.2,
                            style='simple'),
                       Item('threshold',
                            label='threshold',
                            editor = RangeEditor(format='%.1f'),
                            width=0.2,
                            style='simple'),
                       Item('diffraction',
                            label='diffractions',
                            style='simple'),
                       Item('nD',
                            label='max nb Diffractions',
                            enabled_when='diffraction' ,
                            style='simple'),
                       Item('nR',
                            label='max nb Reflections',
                            style='simple'),
                       Item('nT',
                            label='max nb Transmissions',
                            style='simple'),
                       Item('applywav',
                            label='applywav',
                            style='simple'),
                       label='Ray Tracing Configuration',
                       show_labels=True,
                       show_border=False))


    ####
    ### MANAGING GROUPS
    ###



    # LEFT GROUP WINDOW
    Beval = Button('Launch Ray-Tracing')

    GLeft = Group(
                  GLay,
                  GWstd, 
                  G_advanced
                  )

    # Antenna GRoup 
    GAnt_ab = HGroup(spring,G_a,spring,G_b,spring)


    GAnt_Eval = Group(GAnt_ab,
                      HGroup(spring,
                             Item('Beval',
                                  enabled_when='Lay_Enum != \'\''
                                  ),
                            show_labels=False)
                      )


    #### TOP GROUP

    GR_0= HSplit(GLeft, 
                 render3d,
                 layout='split')



   # BOTTOM GROUP

    figcir= Instance(Figure(figsize=(8,20)), ())
    figdoa= Instance(Figure(figsize=(8,20)), ())
    figdod= Instance(Figure(figsize=(8,20)), ())
    figas= Instance(Figure(figsize=(8,20)), ())
    figds= Instance(Figure(figsize=(8,20)), ())


    GExploit = Group ( Group(Item('figcir',
                                  editor=MPLFigureEditor(),
                                 ),
                             label='CIR'),
                     Group(Item('figdoa',
                                  editor=MPLFigureEditor()
                                 ),
                             label='DOA'),
                     Group(Item('figdod',
                                  editor=MPLFigureEditor()
                                 ),
                             label='DOD'),
                     Group(Item('figas',
                                  editor=MPLFigureEditor()
                                 ),
                             label='Ang. Spread'),
                     Group(Item('figds',
                                  editor=MPLFigureEditor()
                                 ),
                             label='Delay Spread'),
                    layout='tabbed',
                    )

    GR_1 = HGroup(spring,GAnt_Eval,spring,GExploit)

    JWidget = JupyterWidget()
    JWidget.show()


    view = View(VGroup(GR_0,GR_1),
                # menubar=MenuBar(Menu_file),
                buttons=['Quit'],
                title="Pylayers GUI - beta",
                resizable=True,
                width=1., height=1.,
                handler=WstdHandler)


if __name__ == '__main__':
    gui = PylayersGUI()
    gui.configure_traits()


