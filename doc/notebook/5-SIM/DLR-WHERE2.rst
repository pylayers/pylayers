
.. code:: python

    from IPython.core.display import HTML
    
    def css_styling():
        styles = open("../styles/custom.css", "r").read()
        return HTML(styles)
    css_styling()



.. raw:: html

    <style>
        @font-face {
            font-family: "Computer Modern";
            src: url('http://mirrors.ctan.org/fonts/cm-unicode/fonts/otf/cmunss.otf');
        }
        div.cell{
            width:800px;
            margin-left:16% !important;
            margin-right:auto;
        }
        h1 {
            font-family: Helvetica, serif;
        }
        h4{
            margin-top:12px;
            margin-bottom: 3px;
           }
        div.text_cell_render{
            font-family: Computer Modern, "Helvetica Neue", Arial, Helvetica, Geneva, sans-serif;
            line-height: 145%;
            font-size: 130%;
            width:800px;
            margin-left:auto;
            margin-right:auto;
        }
        .CodeMirror{
                font-family: "Source Code Pro", source-code-pro,Consolas, monospace;
        }
        .prompt{
            display: None;
        }
        .text_cell_render h5 {
            font-weight: 300;
            font-size: 22pt;
            color: #4057A1;
            font-style: italic;
            margin-bottom: .5em;
            margin-top: 0.5em;
            display: block;
        }
        
        .warning{
            color: rgb( 240, 20, 20 )
            }  
    </style>
    <script>
        MathJax.Hub.Config({
                            TeX: {
                               extensions: ["AMSmath.js"]
                               },
                    tex2jax: {
                        inlineMath: [ ['$','$'], ["\\(","\\)"] ],
                        displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
                    },
                    displayAlign: 'center', // Change this to 'center' to center equations.
                    "HTML-CSS": {
                        styles: {'.MathJax_Display': {"margin": 4}}
                    }
            });
    </script>



Example of a trajectory synthesis in DLR WHERE2 environment
-----------------------------------------------------------

The document describing the measurement performed in the DLR environment
is `WHERE2 Deliverable D4.4 Test and Evaluation of the Integrated System
under Laboratory
Conditions <http://www.kn-s.dlr.de/where2/documents/Deliverables/Deliverable-D4.4.pdf>`__

.. code:: python

    from pylayers.simul.simulem import *
    from pylayers.antprop.rays import *
    from pylayers.antprop.channel import *
    from pylayers.antprop.signature import *
    import pylayers.util.pyutil as pyu
    from pylayers.gis.layout import *
    from pylayers.util.project import *
    import pylayers.signal.bsignal as bs
    from datetime import datetime
    import time
    import pdb
    import pickle
    import numpy as np
    import matplotlib.pyplot as plt 
    %matplotlib inline

.. parsed-literal::

    WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.


This function has to be moved in simulem module. It is a temporary
implementation. Signatures can be handled much more efficiently here. It
run a full simulation and returns a list of channel impulse response.

.. code:: python

    def evalcir(S,wav,cutoff=4):
        """
        Parameters
        ----------
    
        S 
        tx
        rx
        wav
        cutoff
    p.rays import *
    from pylayers.antpr
        """
    
        crxp =-1
        ctxp =-1
        tcir = {}
        tx = S.tx.position 
        Ntx = len(tx[0])
        rx = S.rx.position
        Nrx = len(rx[0])
    
        #for kt in range(1,Ntx-1):
        #print kt+1
        kt=0
        tcir[kt+1] = {}
        t = np.array([S.tx.position[0,kt+1],S.tx.position[1,kt+1],S.tx.position[2,kt+1]])
        for kr in range(Nrx-1):
            if (mod(kr,10)==0):
                print kr+1
            r = np.array([S.rx.position[0,kr+1],S.rx.position[1,kr+1],S.rx.position[2,kr+1]])
            ctx = S.L.pt2cy(t)
            crx = S.L.pt2cy(r)
            if (ctx<>ctxp)|(crx<>crxp):
                Si  = Signatures(S.L,ctx,crx)
                ctxp = ctx
                crxp = crx
                Si.run5(cutoff=cutoff,algo
    ='old')
            r2d = Si.rays(t,r)
            #r2d.show(S.L)
    
            r3d = r2d.to3D(S.L)
            r3d.locbas(S.L)
            r3d.fillinter(S.L)
            Ct  = r3d.eval(S.fGHz)
            sca = Ct.prop2tran(S.tx.A,S.rx.A)
            cir = sca.applywavB(wav.sfg)
            tcir[kt+1][kr+1]=cir
        return(tcir)
Loading the Layout
^^^^^^^^^^^^^^^^^^

.. code:: python

    S = Simul()
    filestr = 'DLR2'
    S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
    try:
        S.L.dumpr()
    except:
        S.L.build()
        S.L.dumpw()

::


    ---------------------------------------------------------------------------
    NoOptionError                             Traceback (most recent call last)

    <ipython-input-4-f57e44fbc7e2> in <module>()
    ----> 1 S = Simul()
          2 filestr = 'DLR2'
          3 S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
          4 try:
          5     S.L.dumpr()


    /home/uguen/Documents/rch/devel/pylayers/pylayers/simul/simulem.pyc in __init__(self, _filesimul)
        827         self.filemat = self.filematini.replace('.ini','.mat')
        828         self.fileslab = self.fileslabini.replace('.ini','.slab')
    --> 829         self.slab=SlabDB(self.filematini, self.fileslabini)
        830         self.filestr = 'defstr.str2'
        831         #


    /home/uguen/Documents/rch/devel/pylayers/pylayers/antprop/slab.pyc in __init__(self, filemat, fileslab)
       1865             self.mat.load(filemat)
       1866         if (fileslab != ''):
    -> 1867             self.load(fileslab)
       1868             self.dass()
       1869 


    /home/uguen/Documents/rch/devel/pylayers/pylayers/antprop/slab.pyc in load(self, _fileini)
       2114         for slabname in self.di.values():
       2115             S=Slab(name=slabname,mat=self.mat)
    -> 2116             S['lmatname']=eval(config.get(slabname,'lmatname'))
       2117             S['nbmat']=len(S['lmatname'])
       2118             S['color']=config.get(slabname,'color')


    /home/uguen/anaconda/lib/python2.7/ConfigParser.pyc in get(self, section, option, raw, vars)
        616             value = d[option]
        617         except KeyError:
    --> 618             raise NoOptionError(option, section)
        619 
        620         if raw or value is None:


    NoOptionError: No option 'lmatname' in section: 'ABSORBENT'


.. code:: python

    S.L.display['ednodes']=False
    S.L.display['nodes']=False
    S.L.display['title']='DLR WP4 WHERE2 measurement site'
    S.L.display['overlay']=False
    fig,ax = S.L.showGs()    


.. image:: DLR-WHERE2_files/DLR-WHERE2_8_0.png


.. code:: python

    S.show3()
We have a list of static Anchor Nodes. Those values correspond to the
actual anchor nodes coordinates of the WHERE2 project DLR measurement
campaign.

.. code:: python

    AnchorNodes = {390:{'name':'MT_ACO_05','coord':[6,0.81,1.64]},
                   386:{'name':'MT_ACO_08','coord':[30.574,2.8,1.291]},
                   391:{'name':'MT_ACO_07','coord':[11.78,-5.553,1.5]},
                   385:{'name': 'MT_ACO_01','coord':[19.52,-0.69,1.446]},
                   387:{'name':'MT_ACO_03','coord':[28.606,-0.74,1.467]},
                   400:{'name':'MT_ACO_02','coord':[30.574,2.8,1.291]},
                   1:{'name':'MT_DLR_RTDSlave','coord':[0.85,0,1.18]}
                  }
.. code:: python

    S.tx.clear()
    S.rx.clear()
    S.tx.filant='def.vsh3'
    S.rx.filant='def.vsh3'
    da ={}
    dm ={}
Vizualization of the simulated scenario

.. code:: python

    fig,ax=S.L.showG('s',nodes=False)
    plt.axis('off')
    #
    # add new points in tx and rx
    #
    #for c,k in enumerate(AnchorNodes):
    c = 0 # first anchor nodes
    k = AnchorNodes.keys()[c]
    pta = array([AnchorNodes[k]['coord'][0],AnchorNodes[k]['coord'][1],AnchorNodes[k]['coord'][2]]).reshape(3,1)
    #
    # To add a point 
    #
    S.tx.point(pta,mode="add")
    da[c]=k
    plt.plot(pta[0,:],pta[1,:],'or')



.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x7fa7ab599dd0>]




.. image:: DLR-WHERE2_files/DLR-WHERE2_14_1.png


In the following a trajectory for the receiver is defined.

``linevect`` function allows to define a linear trajectory from ``ptt``
along direction ``vec``.

.. code:: python

    S.rx.linevect(npt=290, step=0.1, ptt=[0, 0, 1.275], vec=[1, 0, 0], mode='subst')
    ps = S.rx.position[:,-1]
    S.rx.linevect(npt=60, step=0.1, ptt=ps,vec=[0,1,0],mode='append')
Looking what is does

.. code:: python

    S.L.display['ednodes']=False
    S.L.display['edges']=True
    S.L.display['nodes']=False
    S.L.display['title']='Trajectory to be simulated'
    S.show(s=20)

.. parsed-literal::

    Warning : no furniture file loaded




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f71504c4e10>,
     <matplotlib.axes.AxesSubplot at 0x2ac0290>)




.. image:: DLR-WHERE2_files/DLR-WHERE2_18_2.png


Choosing a UWB waveform for the simulation

.. code:: python

    wav = wvf.Waveform(type='W1compensate')
    wav.show()


.. image:: DLR-WHERE2_files/DLR-WHERE2_20_0.png


running the simulation

.. code:: python

    #tcir = evalcir(S,wav,cutoff=4)
Saving the data in pickle format

.. code:: python

    #file = open("tcir5.pickle","w")
    #pickle.dump(tcir,file)
    #file.close()
Reading the data from the above file

.. code:: python

    #del tcir
    file=open("tcir5.pickle","r")
    tcir=pickle.load(file)
    file.close()
    #del ttcir
    #
    for i in tcir[1].keys():
        cir = tcir[1][i]
        cir.zlr(0,150)
        try:
            ttcir=np.vstack((ttcir,cir.y))
        except:
            ttcir=cir.y
.. code:: python

    tcir[1][1].x
    tcir[1][102].x



.. parsed-literal::

    array([  1.01214575e-02,   3.03643725e-02,   5.06072874e-02, ...,
             1.49949393e+02,   1.49969636e+02,   1.49989879e+02])



Aggregated CIR along a synthetic trajectory (line in the corridor)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    plt.figure(figsize=(20,20))
    dmax=150
    plt.imshow(20*np.log10(ttcir+1e-20),vmax=-40,vmin=-120,origin='lower',extent=[0,dmax,1,69],interpolation='nearest')
    plt.xlabel(r'delay $\times$ c (meters)',fontsize=20)
    #plt.ylabel(r'distance along trajectory (meters)',fontsize=20)
    plt.ylabel(r'trajectory index number',fontsize=20)
    clb=plt.colorbar()
    clb.set_label('level (dB)',fontsize=20)
    
    plt.axis('tight')



.. parsed-literal::

    (0.0, 150.0, 1.0, 69.0)




.. image:: DLR-WHERE2_files/DLR-WHERE2_29_1.png


.. code:: python

    tcir[1][10].plot(typ=['v'])



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f7118a6ded0>,
     array([[<matplotlib.axes.AxesSubplot object at 0x7f7118a6de90>]], dtype=object))




.. image:: DLR-WHERE2_files/DLR-WHERE2_30_1.png


.. code:: python

    plt.figure(figsize=(10,5))
    tcir[1][1].plot(typ=['v'])
    xlabel('Delay (ns)')
    ylabel('Level (V)')
    title('Received Waveform')



.. parsed-literal::

    <matplotlib.text.Text at 0x7f7150569810>




.. parsed-literal::

    <matplotlib.figure.Figure at 0x7f7118749410>



.. image:: DLR-WHERE2_files/DLR-WHERE2_31_2.png


.. code:: python

    tcir[1][11].plot(typ=['v'])
    xlabel('Delay (ns)')
    ylabel('Level (V)')
    title('Received Waveform')



.. parsed-literal::

    <matplotlib.text.Text at 0x7f715058d090>




.. image:: DLR-WHERE2_files/DLR-WHERE2_32_1.png

