.. module:: pylayers.antprop


antprop
=======

This modules gather all classes and functions related to antenna and
propagation modelling.

Module antenna
**************

.. currentmodule:: pylayers.antprop.antenna


Functions
---------

.. autosummary::
    :toctree: generated/

    indexvsh
    geom_pattern
    VW
    AFLegendre
    show3D

Class SHCoeff
-------------

.. currentmodule:: pylayers.antprop.antenna.SHCoeff

.. class:: pylayers.signal.antenna.SHCoeff

.. autoattribute:: pylayers.signal.antenna.SHCoeff

.. autosummary::
    :toctree: generated/

    inits1
    inits2
    inits3
    s1tos2
    s3tos2
    delete
    delete3
    put3 
    info
    show

Class VSHCoeff
---------------

.. currentmodule:: pylayers.antprop.antenna.VSHCoeff

.. class:: pylayers.signal.antenna.VSHCoeff

.. autoattribute:: pylayers.signal.antenna.VSHCoeff

.. autosummary::
    :toctree: generated/

    show
    s1tos2
    s2tos3
    strip3
    ens3 
    drag3 
    put3

Class Antenna
-------------

.. currentmodule:: pylayers.antprop.antenna.Antenna

.. class:: pylayers.antprop.antenna.Antenna 

.. autoattribute:: pylayers.antprop.antenna.Antenna

.. autosummary::
    :toctree: generated/

    loadmat 
    load_trx 
    loadtrx 
    checkpole 
    info 
    help 
    polar 
    show3_geom 
    show3 
    pol3d 
    mse 
    elec_delay 
    vshd 
    vsh 
    demo 
    Fsynth 
    Fsynth2 
    Fsynth3 
    movie_vsh 
    minsh3 
    savevsh3 
    loadvsh3 
    savevsh2 
    loadvsh2 
    loadvsh3_old 
    pol2cart 
    cart2pol 

    

Module Slab
***********

    
.. currentmodule:: pylayers.antprop.slab.Interface

.. class:: pylayers.antprop.slab.Interface

.. autoattribute:: pylayers.antprop.slab.Interface

.. autosummary::
    :toctree: generated/

    RT
    loss0
    losst 
    pcolor
    plotwrta
    plotwrtf 

.. currentmodule:: pylayers.antprop.slab.MatInterface

.. class:: pylayers.antprop.slab.MatInterface

.. autoattribute:: pylayers.antprop.slab.MatInterface

.. autosummary::
    :toctree: generated/


.. currentmodule:: pylayers.antprop.slab.Mat

.. class:: pylayers.antprop.slab.Mat

.. autoattribute:: pylayers.antprop.slab.Mat

.. autosummary::
    :toctree: generated/

    R

.. currentmodule:: pylayers.antprop.slab.MatDB

.. class:: pylayers.antprop.slab.MatDB

.. autoattribute:: pylayers.antprop.slab.MatDB

.. autosummary::
    :toctree: generated/

    info 
    dass
    maxindex
    delete
    edit 

.. currentmodule:: pylayers.antprop.slab.Slab

.. class:: pylayers.antprop.slab.Slab

.. autoattribute:: pylayers.antprop.slab.Slab

.. autosummary::
    :toctree: generated/

    conv 
    editgui
    info 
    loss0 
    losst


.. currentmodule:: pylayers.antprop.slab.SlabDB

.. class:: pylayers.antprop.slab.SlabDB

.. autoattribute:: pylayers.antprop.slab.SlabDB

.. autosummary::
    :toctree: generated/

    addgui
    choose
    dass
    delete
    edit 
    help
    info
    load
    save
    loadsl
    savesl
    maxindex 
    show
    showall



Channel
*******


Class Ctilde
-------------

.. currentmodule:: pylayers.antprop.channel.Ctilde

.. class:: pylayers.antprop.channel.Ctilde

.. autoattribute:: pylayers.antprop.channel.Ctilde

.. autosummary::
    :toctree: generated/

    choose
    info 
    doadod
    energy
    sort 
    vec2scal 
    vec2scalA


Class VectChannel 
-----------------

.. currentmodule:: pylayers.antprop.channel.VectChannel

.. class:: pylayers.antprop.channel.VectChannel

.. autoattribute:: pylayers.antprop.channel.VectChannel

.. autosummary::
    :toctree: generated/

    show3
    mobility


Class ScalChannel 
-----------------

.. currentmodule:: pylayers.antprop.channel.ScalChannel

.. class:: pylayers.antprop.channel.ScalChannel

.. autosummary::
    :toctree: generated/

    info
    imshow
    apply 
    applywavA
    applywavB
    applywavC
    doddoa
    wavefig
    rayfig
