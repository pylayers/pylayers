.. _routines.bsignal:

Module Bsignal 
==============


Class Bsignal
-------------

signal with a x axis basis

.. currentmodule:: pylayers.signal.bsignal.Bsignal

.. class:: pylayers.signal.bsignal.Bsignal 

.. autoattribute:: pylayers.signal.bsignal.Bsignal 

.. autosummary::
    :toctree: generated/

    setx
    gating 
    setx
    sety
    stem
    step
    plot
    gating
    len
    save
    load
              
Class Usignal
-------------

Uniform signal 

.. currentmodule:: pylayers.signal.bsignal.Usignal

.. class:: pylayers.signal.bsignal.Usignal

.. autoattribute:: pylayers.signal.bsignal.Usignal 

.. autosummary::
    :toctree: generated/

    setx
    dx
    width
    expand
    min
    max
    truncate
    align
    abs
    energy
    zright
    zleft
    zlr

TBsignal
---------

Time signal 

.. currentmodule:: pylayers.signal.bsignal.TBsignal

.. class:: pylayers.signal.bsignal.TBsignal

.. autoattribute:: pylayers.signal.bsignal.TBsignal 

.. autosummary::
    :toctree: generated/

    b2u
    plot 
    translate
    

TUsignal
---------

Uniform time signal 

.. currentmodule:: pylayers.signal.bsignal.TUsignal

.. class:: pylayers.signal.bsignal.TUsignal

.. autoattribute:: pylayers.signal.bsignal.TUsignal 

.. autosummary::
    :toctree: generated/

    diff
    info
    fft
    fftsh
    filter
    ftshift
    psd
    show
    esd
    shift
    correlate
    corrgauss
    Efirst_loc
    resample
    ft 
    convolve
    Yadd_zero2l
    Yadd_zero2r
    Epercent
    Etau0
    Ewin
    Etot 
    Efirst 
    Efirst_corr
    Efirst_toath
    taumax
    Emax
    tau_Emax
    toa_max2
    toa_new
    toa_win
    toa_max
    toa_th
    toa_cum
    toa_th_tm
    toa_th_tmt
    toa_cum_tm
    toa_cum_tmtm
    toa_cum_tmt
    readuwb
    ecdf
    tau_moy
    tau_rms



Utilities
~~~~~~~~~

.. autosummary::
    :toctree: generated/
    
    info

Visualization
~~~~~~~~~~~~~

.. autosummary::
    :toctree: generated/
    
    show

Processing 
~~~~~~~~~~

.. autosummary::
    :toctree: generated/

    ft
    fft
    fftsh
    convolve
    filter
    ftshift
    psd
    esd
    shift
    correlate
    corrgauss
    resample
    Yadd_zeros2l
    Yadd_zeros2

Determining signal energy content
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
  :toctree: generated/

    Epercent
    Etau0
    Ewin
    Etot
    Efirst
    Efirst_corr
    Efirst_toath
    Efirst_loc

Time of Arrival estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. autosummary::
  :toctree: generated/

    taumax
    tau_Emax
    TDoA
    toa_max2
    toa_new
    toa_win
    toa_max

Class TUDsignal
-----~~~~~~~----

Uniform time signal with delays 

.. currentmodule:: pylayers.signal.bsignal.TUDsignal

.. class:: pylayers.signal.bsignal.TUDsignal

.. autoattribute:: pylayers.signal.bsignal.TUDsignal 

.. autosummary::
    :toctree: generated/

    fig

Class FBsignal
--------------

Frequency signal 

.. currentmodule:: pylayers.signal.bsignal.FBsignal

.. class:: pylayers.signal.bsignal.FBsignal

.. autoattribute:: pylayers.signal.bsignal.FBsignal 

.. autosummary::
    :toctree: generated/

    plotri
    plot
    plotdB
    stem


Class FUsignal
--------------

Uniform frequency signal 

.. currentmodule:: pylayers.signal.bsignal.FUsignal

.. class:: pylayers.signal.bsignal.FUsignal

.. autoattribute:: pylayers.signal.bsignal.FUsignal 

.. autosummary::
    :toctree: generated/

    window
    get
    info
    energy 
    enthrsh
    dBthrsh
    zp
    newdf
    dftresamp
    resample
    symH
    symHz
    align
    ifft
    ift
    iftshift
    show
    decimate


Class FUDsignal
---------------

Uniform frequency signal with delays

.. currentmodule:: pylayers.signal.bsignal.FUDsignal

.. class:: pylayers.signal.bsignal.FUDsignal

.. autoattribute:: pylayers.signal.bsignal.FUDsignal 


.. autosummary::
    :toctree: generated/
    
    minphas
    totud
    iftd
    ft1
    ft2
    ftau

Class FHsignal
--------------

Frequency hermitian signal 

.. currentmodule:: pylayers.signal.bsignal.FHsignal

.. class:: pylayers.signal.bsignal.FHsignal

.. autoattribute:: pylayers.signal.bsignal.FHsignal 

.. autosummary::
    :toctree: generated/

    ifft 
    unrex

Class Noise
------------

.. currentmodule:: pylayers.signal.bsignal.Noise

.. class:: pylayers.signal.bsignal.Noise

.. autoattribute:: pylayers.signal.bsignal.Noise


.. autosummary::
    :toctree: generated/

    amplify 
    gating 

Class EnImpulse
---------------

.. currentmodule:: pylayers.signal.bsignal.EnImpulse

.. class:: pylayers.signal.bsignal.EnImpulse

.. autoattribute:: pylayers.signal.bsignal.EnImpulse


MaskImpulse
-----------

.. currentmodule:: pylayers.signal.bsignal.MaskImpulse



    
