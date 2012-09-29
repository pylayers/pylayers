.. _routines.Layout:

.. module:: pylayers.gis


gis
===

This modules gathers all the classes and function related to geographical
information.

Class Layout
************

.. currentmodule:: pylayers.gis.layout.Layout


This class handles the description of the indoor propagation environment. 
It manages various input and output data format. 

The most basic description of a Layout is given in a `.str2` file 
which contains the points and edges list of the Layout.

Then this raw information is processes to generate a .str file which 
contains additional visibility relationships between vertices and edges


.. autoclass:: pylayers.gis.layout.Layout

Layout Edition
--------------

.. autosummary::
    :toctree: generated/

    editor 
    add_pnod
    add_fnod
    add_nfpe
    add_none
    add_edge
    add_subseg
    add_window
    add_door
    del_node
    del_edge
    del_cycle
    get_zone
    

Layout Vizualization
--------------------

.. autosummary::
    :toctree: generated/

    showG 
    showGv 
    show3 

    

Layout I/O
----------

.. autosummary::
    :toctree: generated/

    loadfur
    loadstr
    loadstr2
    savestr2
    saveGv 
    savelay 
    loadlay 
    geomfile 
    save


Graph Construction 
------------------

.. autosummary::
    :toctree: generated/

    buildGw 
    buildGv 
    buildGi 
    buildGr 

Signatures 
----------

.. autosummary::
    :toctree: generated/
    
    signature 
    points_image 
    showSig 

Functions related to nodes mobility
-----------------------------------

.. autosummary::
    :toctree: generated/

    waypoint 
    waypointGw 
    thwall 
    
Utilities
---------

.. autosummary::
    :toctree: generated/

    info 
    pt2ro 
    facets3D 
    ispoint 
    onseg 
    facet3D 
    get_Sg_pos 
    distwall 
    randTxRx 
    boundary 
    loadGv 

readvrml
********

.. automodule:: pylayers.gis.readvrml

.. autosummary::
    :toctree: generated/

    savestr2
    stretch
    segsplit
    extract
    inbracket
    incrochet
    geomLine
    geomFace
    ParseDirectionalLight
    ParseMaterial
    show
    parsevrml
    vrml2sha

Class furniture
***************

.. automodule:: pylayers.gis.furniture.Furniture

.. autoclass:: pylayers.gis.furniture.Furniture

.. autosummary::
    :toctree: generated/

    info
    set_position 
    load
    save
    position 
    show


Class selectl
***************

.. automodule:: pylayers.gis.selectl.SelectL

.. autoclass:: pylayers.gis.selectl.SelectL

.. autosummary::
    :toctree: generated/
    
    show
    call_editor 
