
.. code:: python

    %pylab inline
    import  pylayers.util.mayautil as myu

.. parsed-literal::

    Populating the interactive namespace from numpy and matplotlib


.. parsed-literal::

    WARNING: pylab import has clobbered these variables: ['concatenate', 'random', 'mlab', 'show', 'f', 'dist', 'copy']
    `%matplotlib` prevents importing * from pylab and numpy


Corser
-------

1. `CORMORAN Measurement Campaign <#CORMORAN>`__

   1. `Motivation for creating a specific tool <#Motivation>`__
   2. `Prerequisite installations <#pre>`__

2. `The CorSer Class <#CorSer>`__

   1. `Load Serie <#Load>`__
   2. `Get information on the Series <#getinfo>`__
   3. `Available Data Variable <#Available%20Data>`__

      1. `Radio DataFrames <#Radio%20df>`__
      2. `Non Radio DataFrames <#Non%20Radio%20df>`__
      3. `Involved Devices <#device>`__

   4. `Acessing the Data <#access%20data>`__

      1. `Get Device Position <#getdevp>`__
      2. `Get Link Value <#getlink>`__
      3. `Get Link Distance <#getlinkd>`__

   5. `Visualizing the data <#vizu>`__

      1. `Native Pandas Vizualization <#pandas>`__
      2. `Plot Method <#plot>`__
      3. `Plot Visibility Method <#plotvisi>`__
      4. `Plot Mobility Method <#plotmob>`__
      5. `3D plot <#3Dplot>`__
      6. `3D plot interactively <#3Dploti>`__
      7. `3D Interrative Visibility <#visii>`__

CORMORAN Measurement Campaign 
==============================

The CORMORAN measurement campaign is the first known campaign in the
WBAN context gathering :

-  3 differents radio technologies (HiKoB, CEA plateform and Beespoon
   phone)
-  Up to 24 radio devices equiped on a single body
-  A precise capture of the radio device and body movement using a Vicon
   motion capture (MOCAP) system.
-  A perfect knowledge of the capture environement
-  58 Series with Capture or Group Navigation scenarios

One of the main advantage of this campaign is the use of a precise
motion capture system which allows to get a ground truth position of any
radio device which make the radio observable values open to
interpretation.

Motivation for creating a specific tool 
----------------------------------------

In order to exploit the CORMORAN measurement campaign, a dedicated tool
has been envisaged. Regarding that one aim of the CORMORAN project
(http://pylayers.github.io/pylayers/cormoran.html) is to provide a
simulation plateform from the Channel to the MAC Layer, the tool
natturally takes place inside the PyLayers plateform.

This specific tool creation has been motivated by the intrinsec
complexity of the the measurement campaign. First, no existing tool are
able to exploit simultaneously the radio and MOCAP information from the
measures.

The co-existence of 3 different radio technologies implies 3 different
file formats which have to be interpreted and combined together to be
exploitable. However, the motion capture, and the 3 differents radio
acces technologies (RAT) operating at different sample rate, which leads
to manipulating 4 different time basis. As well, no automatic
start-synchronization mechanism was availble between the different
technologies which leads to a non systematic time shift between the
different basis.

Finally the aim of a measurement campaign is to easily provide valuable
and exploitable information for the project members and more generally
by people in the research community. The goal of such a tool is to help
and simplify dissemination.

Prerequisite Installations 
---------------------------

Before starting using this tool, some requirements have to be satisfied.

1. The open source platform PyLayers ( http://www.pylayers.org ) has to
   be installed following the installation notes here:
   https://github.com/pylayers/pylayers/blob/master/INSTALL.txt

2. The CORMORAN measurements have to be downloaded from the gitlab
   repository (URL PROVIDED SHORTLY)
3. An environement variable $CORMORAN has to be set at the root of your
   CORMORAN measurements directory (help about setup of environement
   variables can be found in pylayers' INSTALL.txt

Once those 3 steps are satisfied, the CORMORAN exploitation measure tool
is ready to be used.

The ``CorSer`` Class 
=====================

The exploitation of measures tool takes place as a specific class named
CorSer (which stands for Cormoran Series). Once PyLayers has been
installed, it is possible to directly access to the class by importing
it.

.. code:: python

    from pylayers.measures.cormoran import *

.. parsed-literal::

    WARNING:traits.has_traits:DEPRECATED: traits.has_traits.wrapped_class, 'the 'implements' class advisor has been deprecated. Use the 'provides' class decorator.


 Get information on the Series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before creating the CorSer object it is possible to consult the
available measurements series using *cor\_log()*. Then for each
**serie** of a given **day** it is possible to get:

-  The involved subject(s)
-  The radio technology
-  A short description of the serie

.. code:: python

    cor_log()



.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>serie</th>
          <th>day</th>
          <th>Subject</th>
          <th>techno</th>
          <th>Short Notes</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0 </th>
          <td>  1</td>
          <td> 11</td>
          <td>            Bernard</td>
          <td>        TCR</td>
          <td>                           Subject Walk circularly</td>
        </tr>
        <tr>
          <th>1 </th>
          <td>  2</td>
          <td> 11</td>
          <td>           Bernard </td>
          <td>        TCR</td>
          <td>                           Subject Walk circularly</td>
        </tr>
        <tr>
          <th>2 </th>
          <td>  3</td>
          <td> 11</td>
          <td>           Bernard </td>
          <td>        TCR</td>
          <td>                           Subject Walk circularly</td>
        </tr>
        <tr>
          <th>3 </th>
          <td>  4</td>
          <td> 11</td>
          <td>           Bernard </td>
          <td>        TCR</td>
          <td>                           Subject Walk circularly</td>
        </tr>
        <tr>
          <th>4 </th>
          <td>  5</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td>                           Subject Walk circularly</td>
        </tr>
        <tr>
          <th>5 </th>
          <td>  6</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td>                           Subject Walk circularly</td>
        </tr>
        <tr>
          <th>6 </th>
          <td>  7</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td>                           Subject Walk circularly</td>
        </tr>
        <tr>
          <th>7 </th>
          <td>  8</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td>                           Subject Walk circularly</td>
        </tr>
        <tr>
          <th>8 </th>
          <td>  9</td>
          <td> 11</td>
          <td>            Bernard</td>
          <td>        TCR</td>
          <td>     INTERRUPTED  Subject Walk circularly ++ speed</td>
        </tr>
        <tr>
          <th>9 </th>
          <td> 10</td>
          <td> 11</td>
          <td>           Bernard </td>
          <td>        TCR</td>
          <td>                  Subject Walk circularly ++ speed</td>
        </tr>
        <tr>
          <th>10</th>
          <td> 11</td>
          <td> 11</td>
          <td>            Bernard</td>
          <td>        TCR</td>
          <td>                  Subject Walk circularly ++ speed</td>
        </tr>
        <tr>
          <th>11</th>
          <td> 12</td>
          <td> 11</td>
          <td>            Bernard</td>
          <td>        TCR</td>
          <td>                  Subject Walk circularly ++ speed</td>
        </tr>
        <tr>
          <th>12</th>
          <td> 13</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td> Subject Walk circularly without looking BS pho...</td>
        </tr>
        <tr>
          <th>13</th>
          <td> 14</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td>     Subject Walk circularly + Navigation movement</td>
        </tr>
        <tr>
          <th>14</th>
          <td> 15</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td> Subject Walk slowly without looking BS phone h...</td>
        </tr>
        <tr>
          <th>15</th>
          <td> 16</td>
          <td> 11</td>
          <td>           Nicolas </td>
          <td>     HKB+BS</td>
          <td> Subject Walk slowly without looking BS phone h...</td>
        </tr>
        <tr>
          <th>16</th>
          <td> 17</td>
          <td> 11</td>
          <td>            Bernard</td>
          <td>        TCR</td>
          <td> Static subject pointing corners then yoga post...</td>
        </tr>
        <tr>
          <th>17</th>
          <td> 18</td>
          <td> 11</td>
          <td>            Bernard</td>
          <td>        TCR</td>
          <td> Static subject pointing corners then yoga post...</td>
        </tr>
        <tr>
          <th>18</th>
          <td> 19</td>
          <td> 11</td>
          <td>            Bernard</td>
          <td>        TCR</td>
          <td> Static subject pointing corners then yoga post...</td>
        </tr>
        <tr>
          <th>19</th>
          <td> 20</td>
          <td> 11</td>
          <td>            Bernard</td>
          <td>        TCR</td>
          <td> Static subject pointing corners then yoga post...</td>
        </tr>
        <tr>
          <th>20</th>
          <td> 21</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td> Static subject pointing corners (withphone) th...</td>
        </tr>
        <tr>
          <th>21</th>
          <td> 22</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td> Static subject pointing corners (withphone) th...</td>
        </tr>
        <tr>
          <th>22</th>
          <td> 23</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td> INTERRUPTED Static subject pointing corners (w...</td>
        </tr>
        <tr>
          <th>23</th>
          <td> 24</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td> Static subject pointing corners (withphone) th...</td>
        </tr>
        <tr>
          <th>24</th>
          <td> 25</td>
          <td> 11</td>
          <td>            Bernard</td>
          <td>        TCR</td>
          <td>                                      Kung-fu Kata</td>
        </tr>
        <tr>
          <th>25</th>
          <td> 26</td>
          <td> 11</td>
          <td>            Bernard</td>
          <td>        TCR</td>
          <td>                     Kung-fu Kata with lost sensor</td>
        </tr>
        <tr>
          <th>26</th>
          <td> 27</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td> subject open door, sit, type on leyboard, take...</td>
        </tr>
        <tr>
          <th>27</th>
          <td> 28</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td> subject open door, sit, type on leyboard, take...</td>
        </tr>
        <tr>
          <th>28</th>
          <td> 29</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td>   Crossfade Yoga Posture with  phone BS left hand</td>
        </tr>
        <tr>
          <th>29</th>
          <td> 30</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td> Crossfade SLOW Yoga Posture with  phone BS lef...</td>
        </tr>
        <tr>
          <th>30</th>
          <td> 31</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td>     HKB+BS</td>
          <td>                           Subject Walk circularly</td>
        </tr>
        <tr>
          <th>31</th>
          <td> 32</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td> TCR+HKB+BS</td>
          <td>       3 turns  circularly inc. speed sequentially</td>
        </tr>
        <tr>
          <th>32</th>
          <td> 33</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td> TCR+HKB+BS</td>
          <td>       3 turns  circularly inc. speed sequentially</td>
        </tr>
        <tr>
          <th>33</th>
          <td> 34</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td> TCR+HKB+BS</td>
          <td> 3 turns  circularly inc. Speed + muscle-buildi...</td>
        </tr>
        <tr>
          <th>34</th>
          <td> 35</td>
          <td> 11</td>
          <td>            Nicolas</td>
          <td> TCR+HKB+BS</td>
          <td> 3 turns  circularly inc. Speed + muscle-buildi...</td>
        </tr>
        <tr>
          <th>35</th>
          <td>  1</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>        TCR</td>
          <td>                         DATA ISSUE 3 FireMen Nav </td>
        </tr>
        <tr>
          <th>36</th>
          <td>  2</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>        TCR</td>
          <td>              3 FireMen Nav (possible mocap issue)</td>
        </tr>
        <tr>
          <th>37</th>
          <td>  3</td>
          <td> 12</td>
          <td> Nicolas jihad Eric</td>
          <td>        TCR</td>
          <td>              3 FireMen Nav (possible mocap issue)</td>
        </tr>
        <tr>
          <th>38</th>
          <td>  4</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>        TCR</td>
          <td>                        INTERRUPTED 3 FireMen Nav </td>
        </tr>
        <tr>
          <th>39</th>
          <td>  5</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>        TCR</td>
          <td> subjects Random walk + new interfering subject...</td>
        </tr>
        <tr>
          <th>40</th>
          <td>  6</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>        TCR</td>
          <td> subjects Random walk + new interfering subject...</td>
        </tr>
        <tr>
          <th>41</th>
          <td>  7</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>        TCR</td>
          <td> subjects slow Random walk + interfering subjec...</td>
        </tr>
        <tr>
          <th>42</th>
          <td>  8</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>        TCR</td>
          <td> subjects slow Random walk + interfering subjec...</td>
        </tr>
        <tr>
          <th>43</th>
          <td>  9</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td> TCR+HKB+BS</td>
          <td> Subject Slow motion: Indoor Nav then Firemen t...</td>
        </tr>
        <tr>
          <th>44</th>
          <td> 10</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td> TCR+HKB+BS</td>
          <td> Subject Slow motion: Indoor Nav then Firemen t...</td>
        </tr>
        <tr>
          <th>45</th>
          <td> 11</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td> TCR+HKB+BS</td>
          <td> Subject normal speed: Indoor Nav then Firemen ...</td>
        </tr>
        <tr>
          <th>46</th>
          <td> 12</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td> TCR+HKB+BS</td>
          <td> Subject normal speed: Indoor Nav then Firemen ...</td>
        </tr>
        <tr>
          <th>47</th>
          <td> 13</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td> TCR+HKB+BS</td>
          <td> subjects Random walk + new interfering subject...</td>
        </tr>
        <tr>
          <th>48</th>
          <td> 14</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td> TCR+HKB+BS</td>
          <td> subjects Random walk + new interfering subject...</td>
        </tr>
        <tr>
          <th>49</th>
          <td> 15</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td> TCR+HKB+BS</td>
          <td> subjects Random walk + new interfering subject...</td>
        </tr>
        <tr>
          <th>50</th>
          <td> 16</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td> TCR+HKB+BS</td>
          <td> subjects Random walk + new interfering subject...</td>
        </tr>
        <tr>
          <th>51</th>
          <td> 17</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>     HKB+BS</td>
          <td> NO HKB Subject normal speed: Indoor Nav then F...</td>
        </tr>
        <tr>
          <th>52</th>
          <td> 18</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>     HKB+BS</td>
          <td> NO HKB Subject normal speed: Indoor Nav then F...</td>
        </tr>
        <tr>
          <th>53</th>
          <td> 19</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>     HKB+BS</td>
          <td> NO HKB Subject normal speed: Indoor Nav then F...</td>
        </tr>
        <tr>
          <th>54</th>
          <td> 20</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>     HKB+BS</td>
          <td> NO HKB Subject normal speed: Indoor Nav then F...</td>
        </tr>
        <tr>
          <th>55</th>
          <td> 21</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>     HKB+BS</td>
          <td> subjects Random walk + new interfering subject...</td>
        </tr>
        <tr>
          <th>56</th>
          <td> 22</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>     HKB+BS</td>
          <td> subjects Random walk + new interfering subject...</td>
        </tr>
        <tr>
          <th>57</th>
          <td> 23</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>     HKB+BS</td>
          <td> subjects Random walk + new interfering subject...</td>
        </tr>
        <tr>
          <th>58</th>
          <td> 24</td>
          <td> 12</td>
          <td> Nicolas Jihad Eric</td>
          <td>     HKB+BS</td>
          <td> subjects Random walk + new interfering subject...</td>
        </tr>
      </tbody>
    </table>
    </div>



 Load Serie
-----------

As an example, serie 6 from day 11 can be loaded using the following
command:

.. code:: python

    S=CorSer(serie=6,day=11)

.. parsed-literal::

    
    load infrastructure node position: **** Processor coding : Intel-PC
    
    load  Nicolas  body: **** Processor coding : Intel-PC
    
    BS data frame index:  Align on mocap OK... WARNING time-offset NOT applied
    No BS offset not yet set => use self.offset_setter 
    
    HKB data frame index: Align on mocap OK... time-offset applied OK
    
    Create distance Dataframe... OK


Once loaded information about the serie (date, type, ...) can be
obtained just by calling the object itself:

.. code:: python

    S



.. parsed-literal::

    Filename: Sc20_S6_R2_HKBS
    Day : 11/06/2014
    Serie : 6
    Scenario : 20
    Run : 2
    Type : HKBS
    Original Video Id : Single
    Subject(s) : Nicolas 
    
    Body available: True
    
    BeSPoon : Sc20_S6_R2_HKBS.csv
    HIKOB : Sc2_0_S6_r2_HKB_Single.mat



 Available data
---------------

 Radio DataFrames
~~~~~~~~~~~~~~~~~

Data frames are *Pandas* objects which can be interpreted as tables.

-  Each line correspond a given timestamp
-  Each column correspond to a given link between 2 radio devices

Depending on available RAT involved in the serie, different data frames
are available:

-  HiKoB (HKB) data : *S.hkb*
-  BeSpoon data : *S.bespo*
-  TCR data : *S.tcr*

In the example serie chosen, only HiKoB and Bespoon are available.

Here is an example of the RSS values obtained by the HKB sensors for the
120 available links and the 5 first available timestamp :

.. code:: python

    S.hkb.head(5)



.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>AP1-AP2</th>
          <th>AP1-AP3</th>
          <th>AP1-AP4</th>
          <th>AP1-HeadRight</th>
          <th>AP1-TorsoTopRight</th>
          <th>AP1-TorsoTopLeft</th>
          <th>AP1-BackCenter</th>
          <th>AP1-ElbowRight</th>
          <th>AP1-ElbowLeft</th>
          <th>AP1-HipRight</th>
          <th>...</th>
          <th>WristRight-WristLeft</th>
          <th>WristRight-KneeLeft</th>
          <th>WristRight-AnkleLeft</th>
          <th>WristRight-AnkleRight</th>
          <th>WristLeft-KneeLeft</th>
          <th>WristLeft-AnkleLeft</th>
          <th>WristLeft-AnkleRight</th>
          <th>KneeLeft-AnkleLeft</th>
          <th>KneeLeft-AnkleRight</th>
          <th>AnkleLeft-AnkleRight</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0.000000</th>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>0.010001</th>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>0.020002</th>
          <td>-60</td>
          <td>-64</td>
          <td>-61</td>
          <td>-71</td>
          <td>-81</td>
          <td>-73</td>
          <td>-78</td>
          <td>-79</td>
          <td>-84</td>
          <td>-73</td>
          <td>...</td>
          <td>-64</td>
          <td>-88</td>
          <td>-64</td>
          <td>-55</td>
          <td>-63</td>
          <td>-61</td>
          <td>-77</td>
          <td>-60</td>
          <td>-84</td>
          <td>-79</td>
        </tr>
        <tr>
          <th>0.030003</th>
          <td>-60</td>
          <td>-64</td>
          <td>-61</td>
          <td>-71</td>
          <td>-81</td>
          <td>-73</td>
          <td>-78</td>
          <td>-79</td>
          <td>-84</td>
          <td>-73</td>
          <td>...</td>
          <td>-64</td>
          <td>-88</td>
          <td>-64</td>
          <td>-55</td>
          <td>-63</td>
          <td>-61</td>
          <td>-77</td>
          <td>-60</td>
          <td>-84</td>
          <td>-79</td>
        </tr>
        <tr>
          <th>0.040004</th>
          <td>-60</td>
          <td>-64</td>
          <td>-61</td>
          <td>-71</td>
          <td>-81</td>
          <td>-73</td>
          <td>-78</td>
          <td>-79</td>
          <td>-84</td>
          <td>-73</td>
          <td>...</td>
          <td>-64</td>
          <td>-88</td>
          <td>-64</td>
          <td>-55</td>
          <td>-63</td>
          <td>-61</td>
          <td>-77</td>
          <td>-60</td>
          <td>-84</td>
          <td>-79</td>
        </tr>
      </tbody>
    </table>
    <p>5 rows × 120 columns</p>
    </div>



 Non Radio DataFrames
~~~~~~~~~~~~~~~~~~~~~

Extra data frames are also available to acces to non radio information.
In particular, it exists :

-  *S.devdf*: the device dataframe, which gives mechanical information:
   position (x,y,z), velocity (v,vx,vy,vz) and acceleration (a,ax,ay,az)
   of the devices at any time stamps
-  *S.distdf*: the distance data frame, which gives ground truth
   distances between the different radio links.

Here is the 5 last data of the device data frame...

.. code:: python

    S.devdf.tail(5)



.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>id</th>
          <th>subject</th>
          <th>x</th>
          <th>y</th>
          <th>z</th>
          <th>v</th>
          <th>vx</th>
          <th>vy</th>
          <th>vz</th>
          <th>a</th>
          <th>ax</th>
          <th>ay</th>
          <th>az</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>104.2</th>
          <td> HKB:14</td>
          <td> Nicolas</td>
          <td> 0.158588</td>
          <td>-1.574102</td>
          <td> 0.526740</td>
          <td> 0.012375</td>
          <td>-0.005046</td>
          <td> 0.010521</td>
          <td> 0.004119</td>
          <td> 2.241849</td>
          <td> 1.972888</td>
          <td> 0.738384</td>
          <td> 0.767065</td>
        </tr>
        <tr>
          <th>104.2</th>
          <td>  HKB:1</td>
          <td>        </td>
          <td> 0.018552</td>
          <td>-2.749937</td>
          <td> 0.979166</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>104.2</th>
          <td> HKB:16</td>
          <td> Nicolas</td>
          <td>-0.229677</td>
          <td>-1.445404</td>
          <td> 0.175125</td>
          <td> 0.010563</td>
          <td>-0.007414</td>
          <td>-0.006640</td>
          <td>-0.003540</td>
          <td> 0.547761</td>
          <td> 0.122199</td>
          <td>-0.250196</td>
          <td>-0.471711</td>
        </tr>
        <tr>
          <th>104.2</th>
          <td> HKB:10</td>
          <td> Nicolas</td>
          <td> 0.262695</td>
          <td>-1.433168</td>
          <td> 1.143153</td>
          <td> 0.057829</td>
          <td>-0.048329</td>
          <td>-0.030039</td>
          <td>-0.010302</td>
          <td> 0.924303</td>
          <td>-0.697193</td>
          <td> 0.368582</td>
          <td>-0.482085</td>
        </tr>
        <tr>
          <th>104.2</th>
          <td>  HKB:3</td>
          <td>        </td>
          <td> 0.021135</td>
          <td> 3.375590</td>
          <td> 1.003871</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
        </tr>
      </tbody>
    </table>
    </div>



... and the 5 last data of the distance data frame:

.. code:: python

    S.distdf.tail(5)



.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>HKB:1-HKB:2</th>
          <th>HKB:1-HKB:3</th>
          <th>HKB:1-HKB:4</th>
          <th>HKB:1-HKB:5</th>
          <th>HKB:1-HKB:6</th>
          <th>HKB:1-HKB:7</th>
          <th>HKB:1-HKB:8</th>
          <th>HKB:1-HKB:9</th>
          <th>HKB:1-HKB:10</th>
          <th>HKB:1-HKB:11</th>
          <th>...</th>
          <th>HKB:12-HKB:15</th>
          <th>HKB:12-HKB:16</th>
          <th>HKB:13-HKB:14</th>
          <th>HKB:13-HKB:15</th>
          <th>HKB:13-HKB:16</th>
          <th>HKB:14-HKB:15</th>
          <th>HKB:14-HKB:16</th>
          <th>HKB:15-HKB:16</th>
          <th>BS:0-BS:74</th>
          <th>BS:0-BS:157</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>104.159996</th>
          <td> 6.102589</td>
          <td> 6.125578</td>
          <td> 6.135849</td>
          <td> 1.308815</td>
          <td> 1.163639</td>
          <td> 1.131707</td>
          <td> 1.387571</td>
          <td> 1.322510</td>
          <td> 1.350930</td>
          <td> 1.223406</td>
          <td>...</td>
          <td> 1.071233</td>
          <td> 0.990922</td>
          <td> 0.411064</td>
          <td> 0.753501</td>
          <td> 0.910143</td>
          <td> 0.364396</td>
          <td> 0.539795</td>
          <td> 0.445009</td>
          <td> 1.046829</td>
          <td> 0.119864</td>
        </tr>
        <tr>
          <th>104.169997</th>
          <td> 6.102589</td>
          <td> 6.125578</td>
          <td> 6.135849</td>
          <td> 1.309074</td>
          <td> 1.163713</td>
          <td> 1.131587</td>
          <td> 1.387549</td>
          <td> 1.322884</td>
          <td> 1.350486</td>
          <td> 1.223658</td>
          <td>...</td>
          <td> 1.071489</td>
          <td> 0.990873</td>
          <td> 0.410944</td>
          <td> 0.753502</td>
          <td> 0.909901</td>
          <td> 0.364396</td>
          <td> 0.539682</td>
          <td> 0.445027</td>
          <td> 1.046903</td>
          <td> 0.119868</td>
        </tr>
        <tr>
          <th>104.179998</th>
          <td> 6.102589</td>
          <td> 6.125578</td>
          <td> 6.135849</td>
          <td> 1.309470</td>
          <td> 1.163938</td>
          <td> 1.131414</td>
          <td> 1.387530</td>
          <td> 1.323230</td>
          <td> 1.350018</td>
          <td> 1.223874</td>
          <td>...</td>
          <td> 1.071624</td>
          <td> 0.990832</td>
          <td> 0.410933</td>
          <td> 0.753522</td>
          <td> 0.909759</td>
          <td> 0.364316</td>
          <td> 0.539533</td>
          <td> 0.445038</td>
          <td> 1.046936</td>
          <td> 0.119734</td>
        </tr>
        <tr>
          <th>104.189999</th>
          <td> 6.102589</td>
          <td> 6.125578</td>
          <td> 6.135849</td>
          <td> 1.309873</td>
          <td> 1.164064</td>
          <td> 1.131319</td>
          <td> 1.387509</td>
          <td> 1.323601</td>
          <td> 1.349608</td>
          <td> 1.224129</td>
          <td>...</td>
          <td> 1.071955</td>
          <td> 0.990734</td>
          <td> 0.410871</td>
          <td> 0.753529</td>
          <td> 0.909520</td>
          <td> 0.364281</td>
          <td> 0.539368</td>
          <td> 0.445063</td>
          <td> 1.047000</td>
          <td> 0.119982</td>
        </tr>
        <tr>
          <th>104.200000</th>
          <td> 6.102589</td>
          <td> 6.125578</td>
          <td> 6.135849</td>
          <td> 1.310357</td>
          <td> 1.164289</td>
          <td> 1.131228</td>
          <td> 1.387509</td>
          <td> 1.323915</td>
          <td> 1.349214</td>
          <td> 1.224341</td>
          <td>...</td>
          <td> 1.072294</td>
          <td> 0.990736</td>
          <td> 0.410651</td>
          <td> 0.753482</td>
          <td> 0.909291</td>
          <td> 0.364271</td>
          <td> 0.539394</td>
          <td> 0.445110</td>
          <td> 1.046967</td>
          <td> 0.119830</td>
        </tr>
      </tbody>
    </table>
    <p>5 rows × 122 columns</p>
    </div>



 Involved devices (*S.dev*)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The *S.dev* command allows to obtain the complete list of devices
involved in the serie and:

-  the Name of the device used in the radio dataframe
-  the Real device Id used during the measurement campaign
-  The corresponding device Id used on the Body wear description
-  At wich Subject the device is related.

Infrastrucure access point obviously don't have related Subject.

.. code:: python

    S.dev

.. parsed-literal::

    Name in Dataframe     | Real Id | Body Id  | Subject    
    ========================================================
    AP4                   |       4 | HKB:4    |            
    AP1                   |       1 | HKB:1    |            
    AP2                   |       2 | HKB:2    |            
    AP3                   |       3 | HKB:3    |            
    --------------------------------------------------------          
    AnkleRight            |      16 | HKB:16   | Nicolas    
    KneeLeft              |      14 | HKB:14   | Nicolas    
    AnkleLeft             |      15 | HKB:15   | Nicolas    
    WristRight            |      12 | HKB:12   | Nicolas    
    WristLeft             |      13 | HKB:13   | Nicolas    
    ElbowLeft             |      10 | HKB:10   | Nicolas    
    HipRight              |      11 | HKB:11   | Nicolas    
    HeadRight             |       5 | HKB:5    | Nicolas    
    TorsoTopRight         |       6 | HKB:6    | Nicolas    
    TorsoTopLeft          |       7 | HKB:7    | Nicolas    
    BackCenter            |       8 | HKB:8    | Nicolas    
    ElbowRight            |       9 | HKB:9    | Nicolas    
                          |         |          |            
    WristRight            |     157 | BS:157   | Nicolas    
    AnkleRight            |      74 | BS:74    | Nicolas    
    HandRight             |       0 | BS:0     | Nicolas    
    --------------------------------------------------------          


 Accessing the data
-------------------

In order to help people not familiar with the Pandas query format, some
useful methods are provided in order to extract values from radio and
non radio dataframes.

 ### Get device position (*S.getdevp*)

The value of the device position at a specific time or range or time can
be obtained by specifying:

-  The device (Name in dataframe OR real id OR body id)
-  The radio *techno* (Precising the techno is optional except when an
   ambiguity occurs, therefore error is raised)
-  a given time in second or a [start time,stop time]. If no time is
   given, the position for all time stamps are provided

Hence, It is possible to get the positions of the HKB radio node 11 (Hip
Right), between 5.0 seconds and 5.2 seconds with:

.. code:: python

    Positions = S.getdevp(11,t=[5,5.2])
    Positions



.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>x</th>
          <th>y</th>
          <th>z</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>5.000480</th>
          <td>-0.139566</td>
          <td> 0.224905</td>
          <td> 1.016796</td>
        </tr>
        <tr>
          <th>5.010481</th>
          <td>-0.139553</td>
          <td> 0.224845</td>
          <td> 1.016826</td>
        </tr>
        <tr>
          <th>5.020482</th>
          <td>-0.139545</td>
          <td> 0.224825</td>
          <td> 1.016818</td>
        </tr>
        <tr>
          <th>5.030483</th>
          <td>-0.139564</td>
          <td> 0.224730</td>
          <td> 1.016849</td>
        </tr>
        <tr>
          <th>5.040484</th>
          <td>-0.139609</td>
          <td> 0.224642</td>
          <td> 1.016859</td>
        </tr>
        <tr>
          <th>5.050485</th>
          <td>-0.139580</td>
          <td> 0.224613</td>
          <td> 1.016898</td>
        </tr>
        <tr>
          <th>5.060486</th>
          <td>-0.139554</td>
          <td> 0.224586</td>
          <td> 1.016920</td>
        </tr>
        <tr>
          <th>5.070487</th>
          <td>-0.139604</td>
          <td> 0.224492</td>
          <td> 1.016937</td>
        </tr>
        <tr>
          <th>5.080488</th>
          <td>-0.139545</td>
          <td> 0.224452</td>
          <td> 1.016989</td>
        </tr>
        <tr>
          <th>5.090489</th>
          <td>-0.139521</td>
          <td> 0.224391</td>
          <td> 1.016992</td>
        </tr>
        <tr>
          <th>5.100489</th>
          <td>-0.139386</td>
          <td> 0.224397</td>
          <td> 1.016997</td>
        </tr>
        <tr>
          <th>5.110490</th>
          <td>-0.139296</td>
          <td> 0.224315</td>
          <td> 1.017041</td>
        </tr>
        <tr>
          <th>5.120491</th>
          <td>-0.139164</td>
          <td> 0.224189</td>
          <td> 1.017098</td>
        </tr>
        <tr>
          <th>5.130492</th>
          <td>-0.138988</td>
          <td> 0.224128</td>
          <td> 1.017131</td>
        </tr>
        <tr>
          <th>5.140493</th>
          <td>-0.138810</td>
          <td> 0.224048</td>
          <td> 1.017142</td>
        </tr>
        <tr>
          <th>5.150494</th>
          <td>-0.138605</td>
          <td> 0.223969</td>
          <td> 1.017148</td>
        </tr>
        <tr>
          <th>5.160495</th>
          <td>-0.138406</td>
          <td> 0.223877</td>
          <td> 1.017164</td>
        </tr>
        <tr>
          <th>5.170496</th>
          <td>-0.138043</td>
          <td> 0.223803</td>
          <td> 1.017230</td>
        </tr>
        <tr>
          <th>5.180497</th>
          <td>-0.137791</td>
          <td> 0.223654</td>
          <td> 1.017305</td>
        </tr>
        <tr>
          <th>5.190498</th>
          <td>-0.137388</td>
          <td> 0.223580</td>
          <td> 1.017321</td>
        </tr>
      </tbody>
    </table>
    </div>



**NOTE : You may also obtain a classical numpy array instead of this
Pandas object by using the "*values*\ " method :**

.. code:: python

    Positions.values



.. parsed-literal::

    array([[-0.13956557,  0.22490462,  1.01679608],
           [-0.13955284,  0.22484492,  1.01682581],
           [-0.13954524,  0.22482529,  1.01681787],
           [-0.1395645 ,  0.2247298 ,  1.01684918],
           [-0.13960907,  0.224642  ,  1.01685901],
           [-0.13957962,  0.2246127 ,  1.01689801],
           [-0.13955351,  0.22458575,  1.01691986],
           [-0.13960399,  0.22449205,  1.01693719],
           [-0.13954485,  0.22445244,  1.01698865],
           [-0.13952087,  0.22439058,  1.0169917 ],
           [-0.13938625,  0.22439655,  1.0169975 ],
           [-0.13929645,  0.22431535,  1.01704102],
           [-0.13916449,  0.22418907,  1.0170979 ],
           [-0.1389884 ,  0.22412761,  1.01713135],
           [-0.13880983,  0.22404759,  1.0171424 ],
           [-0.13860497,  0.22396939,  1.01714777],
           [-0.1384055 ,  0.22387668,  1.01716443],
           [-0.13804305,  0.22380293,  1.01722955],
           [-0.13779123,  0.2236543 ,  1.01730511],
           [-0.13738791,  0.22358025,  1.01732141]])



 Get link value (*S.getlink*)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The value of a link *a* and *b* at a specific time or range or time can
be obtained by specifying:

-  The device :math:`a` (Name in dataframe OR real id OR body id)
-  The device :math:`b` (Name in dataframe OR real id OR body id)
-  The radio *technoa* and *technob* (Precising the techno is optional
   except when an ambiguity occurs, therefore error is raised)
-  a given time in second or a [start time,stop time]. If no time is
   given, the position for all time stamps are provided

Hence, It is possible to get the HKB values between radio node 11 (Hip
Right) and node 16 (Ankle Right) , between 5 seconds and 5.2 seconds
with:

.. code:: python

    Values = S.getlink(11,16,t=[5,5.2])
    Values



.. parsed-literal::

    5.000500   -67
    5.010501   -67
    5.020502   -67
    5.030503   -67
    5.040504   -67
    5.050505   -67
    5.060506   -67
    5.070507   -67
    5.080508   -67
    5.090509   -67
    5.100510   -67
    5.110511   -67
    5.120512   -67
    5.130513   -67
    5.140514   -67
    5.150515   -67
    5.160516   -67
    5.170517   -67
    5.180518   -67
    5.190519   -66
    Name: HipRight-AnkleRight, dtype: float64



 Get link distance (*S.getlinkd*)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ground truth distance separating a device *a* and device *b* at a
specific time or range or time can be obtained by specifying:

-  The device :math:`a` (Name in dataframe OR real id OR body id)
-  The device :math:`b` (Name in dataframe OR real id OR body id)
-  The radio *technoa* and *technob* (Precising the techno is optional
   except when an ambiguity occurs, therefore error is raised)
-  a given time in second or a [start time,stop time]. If no time is
   given, the position for all time stamps are provided

Hence, It is possible to get the HKB values between radio node 11 (Hip
Right) and node 16 (Ankle Right) , between 5 seconds and 5.2 seconds
with:

.. code:: python

    Distances = S.getlinkd(11,16,t=[5,5.2])
    Distances



.. parsed-literal::

    5.000480    0.845013
    5.010481    0.845034
    5.020482    0.845045
    5.030483    0.845068
    5.040484    0.845090
    5.050485    0.845180
    5.060486    0.845229
    5.070487    0.845235
    5.080488    0.845309
    5.090489    0.845339
    5.100489    0.845353
    5.110490    0.845423
    5.120491    0.845482
    5.130492    0.845559
    5.140493    0.845563
    5.150494    0.845595
    5.160495    0.845602
    5.170496    0.845677
    5.180497    0.845769
    5.190498    0.845785
    Name: HKB:11-HKB:16, dtype: float64



 Visualizing the Data
---------------------

 Native Pandas Vizualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because radio data in CorSer are stored into Pandas objects, convenient
vizualization method are directly available. Most of them can be found
here : http://pandas.pydata.org/pandas-docs/stable/visualization.html

As an example, it is possbile tthe previous obtained values and distance
with :

.. code:: python

    # Ploting 
    ax=Values.plot() # plot values
    l=Distances.plot(secondary_y=True,ax=ax) # plot distances on the right side
    
    ## Labelling 
    ax.legend() # add legend box
    ax.set_ylabel('RSS Values (dBm)') # set left ylabel
    ax.right_ax.set_ylabel('Distances (m)') # set right ylabel
    ax.set_xlabel('time (s)') # set xlabel
    ax.set_title('RSS and distance as a function of time')



.. parsed-literal::

    <matplotlib.text.Text at 0x7f39e3a54f50>




.. image:: Corser_files/Corser_45_1.png


In addition, CorSer also provides specific plotting methods which
includes extra features.

 Plot method (S.plot)
~~~~~~~~~~~~~~~~~~~~~

The plot function allows to display the radio values of a link. The main
parameters are always the same:

-  The device :math:`a` (Name in dataframe OR real id OR body id)
-  The device :math:`b` (Name in dataframe OR real id OR body id)
-  The radio *techno* (Precising the techno is optional except when an
   ambiguity occurs, therefore error is raised)
-  A given time in second or a [start time,stop time]. If no time is
   given, the position for all time stamps are provided

More option are availble, please refer to the docstring (*S.plot?*) for
more information

 Plot values
^^^^^^^^^^^^

Continuying with the same example, it is possible to plot the HKB values
between radio node 11 (Hip Right) and node 16 (Ankle Right) , between 5
seconds and 5.2 seconds with:

.. code:: python

    S.plot(11,16,t=[5,5.2])



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f39e85e54d0>,
     <matplotlib.axes.AxesSubplot at 0x7f39e3aa97d0>)




.. image:: Corser_files/Corser_49_1.png


Plot distance
~~~~~~~~~~~~~

As well, it is possible to plot the distance using the *distance*
parameter

.. code:: python

    S.plot(11,16,t=[5,5.2],distance = True)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f39e3884650>,
     <matplotlib.axes.AxesSubplot at 0x7f39e38e6e90>)




.. image:: Corser_files/Corser_52_1.png


It is also possible to get the same result than with the Pandas
procedure with the following code :

.. code:: python

    # plot value
    f,ax = S.plot(11,16,t=[5,5.2],color ='b',title=False)
    
    # create right axis
    ax2=ax.twinx()
    
    # plot distance
    S.plot(11,16,t=[5,5.2],color ='g',title=False,
           distance=True,
           fig=f,ax=ax2)




.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f39e3884850>,
     <matplotlib.axes.AxesSubplot at 0x7f39e3799410>)




.. image:: Corser_files/Corser_54_1.png


 Plot visibility method (S.pltvisi)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to go further in the radio value interpretation, it is
convenient to have some extra information about the **optical
visibility/occultation** of devices involved in a link.

This information allows to determine the line of sight (LOS) or non line
of sight (NLOS) cases which are crutial for power level and delay
interpretation.

This information can be superimposed to the radio values. To this end,
the plot visibility (*S.pltvisi*) method is used. The **hatched** area
denoted **NLOS** wheras **clear** area denotes **LOS**.

Parameters are the same than those the *plot* method:

.. code:: python

    f,ax = S.plot(1,16)
    S.pltvisi(1,16,fig=f,ax=ax)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f39e3799f90>,
     <matplotlib.axes.AxesSubplot at 0x7f39e368edd0>)




.. image:: Corser_files/Corser_57_1.png


 Plot mobility method (S.pltmob)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As well it is possible to determine and indicate whether the subject is
static or not by using the plot mobility method (*S.pltmob*). The
succession of Static and Mobile sequences are denoted :math:`S_x` and
:math:`M_x` resplectively, where :math:`x` is an index of the sequence.

.. code:: python

    f,ax = S.plot(1,16)
    S.pltmob(fig=f,ax=ax)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f194f27e2d0>,
     <matplotlib.axes.AxesSubplot at 0x7f194f021310>)




.. image:: Corser_files/Corser_60_1.png


The 2 upmentionned methods can also be used simultaneously as shown in
the following example :

.. code:: python

    # plot data in green)
    f,ax=S.plthkb(1,13,figsize=(10,5))
    # plot optical occultation (hatched lines)
    S.pltvisi(1,13,fig=f,ax=ax)
    # plot subject mobility (grey areas)
    S.pltmob(showvel=False,ylim=([-100,-40]),fig=f,ax=ax)



.. parsed-literal::

    (<matplotlib.figure.Figure at 0x7f194f273bd0>,
     <matplotlib.axes.AxesSubplot at 0x7f194f35ff50>)




.. image:: Corser_files/Corser_62_1.png


 3D plot (S.\_show3)
~~~~~~~~~~~~~~~~~~~~

With the help of the Mayavi Library, the CorSer class allows to display
in 3D :

-  The building where measurements have taken place
-  The positions of Vicon Cameras
-  The Multi-cylindric representation of the the subjects involved in
   the selected serie
-  The position/ antenna pattern of the devices on the body(ies) and in
   the infrastructure.

By default, the use of the \*S.\_show3\* method display the complete
scene with body(ies) and associated devices at 4 different timestamp

.. code:: python

    S._show3()
    
    # the following line is only used to display in the notebook a screenshot of the mayavi window
    myu.inotshow('fig1')

.. parsed-literal::

    /home/uguen/anaconda/lib/python2.7/site-packages/traits/has_traits.py:1766: FutureWarning: comparison to `None` will result in an elementwise object comparison in the future.
      setattr( self, name, value )



.. image:: Corser_files/Corser_65_1.png


Specify time (*bodytime* parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to display scene at specific timestamps, the parameter
*bodytime* can be used

Example: to show the body position at :math:`t=0s`, :math:`t=30s` and
:math:`t=90s`.

.. code:: python

    S._show3(bodytime=[0.,30.,90.])
    
    # the following line is only used to display in the notebook a screenshot of the mayavi window
    myu.inotshow('fig2')


.. image:: Corser_files/Corser_68_0.png


display trajectory (*trajectory* parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    S._show3(trajectory = True,bodytime=[0.,30.,90.])
    
    # the following line is only used to display in the notebook a screenshot of the mayavi window
    myu.inotshow('fig3')


.. image:: Corser_files/Corser_70_0.png


 3D plot interactive (\*S.\_show3i\*)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The method \*S.\_show3i()\* allows to display the 3D scene with an extra
window incluying a slider acting like a jog shuttle, to choose the
timestamp to vizualize.

Note : This function is note available in the notebook

.. code:: python

    S._show3i(t=35) # t=35 is an initialization value
    
    # the following line is only used to display in the notebook a screenshot of the mayavi window
    myu.inotshow('fig4')


.. image:: Corser_files/Corser_73_0.png



.. image:: Corser_files/Corser_73_1.png


 Interactive visibility (*S.imshowvisibility\_i*)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The visibility matrix can be displayed simultaneously to the 3D view.

For that purpose a visibility/occultation matrix is computed the first
time the vizualization is called. The following code displays the
tisibility matrix and associated 3D scene at the inital time
:math:`t=35s`

.. code:: python

    S.imshowvisibility_i(t=35)
    
    # the following line is only used to display in the notebook a screenshot of the mayavi window
    myu.inotshow('fig5')

.. parsed-literal::

    Visibility is computed only once, Please wait
    
    processing shadowing from 

.. parsed-literal::

    /home/uguen/anaconda/lib/python2.7/site-packages/traits/has_traits.py:1766: FutureWarning: comparison to `None` will result in an elementwise object comparison in the future.
      setattr( self, name, value )
    /home/uguen/anaconda/lib/python2.7/site-packages/traits/has_traits.py:1771: FutureWarning: comparison to `None` will result in an elementwise object comparison in the future.
      setattr( self, name, value )


.. parsed-literal::

     Nicolas



.. image:: Corser_files/Corser_76_3.png



.. image:: Corser_files/Corser_76_4.png


Using Pylayers Ray-tracing with CorSer data
-------------------------------------------

Coming soon, work in progress

.. code:: python

    import pylayers.simul.simultraj as st
.. code:: python

    #ST=st.Simul(S)
