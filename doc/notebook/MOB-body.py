{
 "metadata": {
  "name": "",
  "signature": "sha256:52c820e7f4e0f4b4da7bd2ec17815679f601d9f37d383be1a526069e44a9ef73"
 },
 "nbformat": 3,
 "nbformat_minor": 0,
 "worksheets": [
  {
   "cells": [
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "# Handling Body Mobility in PyLayers"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "from pylayers.mobility.body.body import *\n",
      "from pylayers.mobility.trajectory import Trajectory\n",
      "from IPython.display import Image"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "display_data",
       "text": [
        "<matplotlib.figure.Figure at 0x894a310>"
       ]
      }
     ],
     "prompt_number": 1
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "The body mobility is imported from motion capture files. This is the chosen manner to achieve a high degree of realism for the modeling of the human motion. Two kind of files exist : \n",
      "\n",
      "+ `c3d` files are a set of point which are evolving in time\n",
      "+ `bvh` files are a stuctured version of the motion capture.\n",
      "\n",
      "Both type of file will be exploited in the following. "
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "## `BodyCylinder` data structure"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "To ease electromagnetic simulation a simplification of the motion capture data structure is necessary. Generally there is a large number of captutred points, not all of them being useful for our modeling. \n",
      "\n",
      "The body model is a restriction of key body segments which are transformed into $K$ cylinders of radius $r_k$. \n",
      "\n",
      "The chosen body model is made of 11 cylinders. 4 cylinders decribing the two arms, 4 cylinders decribing the two legs,  2 cylinders describing the trunk and 1 cylinder for the head. \n",
      "\n",
      "The body cylinder model is handle by a dedicated Python class call `Body`"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "To create a void body, simply instantiate a Body object from the class"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John = Body()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 2
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "which is equivalent to :"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John = Body(_filebody='John.ini',_filemocap='07_01.c3d')"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 3
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "The default body filename is John.ini and the default motion capture filename is '07_01.c3d'. The creation of a Body consists in reading a `_filebody` and a `_filemocap`"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "## Description of a body file\n",
      "\n",
      "An example of a body file is given below. It is a file in `ini` format with 4 sections. \n",
      "\n",
      "+ [nodes] \n",
      " \n",
      "   This section associates a node number to a c3d fils conventional  node number \n",
      "   \n",
      "    NodeId = C3DNODE\n",
      "    \n",
      "+ [cylinder]\n",
      "\n",
      "   This section associates a cylinder Id to a dictionnary wich contains cylinder tail head and radius information \n",
      "   \n",
      "      CylId = {'t',NodeId1,'h',NodeId2,'r',float (m),'name',}\n",
      "      \n",
      "+ [device]\n",
      "\n",
      "   This section associates a device name to a dictionnary wich contains cylinder device related information\n",
      "   \n",
      "       DevId = {'typ' : {static|mobile} \n",
      "                 'cyl': CylId \n",
      "                 'l'  : length coordinate in ccs,\n",
      "                 'h'  : height coordinate in ccs, \n",
      "                 'a'  : angle coordinate in ccs, \n",
      "                 'file' : antenna file , \n",
      "                 'T' : Rotation matrix }\n",
      "   "
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "#%load /home/uguen/Bureau/P1/ini/Francois.ini"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 4
    },
    {
     "cell_type": "raw",
     "metadata": {},
     "source": [
      "[nodes]\n",
      "0 = STRN\n",
      "1 = CLAV\n",
      "2 = RFHD \n",
      "3 =RSHO \n",
      "4 =LSHO\n",
      "5 =RELB \n",
      "6 =LELB \n",
      "7 =RWRB \n",
      "8 =LWRB\n",
      "9 =RFWT\n",
      "10 =LFWT \n",
      "11 =RKNE\n",
      "12 =LKNE \n",
      "13 =RANK \n",
      "14 =LANK \n",
      "15 =BOTT \n",
      "[cylinder]\n",
      "; sternum (STRN) - clavicule (CLAV)\n",
      "trunku = {'t':0,'h':1,'r':0.18,'i':0}\n",
      "; bottom  (BOTT) sternum (STRN)  \n",
      "trunkb = {'t':15,'h':0,'r':0.17,'i':10}\n",
      "; clavicule (CLAV)  - tete (RFHD)\n",
      "headu = {'t':1,'h':2,'r':0.12,'i':1}\n",
      "; coude droit (RELB)  epaule droite (RSHO)\n",
      "armr = {'t':5,'h':3,'r':0.05,'i':2}\n",
      "; coude gauche (LELB)  epaule gauche (LSHO)\n",
      "arml  = {'t':6,'h':4,'r':0.05,'i':3}\n",
      "; poignet droit (RWRB) coude droit (RELB)  \n",
      "forearmr = {'t':7,'h':5,'r':0.05,'i':4}\n",
      "; left wrist (LWRB)  left elbow (LELB)  \n",
      "forearml = {'t':8,'h':6,'r':0.05,'i':5}\n",
      "; knee droit (RKNE) hanche droit (RFWT)  \n",
      "thighr = {'t':11,'h':9,'r':0.05,'i':6}\n",
      "; knee left (LKNE)  hanche left (LFWT)  \n",
      "thighl = {'t':12,'h':10,'r':0.05,'i':7}\n",
      "; cheville droit (RANK) genou  droit (RKNE)  \n",
      "calfr = {'t':13,'h':11,'r':0.05,'i':8}\n",
      "; cheville droit (LANK) genou  droit (LKNE)  \n",
      "calfl = {'t':14,'h':12,'r':0.05,'i':9}\n",
      "[device]\n",
      "0 = {'typ':'static','name':'BeSpoon Phone','cyl':'trunku','l':0.1,'h':0.01,'a':0,'file':'S2R2.sh3','T':np.array([[1,0,0],[0,1,0],[0,0,1]])}\n",
      "1 = {'typ':'static','name':'Movea Accel','cyl':'trunku','l':0.1,'h':0.01,'a':180,'file':'S2R2.sh3','T':np.array([[1,0,0],[0,1,0],[0,0,1]])}\n",
      "2 = {'typ':'static','name':'Optivent Glass','cyl':'head','l':0.7,'h':0.01,'a':0,'file':'S2R2.sh3','T':np.array([[1,0,0],[0,1,0],[0,0,1]])}\n",
      "3 = {'typ':'static','name':'Geonaute Podo','cyl':'trunkb','l':0.1,'h':0.01,'a':45,'file':'S2R2.sh3','T':np.array([[1,0,0],[0,1,0],[0,0,1]])}\n",
      "4 = {'typ':'static','name':'Breizh Watch','cyl':'forearmr','l':0.2,'h':0.01,'a':0,'file':'S2R2.sh3','T':np.array([[1,0,0],[0,1,0],[0,0,1]])}\n",
      "5 = {'typ':'static','name':'Breizh Watch','cyl':'forearml','l':0.2,'h':0.01,'a':0,'file':'S2R2.sh3','T':np.array([[1,0,0],[0,1,0],[0,0,1]])}\n",
      "[mocap]\n",
      "\n",
      "walk = '07_01_c3d' \n"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 5,
       "text": [
        "My name is : John\n",
        "\n",
        "I have a Galaxy Gear device on the left forearm\n",
        "I have a cardio device on the upper part of trunk\n",
        "I am nowhere yet\n",
        "\n",
        "filename : 07_01.c3d\n",
        "nframes : 300\n",
        "Centered : True\n",
        "Mocap Speed : 1.36558346484"
       ]
      }
     ],
     "prompt_number": 5
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "Francois = Body(_filebody='Francois.ini',_filemocap='07_01.c3d')\n",
      "Francois"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 6,
       "text": [
        "My name is : Francois\n",
        "\n",
        "I have a Movea Accel device on the upper part of trunk\n",
        "I have a BeSpoon Phone device on the upper part of trunk\n",
        "I have a Geonaute Podo device on the lower part of trunk\n",
        "I have a Optivent Glass device hea\n",
        "I have a Breizh Watch device on the left forearm\n",
        "I have a Breizh Watch device on the right forearm\n",
        "I am nowhere yet\n",
        "\n",
        "filename : 07_01.c3d\n",
        "nframes : 300\n",
        "Centered : True\n",
        "Mocap Speed : 1.36558346484"
       ]
      }
     ],
     "prompt_number": 6
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "## Loading a Motion Capture File"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "A `.c3d` motion capture file is loaded with the method `loadC3D` with as arguments the motion capture file and the number of frames to load.  \n",
      "\n",
      "The motion is represented as a sequence of framef stored in the `d` variable member."
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "It is possible to get the information from the C3D header by using the verbose option of the `read_c3d` function"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "# Video Frame Rate \n",
      "Vrate = 120\n",
      "# Inter Frame\n",
      "Tframe = 1./120\n",
      "# select a number of frame \n",
      "nframes = 300\n",
      "# Time duration of the whole selected frame sequence\n",
      "Tfseq = Tframe*nframes\n",
      "#\n",
      "# load a .c3dmotion capture file\n",
      "# this update the g.pos \n",
      "#\n",
      "#bc.loadC3D(filename='07_01.c3d',nframes=nframes,centered=True)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 7
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "The duration of the capture is "
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "print \"Duration of the motion capture sequence\", Tfseq,\" seconds\""
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "Duration of the motion capture sequence 2.5  seconds\n"
       ]
      }
     ],
     "prompt_number": 8
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "`d` is a MDA of shape `(3,npoint,nframe)`. It contains all the possible configurations of the body. In general it is supposed to be a cyclic motion as an integer number of walking steps. This allows to instantiate the body configuration anywhere else in space in a given trajectory. \n",
      "\n",
      "A specific space-time configuration of the body is called a **`topos`**."
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "np.shape(John.d)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 9,
       "text": [
        "(3, 16, 300)"
       ]
      }
     ],
     "prompt_number": 9
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "## Defining a trajectory\n",
      "\n",
      "A Trajectory is a class which : \n",
      "\n",
      "+ derives from a pandas `DataFrame` \n",
      "+ is container for time,position,velocity and acceleration. \n",
      "\n",
      "To define a default trajectory : "
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "traj = Trajectory()\n",
      "t = traj.generate()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 10
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "traj.head()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "html": [
        "<div style=\"max-height:1000px;max-width:1500px;overflow:auto;\">\n",
        "<table border=\"1\" class=\"dataframe\">\n",
        "  <thead>\n",
        "    <tr style=\"text-align: right;\">\n",
        "      <th></th>\n",
        "      <th>x</th>\n",
        "      <th>y</th>\n",
        "      <th>z</th>\n",
        "      <th>vx</th>\n",
        "      <th>vy</th>\n",
        "      <th>vz</th>\n",
        "      <th>ax</th>\n",
        "      <th>ay</th>\n",
        "      <th>az</th>\n",
        "      <th>s</th>\n",
        "    </tr>\n",
        "  </thead>\n",
        "  <tbody>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00</th>\n",
        "      <td> 0.000000</td>\n",
        "      <td> 0.000000</td>\n",
        "      <td>-0.675722</td>\n",
        "      <td> 0.061186</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td> 0.539726</td>\n",
        "      <td>-0.000229</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.401367</td>\n",
        "      <td> 0.000000</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.204082</th>\n",
        "      <td> 0.061186</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td>-0.135996</td>\n",
        "      <td> 0.060957</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td> 0.941093</td>\n",
        "      <td>-0.000458</td>\n",
        "      <td> 0</td>\n",
        "      <td>-0.976458</td>\n",
        "      <td> 0.580256</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.408163</th>\n",
        "      <td> 0.122143</td>\n",
        "      <td> 0.408163</td>\n",
        "      <td> 0.805096</td>\n",
        "      <td> 0.060499</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td>-0.035365</td>\n",
        "      <td>-0.000684</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.439891</td>\n",
        "      <td> 1.545150</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.612245</th>\n",
        "      <td> 0.182642</td>\n",
        "      <td> 0.612245</td>\n",
        "      <td> 0.769731</td>\n",
        "      <td> 0.059815</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td> 0.404526</td>\n",
        "      <td>-0.000909</td>\n",
        "      <td> 0</td>\n",
        "      <td>-2.415204</td>\n",
        "      <td> 1.760928</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.816327</th>\n",
        "      <td> 0.242457</td>\n",
        "      <td> 0.816327</td>\n",
        "      <td> 1.174257</td>\n",
        "      <td> 0.058906</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td>-2.010679</td>\n",
        "      <td>-0.001129</td>\n",
        "      <td> 0</td>\n",
        "      <td> 2.788602</td>\n",
        "      <td> 2.217949</td>\n",
        "    </tr>\n",
        "  </tbody>\n",
        "</table>\n",
        "<p>5 rows \u00d7 10 columns</p>\n",
        "</div>"
       ],
       "output_type": "pyout",
       "prompt_number": 11,
       "text": [
        "                                   x         y         z        vx        vy  \\\n",
        "1970-01-01 00:00:00         0.000000  0.000000 -0.675722  0.061186  0.204082   \n",
        "1970-01-01 00:00:00.204082  0.061186  0.204082 -0.135996  0.060957  0.204082   \n",
        "1970-01-01 00:00:00.408163  0.122143  0.408163  0.805096  0.060499  0.204082   \n",
        "1970-01-01 00:00:00.612245  0.182642  0.612245  0.769731  0.059815  0.204082   \n",
        "1970-01-01 00:00:00.816327  0.242457  0.816327  1.174257  0.058906  0.204082   \n",
        "\n",
        "                                  vz        ax  ay        az         s  \n",
        "1970-01-01 00:00:00         0.539726 -0.000229   0  0.401367  0.000000  \n",
        "1970-01-01 00:00:00.204082  0.941093 -0.000458   0 -0.976458  0.580256  \n",
        "1970-01-01 00:00:00.408163 -0.035365 -0.000684   0  0.439891  1.545150  \n",
        "1970-01-01 00:00:00.612245  0.404526 -0.000909   0 -2.415204  1.760928  \n",
        "1970-01-01 00:00:00.816327 -2.010679 -0.001129   0  2.788602  2.217949  \n",
        "\n",
        "[5 rows x 10 columns]"
       ]
      }
     ],
     "prompt_number": 11
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "f,a = traj.plot()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "display_data",
       "png": "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\n",
       "text": [
        "<matplotlib.figure.Figure at 0x9b51fd0>"
       ]
      }
     ],
     "prompt_number": 12
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "### settopos () method\n",
      "\n",
      "Once the trajectory has been defined it is possible to send the body at the position corresponding to any time \n",
      "of the trajectory with the `settopos` method. \n",
      "\n",
      "settopos takes as argument \n",
      "\n",
      "+ A trajectory \n",
      "+ A time index"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "traj.__repr__()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 13,
       "text": [
        "'t (s) : 0.0:9.591837\\nd (m) : 59.4653280432\\nVmoy (m/s) : 6.19957658196\\n'"
       ]
      }
     ],
     "prompt_number": 13
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.settopos(traj,t=5)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 14
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "figure(figsize=(15,20))\n",
      "for t in arange(traj.tmin+0.4,traj.tmax,0.5):\n",
      "    John.settopos(traj,t=t)\n",
      "    f,a=John.show(color='b',plane='yz',topos=True)\n",
      "    axis('off')"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "display_data",
       "png": "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\nAABJRU5ErkJggg==\n",
       "text": [
        "<matplotlib.figure.Figure at 0x9b5ee10>"
       ]
      }
     ],
     "prompt_number": 15
    },
    {
     "cell_type": "code",
     "collapsed": true,
     "input": [
      "John"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 16,
       "text": [
        "My centroid position is \n",
        "[ 0.31956097  9.3877551 ]\n",
        "filename : 07_01.c3d\n",
        "nframes : 300\n",
        "Centered : True\n",
        "Mocap Speed : 1.36558346484"
       ]
      }
     ],
     "prompt_number": 16
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "Francois.settopos(traj,t=6)\n",
      "Francois"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 17,
       "text": [
        "My centroid position is \n",
        "[ 0.97911919  5.91836735]\n",
        "filename : 07_01.c3d\n",
        "nframes : 300\n",
        "Centered : True\n",
        "Mocap Speed : 1.36558346484"
       ]
      }
     ],
     "prompt_number": 17
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "+ 3  : dimension of space\n",
      "+ 16 : number of nodes\n",
      "+ 300 : number of frames "
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "The figure below shows the projection in a vertival plane of the body nodes. "
     ]
    },
    {
     "cell_type": "heading",
     "level": 2,
     "metadata": {},
     "source": [
      "Centering the motion"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.centered"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 18,
       "text": [
        "True"
       ]
      }
     ],
     "prompt_number": 18
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "In order to translate the motion in any point in space-time, a distinction is made between the real motion or topos and the centered motion capture which acts as a virtual motion. \n",
      "\n",
      "Let $\\mathbf{p}^k$ denotes the center of gravity of the body in the (O,x,y) plane"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.center()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 19
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "a = np.hstack((John.vg,John.vg[:,-1][:,newaxis]))"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 20
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "$\\mathbf{v}_g$ is the velocity vector of the gravity center of the body."
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "print np.shape(John.pg)\n",
      "print np.shape(John.vg)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "(3, 300)\n",
        "(3, 300)\n"
       ]
      }
     ],
     "prompt_number": 21
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "print John.vg[:,145]\n",
      "print John.vg[:,298]"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "[ 0.0114987  -0.00026335  0.        ]\n",
        "[  1.08123514e-02   7.24411022e-05   0.00000000e+00]\n"
       ]
      }
     ],
     "prompt_number": 22
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "At that point the body structure is centered. \n",
      "\n",
      "\n",
      "\n",
      "The frame is centered in the xy plane by substracting from the configuration of points the projection of the body in the xy plane. "
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "np.shape(John.d)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 23,
       "text": [
        "(3, 16, 300)"
       ]
      }
     ],
     "prompt_number": 23
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.npoints"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 24,
       "text": [
        "16"
       ]
      }
     ],
     "prompt_number": 24
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "Each frame is centered above the origin. For example for a walk motion the effect of the centering is just like if the body was still walking but not moving forward exactly in the same manner as a walk on a conveyor belt."
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "pgc = np.sum(John.d[:,:,0],axis=1)/16\n",
      "pg0 = John.pg[:,0]\n",
      "print \"True center of gravity\", pg0\n",
      "print \"Center of gravity of the centered frame\",pgc"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "True center of gravity [-1.74251571  0.49373077  0.        ]\n",
        "Center of gravity of the centered frame [  3.74700271e-16   1.38777878e-17   8.94887363e-01]\n"
       ]
      }
     ],
     "prompt_number": 25
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "np.shape(John.pg)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 26,
       "text": [
        "(3, 300)"
       ]
      }
     ],
     "prompt_number": 26
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "The current file contains 300 frames"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "tframe = arange(John.nframes)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 27
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "np.shape(John.pg[0:-1,:])"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 28,
       "text": [
        "(2, 300)"
       ]
      }
     ],
     "prompt_number": 28
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "xg = John.pg[0,:]\n",
      "yg = John.pg[1,:]\n",
      "zg = John.pg[2,:]\n",
      "figure(figsize=(8,8))\n",
      "subplot(311)\n",
      "plot(tframe,xg)\n",
      "title('x component')\n",
      "ylabel('m')\n",
      "subplot(312)\n",
      "xlabel('frame index')\n",
      "title('y component')\n",
      "ylabel('m')\n",
      "plot(tframe,yg)\n",
      "subplot(313)\n",
      "xlabel('frame index')\n",
      "title('Motion capture centroid trajectory')\n",
      "ylabel('m')\n",
      "plot(xg,yg,'.b')\n",
      "\n",
      "d = John.pg[0:-1,1:]-John.pg[0:-1,0:-1]\n",
      "smocap = cumsum(sqrt(sum(d*d,axis=0)))\n",
      "\n",
      "Vmocap = smocap[-1]/Tfseq\n",
      "title('Length = '+str(smocap[-1])+' V = '+str(Vmocap*3.6)+' km/h')\n",
      "axis('scaled')\n",
      "axis('off')\n",
      "plt.tight_layout()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "display_data",
       "png": "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\neeaZZ2T27NkiIjJixAg5ePCgiIhs2bJF3n77bREROX78uLi7u0tiYqL4+vpKSUmJ8RuPiOocD1cR\nkVFoNBr069cPrVu3VoanTp2KHTt2oLS0FBcvXlTK+vv7w8XFBQDg4eGBAQMGQKvVomfPnli4cCEA\nYMOGDdi4cSNiYmIAALdu3UJRURGeeOKJSl9/wIABaNmyJVq2bIkmTZogOzsb27dvR1hYGOzt7QEA\nzz//PAD9IbX4+Hj8+uuv5Zbl4eGBYcOGISIiAikpKfVy6I2IjI//uURkNA4ODsrzNWvW4MaNG0hO\nTkZBQQHatGmjTGvZsqXyvFmzZspw06ZNUVxcDAAoLS1FXFwcHB0dq/36Dy63qKgIGo0GUsFP9JWV\nlUGr1SIlJQXm5ublpqenp8PGxgbZ2dnVfn0ialh4Tg4R1YmsrCw4OTnB3NwcixcvRllZWY3mf/31\n1zF//nwl9FR0Ts79KgoyGo0GL774IrZt24bs7GxcvHgRO3bsAKAPQeHh4ViwYAFKS0shIjh69CgA\nfS/SrVu3sHv3bowfPx65ubk1qjsRNQwMOURUaw9eqXT/cFRUFJKTk+Hl5YW7d+/C0tKy0vnuH39v\n2rhx42BtbY3evXujc+fOWLRoUZXl739+v1atWmH69Ol46aWX8Nprr+HFF19Upk2fPh1Xr17Fs88+\nC09PT2zatAk3b97E5MmT8c0338DV1RXjxo2r8WXqRNQwaKSirz9EREREjRx7coiIiEiVGHKIiIhI\nlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiV\nGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUY\ncoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhy\niIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKI\niIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiI\niEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiI\nSJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhI\nlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiVGHKIiIhIlRhyiIiISJUYcoiIiEiV\nGHKIiIhIlRhyiIiISJUYckg1zp07B61Wi7KyMlNXhYjqSExMDIYPH14nyw4ODsa3335bJ8sm02DI\nIaNwdnbGjh076v01d+7cWa+vCQAffvghHB0dYW1tjeeeew7/+Mc/qjXfiBEjoNVqkZGRoYxbu3Yt\nevbsCQsLC4SEhBiUv3nzJnr16gVbW1u0adMGAwcORFJSkjJdp9Nh4cKF6NGjB9q1a4dJkyZBRJTp\nwcHBaN68OaysrGBlZQV3d3eD5e/YsQN/+tOfYG9vj0GDBmHXrl0G08+cOYOwsDDY2trC3t4e8+bN\nM5i+bNkydOnSBVZWVvDw8MDp06cBAFu2bEHv3r1hY2ODZ599FrGxsSgqKqrWOgPA6NGj4ebmhiZN\nmmD58uWVtueMGTMQFBRUbvyNGzfQrFkznDhxotJ5a+vu3btwd3dHu3btqiyXkJCAAQMGwM7ODoMG\nDUJOTo4ybdeuXQgJCUHLli3Rvn37CuevrG2vXr2K/v3745lnnoFWq8WFCxcM5tPpdIiOjoaDgwM6\nduyIb775xmD66dOnERgYCBsbGwQFBeHMmTMVvn7fvn3LfWGwtLRUtiUrKyuYmZkhOjoaAHDixAk8\n++yzeOqpp9C2bVtERkbi6NGj1W6/kJAQ2NnZoVWrVggLC8P69esra1poNJpKp9XG3bt30bp1axQU\nFECj0Rh9+WRaDDlkFKbYOWg0GoMP9foycuRInDhxArm5uVi6dCm+/vprJCQkVDlPcnIyMjIyyrVR\nq1at8MEHH2DSpEnl5rG0tMSSJUtw7do1nD9/Hv369UNUVBTu3r0LANiwYQMWL16MZcuWISEhAQkJ\nCfjss8+U+TUaDb788kvk5eUhLy8PJ0+eVKaVlpZixIgR6NevHzIzMxEWFoYRI0YoH2p3795F7969\n0b17dxw7dgxnz55FaGioMv+GDRsQGxuLyZMnIycnB1u2bIGtrS0A4Pbt25g6dSquXLmC77//HnFx\ncVi2bFm11hkAfHx88NVXX6Fr165VblPDhw/Hvn37cO7cOYPxa9asQZcuXeDh4VHpvLU1a9Ys2NnZ\n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\np09LYGCgWFtbS2BgoJw5c6bCdazof0lEZMyYMcrVWvdbt26duLm5iaWlpfj7+8uMGTPk6tWr1Wq/\nkydPSvfu3cXKyko6dOggf/nLX+TIkSMVziuif6+HDx+uDH/yySfi6+srf/zxR7lpiYmJBoGppKRE\ntFqtZGVlyebNmyUiIsJg2Q/uL5YtWyYBAQGV1oUaPo2ICW40QlQLEydORHFxcbmb0hER1dS7774L\nLy8vvP3226auCtUhnpNDDdbVq1exd+9e6HQ6bN68GWvXrjW4IR0RUW35+PjglVdeMXU1qI6xJ4ca\nrAsXLuDll19GRkYGunfvjlGjRuHVV1+FmZmZqatGRESNAEMOERERqRK/EhMRkYHRo4HffwfOngVK\nSoC7d4HmzYHSUsPnOTlAixb6x4PT6vN5dethjPp26wasWwf859Y71MCxJ4foAaNHA5s3A3/8UfkO\nsT53qnW9s2/fHnjySWD1au64G5t722pxsXG3pfx84PZtU69dw2Vurv/fqY/w9eC0Rx1+3EIbQw6p\nzsN2/A/7pxcBbt0y9VrUv+bN9WGnOm3FYFQ3ahpa6npbNTMD7v0UmpUVkJdX/nlV0+rzubHLVfbc\nwgIoKKi63RqbIUOAtWtNXYu6wZBDDV5Nu84B4+74Tb1Tra+dfW1UFIyM1ePU0IJUTQJIbXvwHiW0\nGHtbsrEBkpKAqVMBjQaYOxeYMMHweUkJ0KxZxdPq83l162GM+ubnA9u3V789jVHO2MP3P/f1BXbu\nbBj/Y3WBIYfqXUUfFlV1sT5K13ltdiheXkC7dvrnle0Q63OnWtc7+x49gKtX9YHiXjvXZOdb15o2\nBVq3Nl53f22DSX338FXng/HetmrsbWziRGDRIvV+8D2KW7eAP/+5/sLXg9MedfjBOixdqu73mSHn\nEVW3l6Ghfjt9mOoEkpp+2AC1+7Cobtf5w3b8j/s//YNu3dK/z7NmVb+tHhaMjP0Nt6Gpqx68moaW\nx21bJaophpxHFBwM7N5d8/nudfPXxYlpxjqhrT6+vVa3i7UmXefc8de9qoKRsb7h3gtSgPG7+x+l\nh68ue/C47RIZF0POIwoPBxISav6NtjEyxgcKUPGHRVVdrOw6fzwZ67CAsQ4tMoAQNT4MOY+opt9o\nK+rmB+r2m2pNpj04XJ1AUtMPG35YEBFRfWDIISIiIlXiD3QSERGRKjHkEBERkSox5BAREZEqMeQQ\nERGRKjHkEBERkSox5BAREZEqMeQQERGRKjHkEBERkSox5BAREZEqMeQQERGRKjHkEBERkSox5BAR\nEZEqMeQQERGRKjHkEBERkSox5BAREZEqMeQQERGRKjHkEBERkSox5BAREZEqMeQQERGRKjHkEBER\nkSox5BAREZEqMeQQERGRKjHkEBERkSox5BAREZEqMeQQERGRKjHkEBERkSox5BAREZEqMeQQERGR\nKjHkEBERkSox5BAREZEqMeQQERGRKjHkEBERkSox5BAREZEqMeQQERGRKv0/8xhbQAFnR1QAAAAA\nSUVORK5CYII=\n",
       "text": [
        "<matplotlib.figure.Figure at 0x7fdd44b96410>"
       ]
      }
     ],
     "prompt_number": 29
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "plot(smocap)\n",
      "title('evolution of curvilinear abscisse from motion capture centroid trajectory') \n",
      "xlabel('frame index')\n",
      "ylabel('distance (meters)')"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 30,
       "text": [
        "<matplotlib.text.Text at 0x7fdd45827110>"
       ]
      },
      {
       "output_type": "display_data",
       "png": "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\nGjVqFI0aNYqIiK5cuUIdO3YUd9g3btyge/fuEZHwo8Dd3Z2OHTtGxcXFdObMGSooKFBYps+fPycr\nKysxmaWnp9PVq1eJSNhG586dS0+fPlVYR0RETk5OtGzZMsrOzqYDBw6QiYmJuCMrq2w376JFi0hf\nX5++/fZbevbsGT179ozGjx9P7733HqWkpNDu3bvJ1NSUbt++rbBuNm7cSJmZmTR06FDq3LkzBQUF\nUWZmJv373/+mPn36KJ33i/Xs5OREcXFxdPHiRWrSpAm5ubnRvn37KC0tjTw8POjXX38lIqK8vDxq\n3Lgx/fTTT3T//n0KDAwkLy8vhc8qvV2X3T6r046KtrGy7O3tac6cOZSamkpPnjwRE351trdu3brR\npUuX6NixY9SmTRs6ePAgERH98ssv5brVvL29af369eJjZcmw9PabmZlJJiYmFB4eTjk5ObRs2TLq\n3r17tdqRm5srtmPv3r1kb29PZ86cobS0NHrvvfdowYIFCsvVz8+PkpOTKSIigurVq0fPnj0jImE/\nMW7cOIX5lN03vNi2K9rflW73/v37qW3btnTixAm6ePEieXh4UJ06dSrsJq1oOxg2bBglJiaK3/PK\nvtelvxdVLdOy+8uTJ09WuD4/++wz6t27N928eZOOHj1Kbdq0ocjISDEefX19mj9/Pj169IiePXtG\n//rXv+ijjz4S37969WqxmKhIpcnw8ePH1KJFC3r69Kn43KxZs+irr74iIqKff/5Z/MKWlJRQq1at\n6MSJE0REFBgYSPPnzxffd/jwYbK3t1dYyNVJhsr6uUtv2LGxsdSiRQuFuFu2bEmxsbFEJGwcLyo0\nIqHiKF0NlLZq1Sry9/cXH9+6dYuaNm1aYZylPXr0qMJfHlVVht7e3jRp0iSF94wdO5aWLFlCRETx\n8fFkaGgofmkmTpxYYWVYNhna2dmJr505c4YsLCyISPm6DQoKEtdtWZ988gnNmDFDfCyTycSdbXXc\nunWLTExMxDZPmDCBevbsKb5+8+ZN8fWLFy+SlZUVyeXycl9cd3d3+v3338t9fullWlhYSC1atKAt\nW7aIvREvzJkzh8aOHUtJSUkKz8fHx1PHjh0Vnhs2bBj99ttvSttTeh0QCcu6bdu24uOioiJq2rQp\n3bx5U3xuzJgxtGrVKnH6F98HIqItW7aQnp6e2Gtx9+5dql+/foU/NNq0aSN+FhFR//79afjw4eLj\n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\nTkSuQmCC\n",
       "text": [
        "<matplotlib.figure.Figure at 0x7fdd45337210>"
       ]
      }
     ],
     "prompt_number": 30
    },
    {
     "cell_type": "heading",
     "level": 2,
     "metadata": {},
     "source": [
      "Defining a large scale trajectory "
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "A large scale trajectory is defined in the $(O,x,y)$ plane.\n",
      "\n",
      "`traj` is a data structure (Npt,2)"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "v = Vmocap\n",
      "print v*3.6,\"Kmph\""
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "4.91610047343 Kmph\n"
       ]
      }
     ],
     "prompt_number": 31
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "# time in seconds\n",
      "time = np.arange(0,10,0.01)\n",
      "x = v*time\n",
      "y = np.zeros(len(time))\n",
      "z = np.zeros(len(time))\n",
      "traj = Trajectory()\n",
      "traj.generate(time,np.vstack((x,y,y)).T)\n",
      "traj.tmax"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 32,
       "text": [
        "9.97"
       ]
      }
     ],
     "prompt_number": 32
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "fig ,ax = traj.plot()\n",
      "traj.head()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "html": [
        "<div style=\"max-height:1000px;max-width:1500px;overflow:auto;\">\n",
        "<table border=\"1\" class=\"dataframe\">\n",
        "  <thead>\n",
        "    <tr style=\"text-align: right;\">\n",
        "      <th></th>\n",
        "      <th>x</th>\n",
        "      <th>y</th>\n",
        "      <th>z</th>\n",
        "      <th>vx</th>\n",
        "      <th>vy</th>\n",
        "      <th>vz</th>\n",
        "      <th>ax</th>\n",
        "      <th>ay</th>\n",
        "      <th>az</th>\n",
        "      <th>s</th>\n",
        "    </tr>\n",
        "  </thead>\n",
        "  <tbody>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00</th>\n",
        "      <td> 0.000000</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.013656</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.000000e+00</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.000000</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.010000</th>\n",
        "      <td> 0.013656</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.013656</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.000000e+00</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.013656</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.020000</th>\n",
        "      <td> 0.027312</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.013656</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.000000e+00</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.027312</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.030000</th>\n",
        "      <td> 0.040968</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.013656</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 6.938894e-18</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.040968</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.040000</th>\n",
        "      <td> 0.054623</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.013656</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td>-1.387779e-17</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.054623</td>\n",
        "    </tr>\n",
        "  </tbody>\n",
        "</table>\n",
        "<p>5 rows \u00d7 10 columns</p>\n",
        "</div>"
       ],
       "output_type": "pyout",
       "prompt_number": 33,
       "text": [
        "                                   x  y  z        vx  vy  vz            ax  \\\n",
        "1970-01-01 00:00:00         0.000000  0  0  0.013656   0   0  0.000000e+00   \n",
        "1970-01-01 00:00:00.010000  0.013656  0  0  0.013656   0   0  0.000000e+00   \n",
        "1970-01-01 00:00:00.020000  0.027312  0  0  0.013656   0   0  0.000000e+00   \n",
        "1970-01-01 00:00:00.030000  0.040968  0  0  0.013656   0   0  6.938894e-18   \n",
        "1970-01-01 00:00:00.040000  0.054623  0  0  0.013656   0   0 -1.387779e-17   \n",
        "\n",
        "                            ay  az         s  \n",
        "1970-01-01 00:00:00          0   0  0.000000  \n",
        "1970-01-01 00:00:00.010000   0   0  0.013656  \n",
        "1970-01-01 00:00:00.020000   0   0  0.027312  \n",
        "1970-01-01 00:00:00.030000   0   0  0.040968  \n",
        "1970-01-01 00:00:00.040000   0   0  0.054623  \n",
        "\n",
        "[5 rows x 10 columns]"
       ]
      },
      {
       "output_type": "display_data",
       "png": "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\n",
       "text": [
        "<matplotlib.figure.Figure at 0x9b3f2d0>"
       ]
      }
     ],
     "prompt_number": 33
    },
    {
     "cell_type": "heading",
     "level": 2,
     "metadata": {},
     "source": [
      "Trajectory"
     ]
    },
    {
     "cell_type": "heading",
     "level": 2,
     "metadata": {},
     "source": [
      "`posvel()`"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "\n",
      "\n",
      "\n",
      "The `posvel()` function (position and velocity) takes as arguments the following parameters \n",
      "\n",
      "+ `traj` a plane trajectory object.\n",
      "+ $t_k$ time for evaluation of topos \n",
      "+ $T_{fs}$ duration of the periodic motion frame sequence \n",
      "\n",
      "and returns \n",
      "\n",
      "+ the frame index $k_f = \\lfloor \\frac{t_k \\pmod{T_{fs}}}{t_f} \\rfloor$\n",
      "+ the trajectory index $k_t = \\lfloor t_k \\rfloor$\n",
      "+ velocity unitary vector along motion capture frame $\\hat{\\mathbf{v}}_s = \\frac{\\mathbf{p}^g[k_f]-\\mathbf{p}^g[k_f-1]}{|\\mathbf{p}^g[k_f]-\\mathbf{p}^g[k_f-1]|}$\n",
      "+ $\\hat{\\mathbf{w}}_s = \\mathbf{\\hat{z}} \\times  \\hat{\\mathbf{v}}_s $\n",
      "+  velocity unitary vector along trajectory $\\hat{\\mathbf{v}}_t = \\frac{\\mathbf{p}^t[k_t]-\\mathbf{p}^g[k_t-1]}{|\\mathbf{p}^g[k_t]-\\mathbf{p}^t[k_t-1]|}$\n",
      "+ $\\hat{\\mathbf{w}}_t = \\mathbf{\\hat{z}} \\times  \\hat{\\mathbf{v}}_t $\n",
      "\n",
      "$t_f = \\frac{T_{fs}}{Nf}$ is the interframe time or frame sampling period, it is equal to the whole duration of the motion sequence $T_{fs}$ divided by the number of frames\n",
      "\n",
      "\n",
      "\n"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "`settopos` is a method which takes as argument : \n",
      "\n",
      "+ `traj` a plane trajectory (Npt,2)\n",
      "+ $t_k$ time for evaluation of topos "
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "In futher version of the project, this function will be modified to be able to avoid passing the whole trajectory. "
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.settopos(traj=traj,t=3)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 34
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "There is now a new data structure in the Body objet. This data structure is called a `topos`."
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "print np.shape(John.topos)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "(3, 16)\n"
       ]
      }
     ],
     "prompt_number": 35
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.topos"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 36,
       "text": [
        "array([[ 4.19599817,  4.15046572,  4.13081653,  4.07294142,  4.08289287,\n",
        "         4.11542251,  3.94530449,  4.26514518,  3.82039411,  4.12746908,\n",
        "         4.18189918,  3.89824308,  4.28535797,  3.70884964,  4.41212222,\n",
        "         4.15468413],\n",
        "       [ 0.00978643,  0.01039192, -0.08586912, -0.15014439,  0.17454155,\n",
        "        -0.29891008,  0.28336831, -0.31726715,  0.37578539, -0.14822518,\n",
        "         0.12069953, -0.11848735,  0.13268648, -0.04485156,  0.07025805,\n",
        "        -0.01376282],\n",
        "       [ 1.18925612,  1.35442836,  1.57836099,  1.39037819,  1.40039528,\n",
        "         1.07349379,  1.11418117,  0.83620759,  0.89701301,  0.89713179,\n",
        "         0.89777028,  0.44315412,  0.43969984,  0.09704225,  0.05562085,\n",
        "         0.89745103]])"
       ]
      }
     ],
     "prompt_number": 36
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.show3()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 37
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.settopos(traj=traj,t=1)\n",
      "fig,ax=John.plot3d(topos=True)\n",
      "John.settopos(traj=traj,t=4)\n",
      "John.plot3d(topos=True,fig=fig,ax=ax)\n",
      "John.settopos(traj=traj,t=6)\n",
      "John.plot3d(topos=True,fig=fig,ax=ax)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 38,
       "text": [
        "(<matplotlib.figure.Figure at 0x7fdd458256d0>,\n",
        " <matplotlib.axes.Axes3DSubplot at 0x7fdd4582d2d0>)"
       ]
      },
      {
       "output_type": "display_data",
       "png": "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\n",
       "text": [
        "<matplotlib.figure.Figure at 0x7fdd458256d0>"
       ]
      }
     ],
     "prompt_number": 38
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "## Definition of Several Coordinates systems"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "Each cylinder of the `Body` model bears one specific coordinate system. \n",
      "\n",
      "One or several cylinder coordinate systems can be chosen to define the Body Local Coordinates System (BLCS) which is required for motion capture (BLCS) applications. \n",
      "\n",
      "In general, the origin will be chosen on a position which is the most time invariant as on the chest or the back. \n",
      "\n",
      "\n",
      "Those frames of references are all defined in the Global Coordinate System (GCS) of the scene."
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "### Construction of the Cylinder Coordinate System (CCS)"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "The method `setccs()` is used to associate a Cylinder Coordinate System (CCS) to each cylinder of the bodyCylinder model. Notice that those cylinders coordinates systems are not known by the localization application. The localization application will define the BLCS from the position of radiating devices placed on the body surface. "
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "Each basis is constructed with the function from `geomutil.onbfromaxe()` : orthonormal bases from axes. This function takes 2 sets of $n$ points $\\mathbf{p}_{A,n}$ and $\\mathbf{p}_{B,n}$ as input and provides an orthonormal basis as output. \n",
      "\n",
      "3 unitary vectors are constructed :\n",
      "\n",
      "\n",
      "$$\\hat{\\mathbf{w}}_n = \\frac{\\mathbf{p}_B-\\mathbf{p}_A}{| \\mathbf{p}_B-\\mathbf{p}_A |} $$\n",
      "\n",
      "$$\\hat{\\mathbf{u}}_n = \\frac{\\hat{\\mathbf{v}}_g - (\\hat{\\mathbf{v}}_g.{\\hat{\\mathbf{w}}_n}) \\mathbf{\\hat{w}}_n}{|\\hat{\\mathbf{v}_g} - (\\hat{\\mathbf{v}_g}.{\\hat{\\mathbf{w}}_n}) \\mathbf{\\hat{w}}_n|} $$\n",
      "\n",
      "$$\\hat{\\mathbf{v}}_n = \\mathbf{\\hat{w}}_n \\times \\mathbf{\\hat{u}}_n  $$\n",
      "\n",
      "Where $\\hat{\\mathbf{v}}_g$ is the unit velocity vector along actual trajectory.  \n",
      "\n",
      "The outpout of `geomutil.onbframe` is an MDA  $(3\\times n \\times 3)$ of $n$ unitary matrices aggregated along axis 1 $$\\mathbf{T}_n=[\\hat{\\mathbf{u}}_n, \\hat{\\mathbf{v}}_n, \\hat{\\mathbf{w}}_n]$$\n"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "To create the CCS : "
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.setccs()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 39
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "import scipy.linalg as la \n",
      "print \"ccs dimensions : \",np.shape(John.ccs)\n",
      "print John.ccs[0,:,:]\n",
      "print \"Check determinant : \", la.det(John.ccs[0,:,:])"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "ccs dimensions :  (11, 3, 3)\n",
        "[[ 0.94678656  0.05306765 -0.31745715]\n",
        " [-0.06834689  0.99696857 -0.03718026]\n",
        " [ 0.31452173  0.05689898  0.94754345]]\n",
        "Check determinant :  1.0\n"
       ]
      }
     ],
     "prompt_number": 40
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "Create a Wireframe body representation from the body graph model"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "Representation of frames associated with the cylinder"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.show3()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 41
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "On the figure below the wireframe model is shown associated with the 11 CCS (Cylinder coordinates systems) "
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "Image('CCS.png')"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "png": "Q0NTLnBuZw==\n",
       "prompt_number": 42,
       "text": [
        "<IPython.core.display.Image at 0x7fdd45c32a90>"
       ]
      }
     ],
     "prompt_number": 42
    },
    {
     "cell_type": "heading",
     "level": 2,
     "metadata": {},
     "source": [
      "Placing a dcs (Device Coordinate System ) on the cylinder"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "A DCS is refered by 4 numbers $(Id,l,h,\\alpha)$"
     ]
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "+ Id : Cylinder Id\n",
      "+ l : length along cylinder\n",
      "+ h : height above cylinder generatrix\n",
      "+ alpha : angle from front direction (degrees)"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "Id = 4 # 4 Left Arm\n",
      "l  = 0.1 # Longitudinal coordinates\n",
      "h  = 0.03 # height \n",
      "alpha = 45 # angle degrees"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 43
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.dcyl"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 44,
       "text": [
        "{'arml': 3,\n",
        " 'armr': 2,\n",
        " 'calfl': 9,\n",
        " 'calfr': 8,\n",
        " 'forearml': 5,\n",
        " 'forearmr': 4,\n",
        " 'head': 1,\n",
        " 'thighl': 7,\n",
        " 'thighr': 6,\n",
        " 'trunkb': 10,\n",
        " 'trunku': 0}"
       ]
      }
     ],
     "prompt_number": 44
    },
    {
     "cell_type": "markdown",
     "metadata": {},
     "source": [
      "Rotate Matrix around z"
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.settopos(traj=traj,t=6,cs=True)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 45
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.dcyl"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "pyout",
       "prompt_number": 46,
       "text": [
        "{'arml': 3,\n",
        " 'armr': 2,\n",
        " 'calfl': 9,\n",
        " 'calfr': 8,\n",
        " 'forearml': 5,\n",
        " 'forearmr': 4,\n",
        " 'head': 1,\n",
        " 'thighl': 7,\n",
        " 'thighr': 6,\n",
        " 'trunkb': 10,\n",
        " 'trunku': 0}"
       ]
      }
     ],
     "prompt_number": 46
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.show3(topos=True,dcs=True)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 47
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.show3(topos=True,pattern=True)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "ename": "ValueError",
       "evalue": "operand has more dimensions than subscripts given in einstein sum, but no '...' ellipsis provided to broadcast the extra dimensions.",
       "output_type": "pyerr",
       "traceback": [
        "\u001b[1;31m---------------------------------------------------------------------------\u001b[0m\n\u001b[1;31mValueError\u001b[0m                                Traceback (most recent call last)",
        "\u001b[1;32m<ipython-input-48-a9fdbef0bd42>\u001b[0m in \u001b[0;36m<module>\u001b[1;34m()\u001b[0m\n\u001b[1;32m----> 1\u001b[1;33m \u001b[0mJohn\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mshow3\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mtopos\u001b[0m\u001b[1;33m=\u001b[0m\u001b[0mTrue\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mpattern\u001b[0m\u001b[1;33m=\u001b[0m\u001b[0mTrue\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0m",
        "\u001b[1;32m/home/uguen/Documents/rch/devel/pylayers/pylayers/mobility/body/body.pyc\u001b[0m in \u001b[0;36mshow3\u001b[1;34m(self, **kwargs)\u001b[0m\n\u001b[0;32m    787\u001b[0m                 \u001b[0mkwargs\u001b[0m\u001b[1;33m[\u001b[0m\u001b[0mkey\u001b[0m\u001b[1;33m]\u001b[0m \u001b[1;33m=\u001b[0m \u001b[0mvalue\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0;32m    788\u001b[0m \u001b[1;33m\u001b[0m\u001b[0m\n\u001b[1;32m--> 789\u001b[1;33m         \u001b[0mbdy\u001b[0m \u001b[1;33m=\u001b[0m \u001b[0mself\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mgeomfile\u001b[0m\u001b[1;33m(\u001b[0m\u001b[1;33m**\u001b[0m\u001b[0mkwargs\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0m\u001b[0;32m    790\u001b[0m         \u001b[0mbdy\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mshow3\u001b[0m\u001b[1;33m(\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0;32m    791\u001b[0m \u001b[1;33m\u001b[0m\u001b[0m\n",
        "\u001b[1;32m/home/uguen/Documents/rch/devel/pylayers/pylayers/mobility/body/body.pyc\u001b[0m in \u001b[0;36mgeomfile\u001b[1;34m(self, **kwargs)\u001b[0m\n\u001b[0;32m    961\u001b[0m                 \u001b[1;31m#T = np.eye(3)\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0;32m    962\u001b[0m                 \u001b[0mT\u001b[0m  \u001b[1;33m=\u001b[0m \u001b[0mself\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0macs\u001b[0m\u001b[1;33m[\u001b[0m\u001b[0mkey\u001b[0m\u001b[1;33m]\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[1;32m--> 963\u001b[1;33m                 \u001b[0mgeo\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mpattern\u001b[0m\u001b[1;33m(\u001b[0m\u001b[0mAnt\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mtheta\u001b[0m\u001b[1;33m[\u001b[0m\u001b[1;33m:\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mnp\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mnewaxis\u001b[0m\u001b[1;33m]\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mAnt\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mphi\u001b[0m\u001b[1;33m[\u001b[0m\u001b[0mnp\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mnewaxis\u001b[0m\u001b[1;33m,\u001b[0m\u001b[1;33m:\u001b[0m\u001b[1;33m]\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mV\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mpo\u001b[0m\u001b[1;33m=\u001b[0m\u001b[0mU\u001b[0m\u001b[1;33m[\u001b[0m\u001b[1;33m:\u001b[0m\u001b[1;33m,\u001b[0m\u001b[1;36m0\u001b[0m\u001b[1;33m]\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mT\u001b[0m\u001b[1;33m=\u001b[0m\u001b[0mT\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0milog\u001b[0m\u001b[1;33m=\u001b[0m\u001b[0mFalse\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mminr\u001b[0m\u001b[1;33m=\u001b[0m\u001b[1;36m0.01\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mmaxr\u001b[0m\u001b[1;33m=\u001b[0m\u001b[1;36m0.2\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0m\u001b[0;32m    964\u001b[0m                 \u001b[0mbodylist\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0mappend\u001b[0m\u001b[1;33m(\u001b[0m\u001b[1;34m'{<'\u001b[0m\u001b[1;33m+\u001b[0m\u001b[0m_filepatt\u001b[0m\u001b[1;33m+\u001b[0m\u001b[1;34m'.off'\u001b[0m\u001b[1;33m+\u001b[0m\u001b[1;34m\"}\\n\"\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0;32m    965\u001b[0m \u001b[1;33m\u001b[0m\u001b[0m\n",
        "\u001b[1;32m/home/uguen/Documents/rch/devel/pylayers/pylayers/util/geomutil.pyc\u001b[0m in \u001b[0;36mpattern\u001b[1;34m(self, theta, phi, E, **kwargs)\u001b[0m\n\u001b[0;32m    668\u001b[0m         \u001b[1;31m# antenna cs -> glogal cs\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0;32m    669\u001b[0m         \u001b[1;31m# q : Nt x Np x 3\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[1;32m--> 670\u001b[1;33m         \u001b[0mq\u001b[0m \u001b[1;33m=\u001b[0m \u001b[0mnp\u001b[0m\u001b[1;33m.\u001b[0m\u001b[0meinsum\u001b[0m\u001b[1;33m(\u001b[0m\u001b[1;34m'ij,klj->kli'\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mT\u001b[0m\u001b[1;33m,\u001b[0m\u001b[0mp\u001b[0m\u001b[1;33m)\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0m\u001b[0;32m    671\u001b[0m         \u001b[1;31m#\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n\u001b[0;32m    672\u001b[0m         \u001b[1;31m# translation\u001b[0m\u001b[1;33m\u001b[0m\u001b[1;33m\u001b[0m\u001b[0m\n",
        "\u001b[1;31mValueError\u001b[0m: operand has more dimensions than subscripts given in einstein sum, but no '...' ellipsis provided to broadcast the extra dimensions."
       ]
      }
     ],
     "prompt_number": 48
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "Image('acs.png')"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "metadata": {},
       "output_type": "pyout",
       "png": "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\n",
       "prompt_number": 103,
       "text": [
        "<IPython.core.display.Image at 0x892be50>"
       ]
      }
     ],
     "prompt_number": 103
    }
   ],
   "metadata": {}
  }
 ]
}