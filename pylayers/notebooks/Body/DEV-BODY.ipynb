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
     "outputs": [],
     "prompt_number": 55
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
     "prompt_number": 56
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
     "prompt_number": 57
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
     "prompt_number": 58
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 59,
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
     "prompt_number": 59
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 60,
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
     "prompt_number": 60
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
     "prompt_number": 61
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
     "prompt_number": 62
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 63,
       "text": [
        "(3, 16, 300)"
       ]
      }
     ],
     "prompt_number": 63
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
     "prompt_number": 64
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
        "      <td>-0.463661</td>\n",
        "      <td> 0.061186</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td> 1.700309</td>\n",
        "      <td>-0.000229</td>\n",
        "      <td> 0</td>\n",
        "      <td>-4.061282</td>\n",
        "      <td> 0.000000</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.204082</th>\n",
        "      <td> 0.061186</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td> 1.236648</td>\n",
        "      <td> 0.060957</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td>-2.360973</td>\n",
        "      <td>-0.000458</td>\n",
        "      <td> 0</td>\n",
        "      <td> 3.168114</td>\n",
        "      <td> 1.713606</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.408163</th>\n",
        "      <td> 0.122143</td>\n",
        "      <td> 0.408163</td>\n",
        "      <td>-1.124325</td>\n",
        "      <td> 0.060499</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td> 0.807142</td>\n",
        "      <td>-0.000684</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.989101</td>\n",
        "      <td> 4.084166</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.612245</th>\n",
        "      <td> 0.182642</td>\n",
        "      <td> 0.612245</td>\n",
        "      <td>-0.317183</td>\n",
        "      <td> 0.059815</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td> 1.796243</td>\n",
        "      <td>-0.000909</td>\n",
        "      <td> 0</td>\n",
        "      <td>-2.425976</td>\n",
        "      <td> 4.918904</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>1970-01-01 00:00:00.816327</th>\n",
        "      <td> 0.242457</td>\n",
        "      <td> 0.816327</td>\n",
        "      <td> 1.479059</td>\n",
        "      <td> 0.058906</td>\n",
        "      <td> 0.204082</td>\n",
        "      <td>-0.629733</td>\n",
        "      <td>-0.001129</td>\n",
        "      <td> 0</td>\n",
        "      <td> 0.214830</td>\n",
        "      <td> 6.727692</td>\n",
        "    </tr>\n",
        "  </tbody>\n",
        "</table>\n",
        "<p>5 rows \u00d7 10 columns</p>\n",
        "</div>"
       ],
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 65,
       "text": [
        "                                   x         y         z        vx        vy  \\\n",
        "1970-01-01 00:00:00         0.000000  0.000000 -0.463661  0.061186  0.204082   \n",
        "1970-01-01 00:00:00.204082  0.061186  0.204082  1.236648  0.060957  0.204082   \n",
        "1970-01-01 00:00:00.408163  0.122143  0.408163 -1.124325  0.060499  0.204082   \n",
        "1970-01-01 00:00:00.612245  0.182642  0.612245 -0.317183  0.059815  0.204082   \n",
        "1970-01-01 00:00:00.816327  0.242457  0.816327  1.479059  0.058906  0.204082   \n",
        "\n",
        "                                  vz        ax  ay        az         s  \n",
        "1970-01-01 00:00:00         1.700309 -0.000229   0 -4.061282  0.000000  \n",
        "1970-01-01 00:00:00.204082 -2.360973 -0.000458   0  3.168114  1.713606  \n",
        "1970-01-01 00:00:00.408163  0.807142 -0.000684   0  0.989101  4.084166  \n",
        "1970-01-01 00:00:00.612245  1.796243 -0.000909   0 -2.425976  4.918904  \n",
        "1970-01-01 00:00:00.816327 -0.629733 -0.001129   0  0.214830  6.727692  \n",
        "\n",
        "[5 rows x 10 columns]"
       ]
      }
     ],
     "prompt_number": 65
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
       "metadata": {},
       "output_type": "display_data",
       "png": "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\nQmCC\n",
       "text": [
        "<matplotlib.figure.Figure at 0x947e590>"
       ]
      }
     ],
     "prompt_number": 66
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 67,
       "text": [
        "'t (s) : 0.0:9.591837\\nd (m) : 51.826645945\\nVmoy (m/s) : 5.40320336397\\n'"
       ]
      }
     ],
     "prompt_number": 67
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
     "prompt_number": 68
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
       "metadata": {},
       "output_type": "display_data",
       "png": "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\nAAAAAElFTkSuQmCC\n",
       "text": [
        "<matplotlib.figure.Figure at 0x7fbd4021b490>"
       ]
      }
     ],
     "prompt_number": 69
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 70,
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
     "prompt_number": 70
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 71,
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
     "prompt_number": 71
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 72,
       "text": [
        "True"
       ]
      }
     ],
     "prompt_number": 72
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
     "prompt_number": 73
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
     "prompt_number": 74
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
     "prompt_number": 75
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
     "prompt_number": 76
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 77,
       "text": [
        "(3, 16, 300)"
       ]
      }
     ],
     "prompt_number": 77
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 78,
       "text": [
        "16"
       ]
      }
     ],
     "prompt_number": 78
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
     "prompt_number": 79
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 80,
       "text": [
        "(3, 300)"
       ]
      }
     ],
     "prompt_number": 80
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
     "prompt_number": 81
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 82,
       "text": [
        "(2, 300)"
       ]
      }
     ],
     "prompt_number": 82
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
       "metadata": {},
       "output_type": "display_data",
       "png": "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\nmsSQQ0RERJrEkENERESaxJBDREREmsSQQ0RERJrEkENERESaxJBDREREmsSQQ0RERJrEkENERESa\nxJBDREREmsSQQ0RERJrEkENERESaxJBDREREmsSQQ0RERJrEkENERESaxJBDREREmsSQQ0RERJrE\nkENERESaxJBDREREmsSQQ0RERJrEkENERESaxJBDREREmsSQQ0RERJrEkENERESaxJBDREREmsSQ\nQ0RERJrEkENERESaxJBDREREmsSQQ0RERJrEkENERESaxJBDREREmsSQQ0RERJrEkENERESaxJBD\nREREmsSQQ0RERJrEkENERESaxJBDREREmsSQQ0RERJrEkENERESaxJBDREREmsSQQ0RERJrEkENE\nRESaxJBDREREmsSQQ0RERJrEkENERESaxJBDREREmsSQQ0RERJrEkENERESaxJBDREREmsSQQ0RE\nRJrEkENERESaxJBDREREmsSQQ5oRFxcHvV6P/Px8axeFiCwkIiICI0aMsMi0g4ODsWTJEotMm6yD\nIYeqhJeXF7Zt2/bA5/n7778/0HkCwJtvvglPT084OjrCw8MDr7/+OnJzc8scb8yYMdDr9bhw4YLa\nbfXq1fj73/+OunXrokePHibD37x5E507d0bDhg1Rr149GI1GrF27Vu2flZWFf/7zn3B3d0f9+vXx\n8ssvm5QjODgYdnZ2MBgMMBgM8PX1NZn+2rVr0bZtWzg6OqJt27ZYt26dSf/Dhw+jW7duMBgMaNy4\nMSIjI036L1iwAM2bN4eDgwPatGmDc+fOAQB+/fVXdOnSBc7OznBzc8M//vEPpKenm7XMADBu3Di0\nbt0atWrVwtKlS0usz9mzZ6N79+5Fut+4cQO1a9fGqVOnShy3orKzs+Hr64umTZuWOlxUVBRatGgB\ng8GAPn36IDk5We23fft29OjRA05OTmjWrFmx45dUt1euXEFoaCjc3d2h1+sRHx9vMl5WVhbGjBmD\nevXqwc3NDZ9++qlJ/6NHj6JDhw6oW7cuAgMDERsbW+z8e/bsWWSHwcHBQV2XDAYDbGxs8OqrrwIA\nTp06hcDAQNSvXx9OTk7o3Lkzdu/ebXb99ejRAy4uLnB0dISvry8WL15cUtVCp9OV2K8isrOz0ahR\nI9y9exc6na7Kp0/WxZBDVcIaGwedTgcReaDzBICxY8fi1KlTSE1NxYEDB7B582ZERUWVOs7u3btx\n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\n",
       "text": [
        "<matplotlib.figure.Figure at 0x8621850>"
       ]
      }
     ],
     "prompt_number": 83
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 84,
       "text": [
        "<matplotlib.text.Text at 0x8918250>"
       ]
      },
      {
       "metadata": {},
       "output_type": "display_data",
       "png": "iVBORw0KGgoAAAANSUhEUgAAAcMAAAEZCAYAAADrI06XAAAABHNCSVQICAgIfAhkiAAAAAlwSFlz\nAAALEgAACxIB0t1+/AAAIABJREFUeJzt3XlcVFX/B/DPICiyiCIIyKrgAoqsipoC7mlKppZaKoop\nam71aLnUI6bZZplLmZlhmltpT6GhD2qOmgomgkuUuCGb4Io6CAIz398f9+H+GHYU5s7yfb9evmRm\n7sz9nnvvnO+cc+89R0ZEBMYYY8yAGUkdAGOMMSY1ToaMMcYMHidDxhhjBo+TIWOMMYPHyZAxxpjB\n42TIGGPM4EmSDNPS0mBkZASVSvVU79+2bRsGDRpUz1HV7MSJE2jXrh0sLS0RExOj8fWXmj59OpYv\nXw4AkMvlcHZ2Fl/r3Lkzjh07JlVotTZx4kS899579fZ56enpsLS0hKbuFMrNzUVwcDCaNWuG+fPn\na2SdmqLpbcnqz5AhQ7B169ZKX3vWerc2ytZNOockcP36dZLJZKRUKut12YbWt29fWrNmjdRhqDly\n5Ag5OTlJHUadTZw4kd577z2pw3hq77//Po0cOVLqMOqFq6srHT58WOownlp0dDT16tVL6jCe2ZIl\nS2jcuHEN9vk11aXachxItT91ppuUtOBXanp6Ory8vBp0HUSkFWWti5KSkqd6n66Vs6wbN27A09Oz\nyteVSqUGo3k2MplMp/fFs3ra41fTGrpuqOk40JXt9NTfvZqyZVZWFo0YMYJsbW2pTZs2YssoKyuL\nmjZtSvfu3ROXPXv2LNnY2FBJSQkplUpatmwZubq6UqtWrWjChAn04MEDIqr4C8XV1ZUOHTokfk7Z\nX0jOzs4kk8nIwsKCLC0t6dSpUxV+OZw4cYICAwPJysqKunbtSidPnhRfCwkJoffee4+ee+45srS0\npIEDB9KdO3eqLO8333xDHh4eZG1tTWFhYZSdnU1ERG3btiUjIyNq2rQpWVpaUlFRUYX3pqen00sv\nvUS2trbUsmVLmjlzZoXyVFb+kJAQWrx4MfXs2ZOaNm1KH3/8MQUGBqp99ueff05hYWFERBQeHk7v\nvvsuEVVsGZb9dbdkyRJ6+eWXacKECWRpaUmdOnWiM2fO1LhviYgSEhKoe/fu1Lx5c3JwcKCZM2eq\nlVkmk9GXX35JHh4e1LZt20q35ahRo8je3p6srKwoODiY/vrrL/G1iRMn0rRp02jAgAFkaWlJISEh\ndOPGDfH1uXPnUqtWrahZs2bk7e1NFy9eJCKix48f01tvvUWurq5kZWVFvXr1osLCwgrbNDo6mtq2\nbUuWlpbUpk0b2rZtGxERXb58mYKDg8nKyopsbGxo9OjR4jr//vtv6t+/P1lbW1OHDh3oxx9/rLRc\n4eHhZGJiQo0bNyZLS0s6dOgQLVmyhEaOHEnjxo2jZs2a0aZNmygrK4uGDRtG1tbW5OHhQRs3bhQ/\nY8mSJTRq1CgaN24cWVpakre3N6WmptKKFSuoVatW5OLiQnFxcZWuv3Q/f/rpp+Tt7U0WFhYUERFB\nOTk59Pzzz1OzZs2of//+dP/+fXH5X3/9lby8vKh58+YUGhpKf//9NxERjRs3TjyuLSws6NNPP62w\nLWsqR3XHWHkXL14Ut7GdnR2tWLGCiGp3vK1Zs4batm1LNjY2NH/+fFKpVJSSkkKmpqbUqFEjsrCw\noBYtWhCR8J369ttvxfeXrzMqO3737t1LPj4+1Lx5c+rZsyedP3++zuVQKpX04Ycfkru7O7Vs2ZJe\neeUVsY4s3a7ff/89ubi4kI2NDX3wwQdERLR//35q3LgxmZiYkIWFBfn6+orlKFs3XL16tcb6rrTc\nJSUl9K9//YtsbGyobdu2tG7duipbhtUdB5s2bSIXFxcKCQkhouq/12Xrppq2aWX15d9//01NmjSp\nsD/z8vJo/PjxZGtrS66urrR8+XJSqVTivu3Zsye9+eab1LJlS1q0aBFZW1vThQsXxHXl5uaSmZlZ\ntXV/tclQqVSSv78/LVu2jIqLi+natWvUtm1b+u9//0tEQrdh2S/GvHnzaPr06UREtGnTJvLw8KDr\n16+TQqGgESNG0Pjx44moYjJwc3NTa55HRUWJySMtLa3CDix7YN+9e5eaN29OP/zwAymVStqxYwe1\naNFCPABDQkLIw8ODLl++TAUFBRQaGkoLFiyotLyHDx8mGxsbSkpKoidPntCsWbMoODhYfL18nGWV\nlJRQly5d6K233qLHjx9TYWEhnThxokJ5Kit/SEgIubq6UkpKCimVSnrw4AFZWlrS5cuXxfcEBgbS\nrl27iEi9i7F8Miwb45IlS8jU1JT2799PKpWKFi5cSN27d6/Vvk1MTKSEhARSKpWUlpZGnp6e9MUX\nX4jrkclkNHDgQLp//z4VFhZWuk2io6NJoVBQUVERzZ07V/yCEwlfGktLSzp+/Dg9efKE5syZI+7T\n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\n",
       "text": [
        "<matplotlib.figure.Figure at 0x8929a90>"
       ]
      }
     ],
     "prompt_number": 84
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
     "prompt_number": 85
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 86,
       "text": [
        "9.97"
       ]
      }
     ],
     "prompt_number": 86
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 87,
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
       "metadata": {},
       "output_type": "display_data",
       "png": "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\n",
       "text": [
        "<matplotlib.figure.Figure at 0x7fbd40919850>"
       ]
      }
     ],
     "prompt_number": 87
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
     "prompt_number": 88
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
     "prompt_number": 89
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 90,
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
     "prompt_number": 90
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
     "prompt_number": 91
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 92,
       "text": [
        "(<matplotlib.figure.Figure at 0x7fbd40928550>,\n",
        " <matplotlib.axes.Axes3DSubplot at 0x7fbd40919110>)"
       ]
      },
      {
       "metadata": {},
       "output_type": "display_data",
       "png": "iVBORw0KGgoAAAANSUhEUgAAAV0AAADtCAYAAAAcNaZ2AAAABHNCSVQICAgIfAhkiAAAAAlwSFlz\nAAALEgAACxIB0t1+/AAAIABJREFUeJzsvXecFPX9P/6ctvX2CkivChgQkCKIChZiEEFCQNQoUfkA\nxoYxGhPNL1asMX5tCfZEUDGCH/SjICWxEkWQqEQRLGBBjhrK3d622Snv3x/De252bnZ3Znd27xbn\n+XjcQ7mbec972nNe7+erMYQQAg8ePHjwUBawrT0BDx48ePghwSNdDx48eCgjPNL14MGDhzLCI10P\nHjx4KCM80vXgwYOHMsIjXQ8ePHgoIzzS9eDBg4cywiNdDx48eCgjPNL14MGDhzLCI10PHjx4KCM8\n0vXgwYOHMsIjXQ8ePHgoIzzS9eDBg4cywiNdDx48eCgjPNL14MGDhzLCI10PHjx4KCM80vXgwYOH\nMsIjXQ8ePHgoIzzS9eDBg4cywiNdD0VDURTIsgyv3Z4HD/nBt/YEPFQmCCEghECSJKTTaciyDIZh\nAAAcx0EQBHAcB5ZlwbKs/jcPHn7o8EjXgyMYyTYej4NlWfA8D4ZhwLIsRFGELMtQFCVjP5ZlwXGc\n/uORsYcfKhivBbsHOzCSraqqAIBEIgFVVaEoCgghOoEyDANBEHRiNY9hhEfGHn5o8EjXQ04QQqCq\nKmRZhqqqYBgGqqpCFEWkUilwHIdgMKhbtul0WidgVVX1/6dkSonVSKrG7QBAlmXwPK9LFB4Zezic\n4MkLHixhRbaEECSTSaTTafh8Pvj9fp0MVVXNIEafz5cxDv2hljIhRN/euB/DMFAUBSzL6nqxkWg9\ny9hDpcMjXQ8ZIIRkRCMYLdt0Og2/34+amhqwLItkMplhoRrHoGAYRidI8zbZyBjQrF2jVUyPYXTe\nmcmY53l9H47jMvbz4KGtwCNdDwCayTaRSIDnebAsC0IIEokEJEnKINtcoBZxPuQi40QioR/fjmVM\n9zNKG5Ikged5nYg9MvbQVuCR7g8cZss2kUigqqoKyWQSkiQhEAggFArlJVu3QAmRkqVxnk5kClmW\nIQhCC0cf/a9RorCyqD14KBU80v2BIpuMQAhBLBZDIBBAOBzOS0Ll8sM6lSkAQBTFFmRMPx5WURcA\nPDL2UHJ4pPsDAyFEj6M1km0ymdQTHCKRCHg+/6PRFogoGxnHYjEIggAAth14dFtjJAaVWugx6L89\nMvZQKDzS/YGAkq0sywCgRwmkUikoioJAIICqqio0NjYeFkRCydgsi1ALn1rGNDoDaI6MsIrIoGSc\nTqczxvPI2INTeKR7mMNMtoBWKyGZTEJVVZ1sjYkNxUoGbkgOpZItqF5sPpZRplAUBZIk6SFyRiI2\nEivdx0zGVMbwyNiDFTzSPUxByTaRSIAQAr/fD1mWkUqloKoqgsEgfD5fUSRgRdBukEq5iYmSodkq\njsfj8Pv9AJBBxjRu2UzEVmQcj8d1ovfI2APgke5hB7Nla7bEAoFA0WRrxuFKGJQMrZx3VpaxFRkb\noyXyWcbGOGOPjA9feKR7mIC+zEYZQZZliKIIAAiHwxAEIe9L7Ia8cLjAKvEDyG4ZW5Ex/ZFlOa9l\nTO+V8ThW2XflCt/zUBp4pFvhoM4gY1UvSZKQSqUAQPfg07RcN+ERdCasyDiVSulWrB3L2Eyq1CpO\np9MIBAL6ccwSBU368ND24ZFuhSIb2SaTSbAsi2AwCEEQkE6nIUmS7XE9InUfRmI1wq5MQe+JmYyN\nqxr6dyMZWxUX8tD68Ei3wmD0rFOk02ndogqHw7rTpjXgkbZ95JMpaGgbJVfqlMtVsc2c8GImY69I\nUOvDI90KgLGWrSiKEEURVVVVOtlyHKdrtmZ4JFh5MJOxoigQRRHBYNB2KrQVGScSCb3WMQB9O4+M\nywuPdNswrAqH0wyyxsZG8DyPqqoqW9ljduGEpD1Cz49szrhCUEjFNnMkhTlhxCxTAF75zFLDI902\nCHMtWwpaOJwQgurqattkWw5i9F7I0iEfcdslY1qJLZ1O57WMZVlu4QswShrmOGMP9uGRbhtCtsLh\nVFIQBAGhUAipVMo24Tp9ITzr9fCBmYxVVdW7cRQqU8iyjGQyCZ/Pp1vMnmXsDB7ptgFYVfwihCCV\nSulkW11dDY7j9EI1TscvBczZV2ZLy0PbhBsyhdkBmMsy9sg4Ex7ptiKoXitJUkZoECVbn8+nky2F\nU0u0lA84fUmbmpp0rZnCqqxipcFNPbYS5mGXjAHoceB2LONsLZd+qDKFR7qtAKNlm06n9WiERCKh\n9x+z06WhFLBD6sbGlAAQiUT0eGHjS0YD+3MtWw/3F62tEDdQ+FzMZCxJEsLhMIDmTEirim3Zsu+o\nsZFMJhEIBPQ5/VAsY490ywgrGQHQPMiNjY22WuIUYum6JS9QKzyVSmXoy1T2oEtOhmEyMuCcdH0w\nvpwe2h7y9b8zJ3yYfRTmrDujZJGv/53RQq5kMvZItwzIptnS/mMAHFu2pbSgzIRnJlun+rLdZWs6\nndYtJWo9VbpE4Qbcutduh69l+73duhQA9NjhfJaxsf8dhTnGuFL633mkW0JYdWkwkm0gEEAgEEAs\nFrNNuKWORjBun41srVDIC53NUqLZdVSiyFanoFJeMg8tyZi+G+Fw2HbFtmxknEql9GiKTz75BFu3\nbsXs2bNb+YyzwyPdEsCqS4OxJY6x/5ixdoLTY5TS0k0mkzrZ2mnf46YFRV82Y4ad0VIy6uHZJApP\nnigN3La6nVjGucpn0o/wjh07sH///qLnV0p4pOsiaEytLMu6FUa7NFCyNXZpAArTXJ089E7Gpx8L\nSZJskW22sUtBeMaX0zinbHox0FzhK1udgnLgcHCktQacyhTJZBLnnXceWJZFKBRC165dMXDgQAwc\nODDDv7B9+3ZcfPHF2Lt3LxiGwaWXXoqrr766xfGvvvpqrFy5EqFQCAsWLMCwYcNcOzevMKcLoA6A\nVCqlSweKoiAWi+kNEmtraxEMBls89IWSrpuxunRJ39DQoAfQ221OaTW3coJKFIIgwO/3IxgM6p51\nQRDAsqwebRGPxxGPx5FMJiGKon6f2rpV3JbIsrX1ZeOH1+fz6Z09QqEQHnroIQwdOhTV1dVYtmwZ\nLrroImzatCljf0EQ8OCDD2LTpk1Yt24dHnnkEXz++ecZ26xYsQJbt27Fli1b8OSTT+KKK64o/EQt\n4Fm6RcCqcLjRA2tl2WZDqV6sXGNSyzyZTILneUQikRblIisVlIyNlpLdJevhWJu2rX1Y3H7eWZbF\n0UcfjVAohFmzZmH8+PGW23Xu3BmdO3cGAFRVVWHAgAHYuXMnBgwYoG+zdOlSzJgxAwAwatQoNDQ0\nYM+ePejUqZMrc/VItwBk69JAO+tyHIdIJGLroSr0a1/MS0TJloZ7Ga3aSrD8CkW+JSuNmDBKFLQ+\n8eESX1zsvNuS1Q20nE9TUxNqa2tt7fvdd99hw4YNGDVqVMbvd+zYgR49euj/7t69O+rr6z3SbQ1Y\nFQ6nueiA1n+M6qKF6K6lepiNJG0k22IrlLXleFqnERtmMiaE6I0ps8UXm1volNKx2RqJMqVEqWSK\naDSKmpqavPvFYjGcc845ePjhh1FVVWU5rhFu3luPdG0gH9nSLg0Mw+gedScodcKDsWhOPrJty0Tq\nFG68KFbxxUbHXb4Qp7ZiGba2FlsuNDY2oq6uLuc2kiRh2rRpuPDCCzFlypQWf+/WrRu2b9+u/7u+\nvh7dunVzbY4e6WYBXXImk8mM5aSx/5iRbM37lhJ2iZHGMQLNqZtWhc7dxuFC2tlgzKiiyKYXA1oC\ngJVV3JbJq9QopaWbS14ghGD27Nk45phjcM0111huM3nyZMybNw/nn38+1q1bh9raWtekBcAj3RYw\npiOqqqovV6z6j1k9NK2h0ZpByZbOl2EYhMNh16uA0XM1p4b+EJFNL47FYggEAo5SoK0+4m1Fi23r\n49Bwx2xYs2YNFi5ciGOPPVYPA7v77rvx/fffAwAuu+wyTJw4EStWrEDfvn0RDocxf/78oudphEe6\nh0CtFHPhcEAT51nWXv+xcoSA5YqPNZIttWwbGxttj384yQttBVT/NcIqvpjKV2ar2Lsf9kCvU673\nc8yYMS3ebyvMmzfPtXmZ8YMnXTPZ0htG0wsBTUag8YD50BqklY1sKxmHA9HkOodsKdDZQtqo8zNX\nGcVywE0L1Q3noNV82vpq6wdLutnIltYaoP3H4vG4o2V5OS1dK7K1ssQLcby1Jtr6S+MUds8nm0QR\nj8czEj0KKZnZ1h1ghcJ4Xq393NrFD450aTymubyikWyNcavlsFwLIUXq0LMreziZi9PtD9cXuq3A\nyioGnJXMdOsZLqWFWixSqRQCgYCrY5YCPxjSpfGzxiI0APSU0GxVtIqxQt1+qCjZUq84lRHyHae1\nJI9KJuJKmL8dMqYGRr5kj3KjFDJFY2OjrRjd1sZhT7pGyzYWi8Hn80EQhLxkax7DLkoRvUDJlsYF\nU6I1FvJwC/nnAmzbBqRSDLp3b7ldWyeqHwKMZEy1/XQ6DUVRIAhC3vjiXCUz29rHyDgfu4kRrY3D\nlnSzyQiiKCKRSFj2H7NCMSTqRMvLFo1gJFsaqkbbsDudT7FQVeCpp1i8/TYLjgOqqoCrrmLRt286\noyeatq36gy48DrSt8CpKouakGKPzrpwlM0txbRoaGmynALcmDjvSpTKCuXB4MpmEJEnged5Rl4Zy\nOMbMyEa29OFiGMZW2Ivb+PRTBm+9xaJrVwJFARoaVCxYEMQttyR1i8q4jAUymxD+0DtAtCaykZzR\neWenZCZ13sqyXFTJzFLIXZ6lW2ZYabaUbNPpNPx+P/x+v2U2US6U05GWj2yLHb/YbQ8c0OSFPXsI\nli3j0bEjgw4d/GhoYNG5s6LvT7PfsoU/0Zfc3O/KI+O2g2x6cTwe18k5V/8zOx9Yty1dT9MtE6zI\nVlVVva6tsdljMpksi9VayD6KoiAajQLIT7atlcDQqZMEWSZo317B//wPwaefcojHGZxwQghHHqni\nrLMk/PSnCmiaerZ02XwedypNtDX9sLXQljLSAOhNIo1j5/vAGj+ybsJMup68UELkI9tAIIBQKJRx\nkwtZlhe6j93aCLIsQxRFEEIQCoXg8/lKEo1QjKUryzISiQS6dVMxa1YVXnghAFUFxowhuPTSBnTs\nGMK77zJYvpzD2WcHIQh+TJqk4Kc/VTBypAKjscQwDGSZQyzGIRQCgsHm+dHSioqiQFEUxONxy5fW\ns4rbFrLFF+f6wALI8AMUcl/Nz2k0GkWvXr2KP6ESo+JIlxJVLBYDx3Hw+XwtWuLQ/mNmlMtqtRON\nQKuU0U4NqqraznpzikIJis5RURQEAgH4/X5MmsRg3DgZySRQUwNEowSCAJx2moJTT5Xxpz+lsW6d\niDffrMJvfuPHnj0MJkyQcdZZMk47TcGBAwyWLeORTgMcB0yYIOOoo0iGk4daL4IgZDh5ssWhtlYr\nnnw4HC11pw5iK4mCvq8cx7lSMtMYvZCvwlhbQMWRrjHu0Cjq2+nSUE6nWLZ9qGarqiqCwSB8Ph8k\nSYIoiiWdk9PtY7EYJElCMBhscV39fu3Hem7AscfKGDVKxE03SfjuOwbLl/P4y198uOQSDr16qRg6\nVMGQIQrCYWDFCh6zZkkIhazGyu/kMT4PRqsYKC6Coq0QZltKanBT0qIfVSOclMy0kp88R1qJwLJs\nRiiYFSlkQzktXTOsyNa4XSk1WrvnoCiKXm+CZVnU1tYWLXX07k0wZ46EOXM0Ar7tNj82buTw/PMC\nbr5ZhM8HxOMMQiFnIXDmOFSzrghAD62zemErrSh4W/gAUJQq89FJyUxKxoQQHDhwAO+++y5isRiq\nq6uLnlupUVlPHrQgb1oPwefzIRAIFB0Pmw/FELUkSYhGo3oXgpqaGvj9/ow5O32IjeN//DGDK67g\ncNllHNaty90PLRtUVUU8Hkc0GtXnYtVEs1h0705w0kkyZs9OY8oUGUuXCuA4IBwuvqOw0SKmSSOh\nUAjhcBh+vx8cx+kro0QiUZENKotBW7HcAecShfG+Ul9NOBzWDZeDBw/imWeewT//+U8MGTIExx9/\nPH71q1+1GGvWrFno1KkTBg8ebHmsd955BzU1NRg2bBiGDRuGO++8s6jzzIaKs3QFQUBNTQ1EUXTc\nQNEtq9XOPrIsIxqNZrVsi50XAHz0EYOf/UxAMqkt6197jcMLL0g45RR7rUZo9SpRFOHz+fQoD2rt\nug2eB376UwVLl3I49lgFb7/No0sX1VJacAtWumIu68lobdFtPMedhrZE3EbpqU+fPnj55ZcxceJE\nvPLKK/jyyy+xa9euFvvMnDkTv/rVr3DxxRdnHffUU0/F0qVLSzn1yiNdo5ezLTjFzKAZPaqqIhQK\ntbBqs6GQDLNHHuGQSmlkpqqAKAIPPcThlFPknPsTQvQCP1aZeXT8Ukg2XboQzJwpIx4HOnUimDfP\nj6lTEzAeqtQhcXa87VS+SiQSAFrWuLWjFbe1UK+2glJl6hFC0K5dO5x00kmW25988sn47rvv8o5Z\nalScvEBRDIGWIoVWlmU0NTVlRFXYlT4KlRcOdeIBIUAqBcRiwLp1LP78Zw5btjAttlcUgvXr01ix\nIo7du1VUV1e71lHCLJfkumY+H1BXB0yfLkOWgf/7v7bx7adWsSAI8Pl8egW3UCikx00rigJRFBGP\nxxGPx5FKpXRnrqqq+P57YMkSHv/3fzz27m07ZNmWPgClIDY3xmQYBu+//z6GDBmCiRMnYvPmzS7M\nrCXaxtPuAPSml1MqoERttb85rKqqqsqx9FGoZTdrloK33mIhy1q8a1WV9rsvvmDw0EMCQiGCCRNU\nnH46i2HDFNxyi4x33glCEEIQBAbz5skYOtS6aE1ZvvgscMcdIq64IoBu3VRUVwNHHaWiLRl2TmJQ\nv/6awe23VyGdZgAwePnlAObOlXDkkYeHtep29EKxsHonixl3+PDh2L59O0KhEFauXIkpU6bgq6++\nKnaaLVCRli59CQqpP+CUULLdRGrZNjU16ToztWwLOUYh2//4xwTz58s47jiCoUMJHn1Uxp13Knj0\nURlff53Gc89JqK6WcfvtPgwYcASeey4MSWLBcQxUFbjpJnd7phWCXr1URCIE8+b58K9/sVixgtct\n+LYMo1Xs9/sRDAaxalUEHMcjFGJRUwMkEsA//ynoVjF13BlDHvOhVEvxQtGW5mIch+rxxSASiSB0\nyMEwYcIESJKEAwcOFD1PMyqSdIHCrbFidV0z2dbW1raQEcoVJUEIwfjxKlaulPDPf0qYOlXVfy9J\naRx1VCOuvjqKN96I46674qiuJvjvfxn85z8MQiFgzx53XsJiLKAPPuAwfbqEt97iUVsL7N0L7NpV\nmY9lIsHA5wOiURYrVviwYwePVEqTKILBoB5rTMMH4/E4EomELlEc7hEUpUQ0GkUkEilqjD179ujX\nf/369bpG7DYqTl6gKDfp0nRdWZZtxQa7YU0Xsr0kSbrzh9ZwkGUZxx2XRl0d0KMHwSefMNi5Exg9\n2nqO5ZIXtPkCPXsSTJok4+BBBiyrVTBrbRRijZ12moING1h07Kji1FNVrFvH45VXgjjnHBHdu7dM\nTjBm22UrHkO3awsONTctVLe6T9BxotFo3hjdCy64AKtXr8a+ffvQo0cPzJ07V28IcNlll2HJkiV4\n7LHHwPM8QqEQFi1aVPQcrVCRpGv2WDp5EAqJRiCEIB6P28p6M8/P6byc7GvcPl/yRf/+Cm66ScE9\n93CoqiJgGOCuu3JHORQDu9d44EAV69ZxOP10GaKoRWF07Fj+spVu4OSTFcgysHw5jw4dgF//uglr\n1vA45ZQQ7rhDxPTpcosojXzhbMbnzxw9QYuN24FbjjQ34OY49JzsFLt54YUXcv59zpw5mDNnjitz\ny4WKJF2g2cFRKtKl+eHG1jh2OzWUw1Kk52x05GWLB6bzmTJFxVlnqdi0CZg82Ydcaep25291rk7u\nx5AhKnge+PZbFrW1BCNGKKiq0iIyKg0MA/z4xwp+/GPNVJckBSNHpnHWWcBllwWwbBmPhx8W0bEj\nwd69DNJpoGNHkpFSbXbcsSyrF3AyOu5oWKLRKi5HzWK3xnVb062Uso5AhZJusREMufYxki0tntPU\n1OTqMdzaJx6P25Y7KAQB6NVLi+097TQew4cD110no3fvzLmUCywLDB6sYvDgZuu2EhxpTjB4sIp3\n3kngj3/0YfToEM49VwLLaudeW0swc6aEfNKh3SQPq64P9H4Wa+2WwkJ1Cx7plglukpsV2RrJvbXm\nlW2e9MVyWh9BVYE//IEDzxPs2MHA5yP49a95LFwooxg/RGvrjW0RRnLx+YBbbkljwAAV113nR69e\nBOPHS4jHGaxYwePCC62lnlwEZTecjZY/TSQSBSV5mI/ZVmC8NtFotCJq6QIVHL0AuENuiqIgFosh\nGo2C4zjU1NS0qDtQyHGKtQriceDBBzn8+tccXnqJhaI010egL4vdbDcjGhqATZtY9O9PsGcPg6oq\nIBplLJMp7EJVVSSTSd0L7zQBpdRYu5bDs88K+Ne/uFaXLXr0UHHBBRJkGXjoIT/SaWD3bndfQ6tw\nNgCOkjys7l9bs3QrlXQr0tJ1Q14wWrbG7hK59nE6P6eFPegxRBEYN07AV19put/f/07w73+LuPlm\n6POkWrPTsQMBTXvkee3/UylAUZqLidvF118Dq1cDySSPESNSGDSI14+jqipEUUQ6nbbdZbZUuO8+\nAY8/7oOqasv5Cy+UMHdu6+kXHToQcByDX/xCwscfc/jrX3249lr7ZT0LAX0O7XTy+OQTFWvXsggG\nVZx+uoQuXZr3a0sZaeYxGhsbMWjQoKLHLQcq2tItJEFCS59N6xZjTU1Niw4TZriVUJFvH3qMt99m\n8c03zKHQKQJZZvDEE2EIQlifZ6HOulAImDlTwd69DAQB2LWLwZgxKn70I/tjffstwWOPqfjmmzT2\n7GGxeHEE27YFdMuK4zjdwqKWFXX4GWNTc1X3cuPF3LOHwWOP+cCy2oeGYYDnnhPw3Xett0Tu04dg\n4kQJe/cy6NpVxcUXS3j8cR+WLWsd+8doFW/cGMQDD1Tj449DWL06hLvvrsa+faxuFVPNmPYddJLk\nYXVct+YPVE4tXeAHZOnSWrG0RUi5OgI7fbgIIWhsTAPgIMsMVJXR291IkqYNOp2TedsZM1T0709w\n0008hg5VcfvtCoyXItvYhGiFctasIfD5/IhGA0ilVHTurOLDD1m9AwQ951yWFW3JY6yNSvVG+iIX\na1kdPKhdO4YBGhsZyLJWfP1//5fHJZdIltEb8bi2X10dYPPxcIwxY1SMHJmGJAHhMDB9uoSf/zyI\nbdsYjBypoKGBQefOBIMGqWWtmbB0KY/qaoKGBu0eNjUx2LAhgEmTtGiMRCIBntdWNNn625mbjVrN\nxQ2Yz8mTF8oEO8RD9UbaETgcDuvLXifHKVWfNCNordeRI3kQEoQkaRIAywInnkgQDjsaLsfcgBNO\nIDj1VBWqqkkNuUBXB7TFSnV1FTiOw/r1HPr2VaEo+cdoPra1F97cNYDGphbj+OndW0VVFcGBAwza\ntydIJrUWQWvX8njoIT+GDlUwYYKMCRNk9OlDsHIlhzvv9IMQIBz24777mjB0qL3zskIusjN23xg6\nVMXKlQlMmBDCMccoOP98Cd9/z6KxETjhhMKP7xRUgkmnGWzcyGL/fgbRKAOeT2PsWI14zV08gNwf\n0mxNKd2OGa4k0q1oeSEXsdHC3I2NjWAYJkNGKEc4F2Dvq04I0Qtpy7KMqqoq1NZWIRhkMHAgQdeu\nBFOnqvj73zM1XDece927E9TXt3z46dhGKUYURYTDYUQiEZxwglbQ5fPPGdTWqpBlYPTowhMazI4f\nn88HjuNaOH5SqVSL1FmaPGCFQABYtCiFI48kSKUY9OpF8NprSbz8chJbtsRw1VVpfPklizPPDGHo\n0DDmzAlAlglqagiSSQY33FCFAsp7FITqamDOnDS2b2exYIEPnTsTbN3KuRI+Z/c5GT9ewcGDDI44\nQsWYMQrOOEPGyScreOghH/r1q8JVV1Vj5UoB5nLL9P6Zi4xTeQloTt6h2ZJ27l8+GIm7UjoBAxVq\n6eaSF8yWrVlGKFcMbb4vuVYfQXsQGUbrFyUIAjiOx5w5PM49V8F997mTD5ttLt27a2FjVlAUBU1N\nTSCEtGgJ37EjcPzxKrp3Z3HGGQoGDZLRrZsA2cUEN7uOH6skAWpFE0Lwox9pMbKEaBY+RSgETJig\nYMIEBaoq4plnBNx1lw/btnH44gvgpJNk7N/PoqEBeWNo3TlfIBQiuPlmEWvWcIa5lq9QzZgxCgSB\n4N13OQSDwKRJMnr3Jvj97zXt/6WXVDz2WBBXXsnjJz+RMXWqjJ/8REYopFnJn33GYutWFpEIwciR\nCmprW/a3o+n0HMdl7W9nZ1VjXkWIoohAIFD0dSoHKpJ0KTLjT5vJ1tgFIdc+TuAWURNC9JbmQHN9\nBBp7u2ABi61bGcyfn5vBCnHumR/U7t2B+vrM7ejykJDcLeHXrWNx3nkqxo1rfmnKgXxJAnT+gKZB\n2ukuy7LAqafKeOIJAUceqeLDDzlEoxoJlqvlVl0dQa9eBN99x2DwYBU7d7IYPlyBqXdjSaHJTipO\nOKGlea8Vn0/gyisJ9u3jsGwZj7/+VcCVVwZw+ukyBg5UoChaYfrdu1ls28Zg+nS5RVcQqvlTIwOw\nn+SRqyklHbsSUJGka7R0VVVFIpFo0XIm175uW612j5OrPgLDaLGyN9/M4/XXJdj5aBfrlOjWTZMX\nCAFUtTmEjuM48DyfsyX8ypUsHntMbhMPujFJgFpVsVgMoVAoa0NDMxEfdRSLOXPSePRRH4JBIBZj\nMW9etIV+WSqwLDB2rIItW1hEo1p6cO/exDV5wc1Qr44dCWbPljB7toR9+2i3ZwHff8+iro5gxow0\nwmEGu3Yx6NOHtBjDiiztJHnQaBfj/t999x0OHDhQUY1GK5J0gWYtVFEU8DxvOxrBqFcWEkNrF8Z9\n7NRHSKc1zFEyAAAgAElEQVSBX/4yjJtvljFgQP5jOX2JrM4hmdT+u3SpiBEjEqir0+SYVCqVc/xt\n24D//pfBiBEEVuHCha4m3AbDMFmdPlYv8jnnsBg9mse994ZwxBEqRoyQUMwrQi01u+B5YMCATCvT\n6Rilhvm5OOIIghkzJKRSWhGlBQt8ePddHmecIRcd/ZFrVUNlpY8//hj33HMPvv76awwePBjHHnss\nzjnnHJx99tn6PrNmzcLy5cvRsWNHbNy40fJYV199NVauXIlQKIQFCxZg2LBhxU0+B9rO3XQAQgii\n0aj+73A4bPvBdMtqtQNVVRGLxTIKnWfLIrvnngA6d1Zx6aX2PDfFEts33wC/+hUDhiF44IEg7r33\nCBAS0pffucb+xz9YjBungmv9GuiOYZWtZXT6dO9OMHasiP/8h8PevWqbqncbjQL79jGOs+pKkYxg\nxpgxKmIxBlOmaEkfPK+tpKzGKcbqphYx/aCeffbZeO+99zB69GgsXLgQZ555JsKmMJ+ZM2di1apV\nWcdcsWIFtm7dii1btuDJJ5/EFVdcUfD87KAiLV0ajUC/eIXsX0pLlzoI0um0/lLnOtY77zBYtEjA\n2283gGFcigszwRyR8MQTLFSVR10dUFvL4rvvtKSMSZPyk/7KlSwuuEDNGLeSYVze7tnDYPFiP9at\n4/GLX7TH+edLuPTSFFS1pdPHHJNaKqlFVYHrrvPjxRc1DXTECAXPPpt0VCvD7WQEMwYOVBEOa3LV\ntm0SNm1iYVWUzy2pw4jGxkbU1dVh6NChGGoR45evIeXSpUsxY8YMAMCoUaPQ0NCAPXv2oFOnTq7O\nk6IiLV0AGZZtqaMR7G5Pw9Si0SgYhkEgEGhRx8GM/fuB2bMFPPpoCkccYT8+yek50GiJxsZGiKKI\nRMKPSERA587NacENDfnHSSaB997TLN22gKYm4JtvmBZhTIXi3nt92L9fsyarqlQsXuzHhx8GWoRC\n8TyfEYFCQ9lo+J9bVjEhBIsW+bFkCQ+WJeA4gg8/5HDLLdn19lLAzrn07k0wZoyKO+8UsWoVj61b\nS6f3G8m7sbExbwHzXNixYwd69Oih/7t79+6oN3uYXUTFki6QWVPX6X5uRjBQZ54xJthO7V1CgDlz\neJx9toKf/KR07RKoA0kURYRCIUQiEZx0ErBvHzBoEEH79lr9hSFDtPPLdX3+9S8Wxx5LctbiLRde\nfJHHoEFhnH66Fmf70UfFP85ffMGitlb7EO3erRXI2bYtsxgQy7It5AljpIc5bbbYmNT167W+cem0\nllmXSpGcxXsIyaxH7JZ1aXeMujrgiisk/PGPLT8Mbjr1jKRbbIyu+b6U0kFckfICUNqaulbbW0kS\nNC02lUpZRk7kO8b8+Sy+/prBggVySaxvowOPZVmdGADg3HNVJBLAG2+w8PuBa65RMHhw/uOvXMni\nzDObrdzWkhe+/ZbB9df79eSFpibgoouC+PTTuO3sOCv06aPi0081sbpLFxVNTZylNmmGVUyxOW02\nW+HxfIWAjjpKAc8LYFmA47TEjfp6BsOGhTFpkoxJk2SMHKmlcj/9NI9nntEK/EybJuHqq+0XRsoF\np2R5ySVpDB1ahXvuEdCnD8Eppyjo3Nm958T4zBVbS7dbt27Yvn27/u/6+np069atqPnlQkVbukB5\nkx3oPoQQJJNJNDQ0QFEUVFdXt3Dm5XtAv/pKCw979ln5UOUv90iXlqs0OvDMhMDzwOzZKl54QcaC\nBTLGjrWTPaeR7oQJLaUF4/mWg4i3bGHB89qctFRVIB5nsG9fcRbK73+fRjisnV8sxmDSJAljxuRe\nhRACbN7M4v33OXz7bXar2JypxTD2CgFdcomIY45RwbJa6nD37gQffBDHM88kEQwS/OY3fhx9dBhT\npgRx//0+BAIEkQjBCy8IeP55vlU+ip9/zmL0aBlvvimAEOC113g0Nrqr6dJxik0Bnjx5Mp599lkA\nwLp161BbW1syPReoYEuXopykS0OMkskkeJ5HJBLJGseZ6xjpNDBjBo9bbrEXHmYXqqrqRX3M2XhO\nHvRsc//yS63y2cCBpX2J7dyfnj21GNZ4nEEgoDWz9PkI2rUrPnb5ggtk8DzBAw/E0KePgFyXjhBg\n4UIeb73Fg2G0f8+aJeHkk62JOl9Mqrl+AaAVOXrppUZ8/LEPsszguONU1NQwAAiGDEnjxhvT+PZb\nBlddFcDBgyx27qTJFirWrOFx3nnu1zrIh61bWZx9tozf/jaAffsYCALB/v0MunQpjbyQy9LN15By\n4sSJWLFiBfr27YtwOIz58+cXPb9cqFjSLae8QF+IWCwGjuNykq2dY9x2G4euXUlGeFgxlm4+mcMt\nUGmhDeRDoH9/Fccdp2D1al5PQ3388ZSlx9wpPv2Uw2mnSejRI/+51tczePttHl26qNi/n0G7dgTP\nPSdg1ChnGn22mFSaVceyBCNHJvWqXslkZnJH794sxo2TsW8fA4YBtm1jIYoMjjjCPV+BE7IMhbQP\n0A03iOjaVSuY78a9oTCSbjQaRdeuXbNum68hJQDMmzfPtbnlgycv5ICx4AshBIFAANXV1baylLId\n4623GCxaxOHxx1t2hnV6Hp99xmLKFBYnncTi97/3A2gpcxQyfrZtV61qKS20lqb71VcsNmzg8Oqr\nCSxYkMLatQmccYY7BPPhhxyOO85eIYlkkgHLao7IJUsE8Lz2AUilio+PpVYxzQ40Ou2sCgGdc04U\nHTsqh0gZCAQILr9ccmVJ73SMk09W0NSkWbh79jDo21dF167utJI3X9dKqjAGHCaWrttlF42hQIBW\nHyGVSmVYIXZgPsb+/cAvfyngiSckdOjgaKgW49bXK/j1ryNQFCAc5vDOOzxuvFHFo4+WJgqisRH4\n6CMGp51mfa1pDHA5oCjAnDkB/OEPaZx8sruha6mUpkcOGWKPdLt0UREKESQSDCIRgs2bWQwdqiIS\n0cLrSkF22QoBBYMqFiyIYd06Dg88EMSUKQm0a5eCopCMfUsZU0zRpQvBuedK2L9fs3C7diWu1ye2\nKy+0NVQs6VIUWqoxG1Eb6yMYLQpagMPJMQDg4EFg/XoWgQDBX/7CYdo0BePGtRzHrsUoSRISiQQ+\n/liAJDEIhzls387gqKMI3nuPRTqtWC7jirFId+0CXnmFxYgR1jV9VVVFU1MTVFXNIAMaNeH2C/7k\nkwI4juCSS9zxzBvx6acs+vZV9eVxPkQiwG9/m8bTTwuoq9Pq9V51VbrsEgyVJ9q353DWWcD27Qpe\nfbUKgwbxGDgwjkgELYqO5ysEZEQhFmpNDVBTk7/2glOYx/As3TLBbUs3X32EQjTXLVtYTJ7sgyxr\nfc84Dnj66dwZdNkeSvP82rcXoKqALGvpoem05tl2oz6L8VzXrmXw0UcMXnxRC536+mugT5/mOSUS\nCf0DRfeh3ndRFLMWmClUc/72Wwb33uvD668nStLZ4cMPOYwc6SyxoWdPgttuS6NfPxWvv863egzz\njh0Mli0TsG4dh9tvD6BTJx6PPCKifXvtghGi9Qi0UwjIbd9AKeSFSqqlC3iarl43Nl99hEJI95pr\nqnHwoBbSJIpau5j5860limwPolX4l9/vx4knAsceK6OxUSPeeBy44QYlKxEVcp327wc+/phB587A\n5s0MTjpJxdtvs5BlLfOuqakJPM+DZVn4fD7deqI/NEQqEAiA4zhdtkkkEojH40gmkxBF0XavLUKA\nq68O4Npr0+jXz/2iQIBGuiNGKAXtf8IJCtatc7fjcCEk9fTTAtLp5majO3dyWLKkuUYk1YnNRcdz\n3Sd6jwrtieY2jNckFot58kI5UUz4VywWgyRJCAQCqKqqyvlwF0K633/PgZDmMdNpBlu35t6HvmTm\nYuy1tbUZ8xME4P77G/HSS+3xpz/xeOopGSNGuPMy0HmIopaZ9f772gejWzeC7dtl7NsXRXW1FpKm\nqirkHNXLs3nljSFS5mUvPX+zPPHMMwKiUQZz5uSWFYohhQ8/5HDDDYXVUzzqKK3q2vbtDHr2tKdx\np1LAs88K+OYbBscfr2Lq1EwHayHYu5eB30/g8xHs2sWgWzcV+/axALLr/fQ+qSqHjz7SalB06aJF\niDBMc2IHrftsTu6wIyO5RdbmDxGtNFgpqJyZmlBoyBiNZVUUJWviQLbjOZUxhg6V8NZbLGRZm2so\nRDBqVPa50mOIopg3/IthNM/w6acTPP44yUu4hXyc6uq0GNGlS1lceqmIbdtEdOnCoEOHap1EnV4T\nOhcaq2p8Weiyl5J4KpXSX/Ddu3nMnRvGq6/GDlU3c1c0/fBDFkuX8ti9m4HPVxg5MEyztduzp3zo\nd9nnKcvAlClBfPqp1pZn4ULgo4/SuOuu4oronnSSgo0bBdTVUbkHGDkyv4OVEK1p52efsQiFgPXr\nOdTXszjnHO2e0Ops9INJP4xOCwG5qem2BavbKSpaXjCm5+aDuT4CwzB5W6+bj+XU0r3//ga0a6el\nbgoCwfTpCn7xi9ze/6ampqxZbtYgJevjFQgQqKoCjlPRpYuEPn18mDzZn2G1muOFiwFd9lLJIhwO\nIxwOw+fz4/rrI5g1K4m+fbXsrXg87kqfLQBYs4bDH/7gx+uv8/D7Ca69NoAdOwp7NSjp2sHatRw2\nb+agKBphSxLwt7/5EItpfy/0nM49V8a0aTJEUQtpu/jiuK3aHvv3M9i8mUWPHgRffMGiWzeC//yH\na1EIiVrFgiA4LgRUzHnlQqmjMdxExVq6FPnI0JilRS1HQBPfnaIQS1EQCJYtkzF8OIGV1m/Uzwgh\nCIfDOTs2UNCHjGXtedmdWOp024MHm/Dww3V46KE0zjxTyPtgl+LBZxgGL73kx44dHP7+dxE+XyjD\n0spW0wCA/rt881q8mEc4THDwoNb6PB4H3n5bwEUXObc4R41SsGiRvR47iQTAsgSyrGX6+f0a+Yoi\ng6qq5pvq9LryPPCrX2kEV1/P4vzzU7BTMpQ+R/X1DNau5TB6tKL/Pp+2nCvTzui0A4B4PN7CWUfr\nT9iB2dKtJMIFKpx06Y22IhNjlpYgCKiubl4SU6vSyQ0r5Ma+844fnTsT/PjH1qxIyRbQCrEnk0nH\n3mKGIVBV9x46uiIghODll8Po2pXFmWcaGyVmm0dpai/s3cvgD3/w43//N6mHwtnRiQHo+qNxuWtV\nXIZamjU1WkhcMqn9rhAMHarim2+0ljv5wrpHjFBAiGbhUsIdNEgpOpWZIhzWiN0u2rcn6NOHYP58\nHv37q9i+ncXAgQpqa1Fw6Uy6egGgR7QEg8G8H81chYCM720sFmtRtLyto6JJF0CLL57mABL1+ghG\nsjXv4/Q4Tonk738P4n/+p6WTKVt4Gk3GcDInhoEtecFOQgiNJvD5fEinGfzpT0EsWCDlJNxSZ6T9\n7nd+XHihhOHDc5+k0dLiOA6SJCEcDme83LSwvNERxHEcfvYzBvfdF0A4rOmsPh9w6qmFxQD7fMCQ\nIQr+/W8Oo0bl/qhXVQHt2wPt26tIpxkMH67gwQdTrsX4BgJEJ0s7zzzLAuedJ2HuXB+uvDKNkSNV\njB6tgGHcja/N99HMda/MTSmLraXbGjgsSJdau5TM7NRHMEYK2D2OE3LZvRtYs8aHv/0tCXqZFaW5\n+WMwGGwRMVGIbswwpKgQJeNHiq4IGIbB44/zGDBAxUkntY6jghCCZct4bNzI4fHHC69QnquNOy0u\nM3p0HKoq4s03AwiFgPPOS6NnTxmEFKfrjhqVe7s77/RjwAAFCxdaE22xRBcMapquE2zYoMVjX3+9\n+4knuZBLnrAqBMQwDJ599ll8++23EEVRL8dodb1WrVqFa665Boqi4JJLLsENN9yQ8fd33nkHP/vZ\nz3DUUUcBAKZNm4abbrqpZOda0aRrJKmmpibd+SLY6FtdCME52f655zhMmiSiqopkhH8FAoGs7XsK\nC38rzNI1asksy2Z8pJqaVPz5z2G88krpCqvnm+vBgwx++1s/FixI6fGmbo5vtrTGjycYN04yNKxU\ndM+809Y8o0YpeOyx3NVd3nuPw+LFPN5/P+GaZWsGtXSdEPdLL/GYNq3l6qw16jcA1veK1inu2bMn\nPv30U3z22Wc47rjjIMsyFi1ahHHjxunbKoqCq666Cm+88Qa6deuGkSNHYvLkyRgwYEDGcU499VQs\nXbq0qPOzi4om3XQ6jVgsdijvPJi16aMVSkm6hGhJEI880gRR1MpB2q3+VQ5Ll2aSEUIyUp0pHn+c\nx8iRaQwb1nrBLbfeGsHkyTJOPLE8xG9+uemylud53RFk7BycS3scNUrBJZdwyBa+HI0CV1wRwJ//\nnMIRR5RuJREKAYmEfZKTJOCVV3i8/bYDIdgB3HR6sSyLM844A4qioG/fvrjxxhuxe/fuFvru+vXr\n0bdvX/Tu3RsAcP755+PVV19tQbrlDD2r6JAxQCtGQ18WJze0lKS7ejWDQEDFscemHIV/FaY1O4te\noNltPp8P1dXVLdKdGxuBP/+Zx+9+F7N5fPc13Tfe4LF2rQ+33iq6Oq5TUHLN1jmYYVoWIRdFEZGI\nhK5dVWzebG3T/P73AYwdK+PMM/MVRy+cpPbvB1as4LBlC4uXX/bDTv/W1as5HHkkQe/eLe9nW4sS\noHOJRqN6RFLnzp0RMXXrtOp/tmPHjhZjvf/++xgyZAgmTpyIzZs3l3TuFW3p+v1+yLIMURSLSgV2\nglwPHyFaKcgnnxTwi1+kIQhaqqXd6mSl0nRpwgVNCDFntxnx8MMcxo9X0a+fXPCLVgwRR6PAtdeG\n8MADjaiqsv94/vvfLF58UWtpM316Wq8P4TbshEaJooJevSQ880wIipLCsGEqfD7NIl6xwoc1azis\nWRMvzQShRSzcdpsfW7ZoDTuXLAmhsZHgqqty67QvvSRg2rTSabluEbfZkda5c+es29o53vDhw7F9\n+3aEQiGsXLkSU6ZMwVdffVX0PLOh4i1doLQ1dY3b59pHkiR88UUMb7yh4I03Apg5U3BcChJwLi/k\nSo6gYXONjY36cjkUCmV9EP/7X+DxxznceKO9soalwK23+jF2rIRTTrEfI7tuHYdLLw3iH//gsHIl\nh5kzg/jss/LaEzQ0ShB8WLUqDEI4fPyx9v8rVgShKCrq6yVce60fDz98EDzfnNhRSFZfLmzdStN4\ntQ9yt25asfdcYV+iCCxfzmPqVGs9ty3BSLrRaDRn9IK5/9n27dvRvXv3jG0ikQhCoRAAYMKECZAk\nCQcOHCjBzDVUNOnSC18O0qWwqk4WjUYxdy4wenR7XHhhLVIprRpWIcTuFFbJEdTibmxshCRJiEQi\neVvBA8D993M491wVRx7peBquvJjvvsth5Uoed9zhLFph/nz+UDyuVpNCkoAXX3TZ+2YTBw8CGzdy\nGDZMwbZtHHr0ADZtEpBIBHD99XW48EIZJ5/M2SoARL30TsFxACEMgkFg8GBavCd3Pds33uAxaJCC\nrl1zp6kXg9aoMDZixAhs2bIF3333HdLpNBYvXozJkydnbLNnzx593PXr14MQgnbt2hU1z1yoaHmB\nwkm2lXGfQoiawhj+tXFjFZ54Iox0ujnQ/uc/F7Bxo3PSdXIeVF4w7mJ2ktHuv8ZGh1bYuRN49lkO\nH36Y1se285IU+xK9+SaHuXP9iMWAAwcYPPJICjU1WrEdu1BVLV65sVErmO3zkYKTG4oFvcSRCIEk\nMSBEex4WLxawYweD555LOyoABACiKNquewsA/fqp6NtXwVdfcaitVbFjB4dzz5VztstZssQ6aoHO\nz21ZoFgYLd1cpMvzPObNm4fx48dDURTMnj0bAwYMwBNPPAFA65G2ZMkSPPbYY+B5HqFQCIsWLXJl\njlnnVNLRSwx64VmW1bOQnOxbCOnS9ijG8K9vvuEO/V176TgO2LNHqyrG885ItNCQMeNHgJKtkwf8\nj3/kMWOGghytplzH+vUspk8P6rGkHEfw+ecszjzT2TjTp0tYv96vE3UgAEyblgKQP53abdTVAX37\nqli7ViO83bsZ1NUR3HefD8uXJ7MWmKc6sTG2XJIkSJKkP3fZkgXMROzzATffnMaqVTx27SLo21fE\nuHHZX/V4HHj9dR733de6jks7MBO3nQLmEyZMwIQJEzJ+d9lll+n/P2fOHMyZM8fdieZARZMuRaEE\n6sSqpJZILBZr0Wn3Rz9qPrYgaMTbsaOW2llKi0t7GVUoCsHBg1GEw9ljgOk5WOGbb4CXXmLxySfF\nVbdyisWLhYzgfUVh8MwzPlx7rTN5YexYBQ89JOLKKwPo3VvBLbekMXiwhNYgXS2rS8auXUD79gzG\njFFw770+XHddGscc43w1xjCMvloBcicLGEk4GORw9tmyHubGstlf9VWreIwcqWQNX2uLli5FpRUw\nBypc06UopaZLnVENh0ot0YpKRu/1iScSzJmjtckJh4HqamDRIqnEscAEH38MzJzpRyzG4IorOqC+\nPruTLNfDftddPC6/XMERRxQ2F0CLmaaFhagXP9/+gQABy2ZuU2hZxbFjFRx9tIrrrpNwwgnFf+mK\nIQi/H+jQgeC44yR89hkHvx+48krnUQFWc6DkmqsAeTqd1iuxUdLNVSj+pZd4nHNOeTPQCoX5miST\nSd0JVimoaEu3lI40+vAa04pTqVTWF/G22xTMnKngv/9lcPTRBNXVmrzgNulS58uuXUnccUcNeF6b\nTyIB/P73HF54Ibt2Zxx7/35gwwatqPqKFSy+/LIwK9dY+9ZY2ATQqknlKjYze7aEBQt8aGoiABgE\ngwQ33li4td3UxKC6um142uvrWfC8ir/8RcA775SmtRCFVdYW1YlpMRmjTmy8J7EYh3/9i8ejj2Zf\nXbQlS9dqDLdbCpUaFU26gLOauub9rPahpEaLzxjTivMdp1cvoFevzJJ8bobbGKuSNTRUQVG0hpeK\nomUfNTQA+/bBUpc1J0AsXcqC54GFCzmMGaNi505nhEWrkdEaqeFwWNcctQ4EKgKBgL4cNhcw4TgO\nPXqwePttFddeG0R9PYcHH0zhxz9WCpZkolFtldEW8P33LFavFnDnnSJ69iz8Gdiwgcdbb/kQDBKc\nc45seyyjTswwjF4u1FwA6JVXWJx4ogifL4FUKrtO7AbcJt22FspmFxVPuoB7lq650aI5PbaQ47hh\n6SqKgkQikVGVrGNHBoQw4HmCo44i+OADFv37E+RrFUUIwc6dLGQZiMWAb75hcPnlCjZtAvr3zz8X\nKrekUild225sbMxYBVAipimzDMNAEAS9QpRRl+zWLY3Jk1V8/LEfo0cnIEnNzROdoqlJa4PeFrBm\nDYeePRVccEHhMc/vvsvh+uvDh+4FsHSpgIULk+jevfBzNBcAeu21IM4/X2sJZdaJ6bb0XhRDmqUk\nyLaUKWcHlWWXW8ANS9fY/JHWSLDy/pc67tY8vqpqDSCj0Sh4ns9omtmlC/DLX4o4eFALgk8kgMGD\nVcsW6ea5aHGcmv581VUKeD5/F2Fj7K8sy4hEIvqLSmOAzdYsPRd6bKPWS5MJAoEABMEHnmdbeOkT\niQRSqZSuS+a79m2FdN97j8OuXSxuvLEJxfDBX/8agCBo9ysUImhqApYudWYnWRGlJGmRI3/7G481\na7Ri5VY6sd/v13ViRVF0nZgWb7LbUNQINy1d2kKo0nBYWbpOvsQ0eiEej+et/mVEqRxj5mNQazJX\noZypU2UMGyahqSkEUVRw/vkCLrlExYABuY/Zo4dmEUejGmE3NgJjxmTuY5y7OfaXygf071QzpDUK\nKHkaf1RVtVy6aktdAqCZiDmOgyiKOa0v81iKounaVVWOL7cr+OgjFv/4h9bu529/84HngYEDi3Po\nSRIDjiPYtYtFTY0WHminhkI+rFnD4vPPWWzaxKF/fxVvvsnj3HNlBALN2xh1Yvpu+Xw+nYCdFABy\nE+YU4ErqAkxx2JAuYH/5Q0mNwklzSqfzckq6qqqisbERHMdZFmA3j9+9u6y3drn9dhkzZvB4910J\nVh1/6HyCQQaTJxNs3aq9xD16EHTsaD0X2jGZyhpUDwSgxyzTkprGudJeZxT0ZTUTMcuyUFWfbhFS\nTZ1+FCkRG6+9OYmAOoRCoTAIUQCXm1bmw/LlPGbNChwqo6jFCft8yGi5UwjOPjuF//f/Qkintc4N\n1dXA+PHOiJyuOihUFfjyS61m7r338pg6VUIiwWDfPiavbGHUiY3jU7kom3ZfCn0YyCx2U0moeNI1\nJkjY8fwbC3YDyFmLwOpYpbJ0za17CqkJPHOmin/+k8XNN3P4059yv5xaimj2uEy6CggEAqipqdFJ\nkx4zlUrpzjIzKWabazYi1vYl+stKrSz6//TFNu5nPuaePUB1tXZ/6UeBfhCc9uByiuuv9+vxxoQA\nokjQvn3xMsfUqWmwLIObbgqjTx+C224TMWBAcXUaGAbgeYKNG7Uu1SecoOLAASanvJTLmLEiYrqP\nsTcajWih9ZsLJWTjR8SzdFsZuQjOKvyL53kcPHjQsSRRqFyQ7RhmJ1k8Hs/Z8SL3/IBHHpFx/PE+\njBunYty47JJBtnnS60SlBEEQ9BeGfrQkSYLf73ec9dZyvrR/lka6HMfpXnb6olLd0JwGaybiaJRF\nJKJVnqMfBSqDGHtwWYWwFYumpsx/KwoQDhdPugyjEe+zzwbxu9+lcdxxxRfGYRitRfu55/oxaZKE\nffsY9OmjomNHd7Xw5nurga6K/H6/Zet2432x6vZBYZYXKi0xAjgMSDdXrK4x/IthmBYWZCktV+Pc\nrGDuJkFb98Tj9kv+Wc2nfXvgr3+VMHu2gA8+SKNDB3tjmXVb6sCiIWC0hKYgCKiqqnIlNpJeA1H0\nQxD8GasOs0VstJhkWT6UZdX8gsZiLCKRZuuKguM4fSyqNVM9kkZcmDVicx+ufDjrLAUvvshAkrTt\nBQHo318pmtCpVZdIMAgECiNFq/P49lst1PDyyyVUVREceWTuYjhmiaLQedDrTNu3098bw9isdGKj\nRXjU5JkAACAASURBVJz5ofXkhVaF+YbkC/+y2sfpMZzsY9SdcznJzNsXgrFjCS64QMFll/F46SVZ\n10ut5m+MtzXqtoIg6ERLSYzn+UMarL3W5tlALWZRTOOJJ6rxpz8FIcsAx7G47z6pxVI3W/C/kYgb\nGlSEw1rkA52r8RyNMNY4MBOxUUKhH518LXoefDCF1avD2LcPqKsjGDFCwaBB7uV/J5NaHLYbUBRg\n7lw/7rhDxNCh7paULATZ7q3x/poLAEmShDfffBPbtm1DXV1da029YFR8yBgFfVEURUFTUxOampr0\nONJsy+Byki4lGhpyZbebRKHzufVWBbt3M3jySRaEWJd/TCQSaGxsBMuyqKmpyZASqKYKQO+YIAgC\nVFVFKpVCNBpFU1OT3i1BlmVbmno6nUZTUxNUVcXKlbW4//4gRFEjg+ef53H33fbsAPqyah71ADZu\n9EGWAVnWwp5YloUsy5YhZ/S6Gdu1sywLn8+HYDCYUQaTrpQSiUSL0ot0vAMHGMRiDL7+OoYvv4wj\nHAa6dXOP0JJJ9yzdxYt51NSQvF0rco3hxjxygVq5PM9ndOygrXhYlsXSpUvx5JNP4pprrsHIkSNx\n6aWXIhqNZoyzatUq9O/fH/369cO9995reayrr74a/fr1w5AhQ7Bhw4aiztEuKt7SNd5IURSRSCQy\nluv59i0H6UqSBPFQCax8TrJCdWMzfD5g/nwZJ58s4K23FEQiDE480Y8LLlChKM0t6iORiG7tUVCi\nyqXbmi3N5sIqzUtIYxslWgWNyhc8z2P5ch7JZHN1NlpI+5Zb7CUUEELw3//KuPZaHz75REA8zuLa\na2vw0ENp0HKo2aQJo45IrVj6cTSGHxqzEa0sL47j8PTTYUyZIiIcVkEIgx07GJdJ1x1LVxSBu+/2\n46mn3GvzbhduEDfd3+fz4ZFHHsFdd92FE088ER06dMCGDRsy+qPZaUi5YsUKbN26FVu2bMEHH3yA\nK664AuvWrStqjnZQ8aRLl8fpdFpPILBrPZaadOmLnkwmdbJ180OQb9tkEujXj2DNGg7Tpyt4/XUf\nwuEUJkyQdCvbvGyzq9vaWfJTIqYQBAGC0NxRo2NHAp5v7mZMCGw3aqRW7PPPB7Frl4BgUMva2rGD\nwcKFPK6+WrY9TyMRU3I1arvGa8wwjP4h0iQIgueeC+D55xuwe7eM1asFbN7M4NtvVQwfTuD3Fy7F\nUKJKpbS6FMXi6acFDBigOm726QZhugHzs97U1IQuXbrg+OOPx0knnZTxNzsNKZcuXYoZM2YAAEaN\nGoWGhgbs2bMHnTp1Kul5HBbyAiFEz55xslwvJhohF4yZZDSG1Ymn3y3S3byZwYABKurqCNauJair\nU7Fpk1+3bqk2S7ONJEnSGy8WInsYl/xapplmJdLMM4ZhIIqiLk1ceWUTamsJBEHLiAuHgbvvzh39\nT51viUQCPp8PBw74EQwCBw8yOHgQCAYJdu/O/2Gj8wwGg6iqqspIjKEpy+l0uoV0Qj9UVAdeuZJH\njx4qjjmGwyuvVKO+PojGRhbbtnFYtcrXomml3ew6IxIJFNyGnhJmUxNw//2t1+zTTeKm4+RypNlp\nSGm1TX19vStzzIWKt3Q5jkM4HNYfaCcoNBoh2wNkdJJRPdlJNILxGG6gro4glVIxdmwKhHCIxVi0\na6dClhXdiksmk47ibe2ARgawLIuqqqoWCR7U0uzRQ8Hq1Q149VUBoqhi/HgJRx0FiGKmNEH3oSQo\nCIL+4Rg+nOC99xgMH67izTc57N7NYMYMZ3ollVPoh8J4TKukDgC6jDJ/fgCzZonYswc4eJCgulqF\nzwf07q1i61b+kFac6RiyW4wc0LRuSYJlsosTPPKID6edpmDQIOeyR7k1XbtjNDY2ZnWkFWrglMOi\nr3jSpSjGyeV0HzMoISQSCQiCkJFJVkoJI9u2dD7HHpvAv/5VjZ07AwAIwmEFEyfGEYs1e+gFQdDb\n2Bf7wFErNB+JG5f8PXoAV10FEMJCVVlLjdholVM9mGLyZAX19QxeeYVD9+4EsRiDqVPzky6NTqAl\nKa3kFCoL0DbsdD9KxF99RfDZZxzOPDOGgwdZKAoPnieYPj0JSVLBcQIASZdPaOwqTacFkHG+5jRn\nrbOwZuUWemsIIThwgMXjjwt4++1EYYO4BLcJLZela6chpXmb+vp6dOvWzdU5WqHiSTdXnK6dfYvt\nrUYzyRiG0ZMucm3vdHw7MFoAxsy2jh2rcNNNwOefp6EoQJ8+KkIhAalUs8OLRjFQp5Dxx262EI3M\nSKfT8Pl8jrL8jOdt1l5pQL2xsEk8Hm/hrLvqKoLLL5ehqsDYsQG8+iqHadOyE6+VU8/JPCkRL1wo\n4KKLFHToUI26Oq3mxWefAe3aKdi1i8PppyegdWvOTOow6tyUiI0ZldQiJoSgoSGNQIDoK4dCsrge\neMCHadNkHHmke7G+bWEMWZazOqWNDSm7du2KxYsX44UXXsjYZvLkyZg3bx7OP/98rFu3DrW1tSXX\nc4HDgHQB92vq2tnHThxwqWE8HiUSWZYRDAb18C6fT8XQoc3RA5LEWC75zQ6wVCqVQcS0EI3RKWTH\nWiwERinBTOL5oibmzlXw61+HcdZZMgIBpsW4diIz7CCZBP7+dx7vvJM6lPWWwvjxMo45JoxEQkCn\nTgS9egVbSBP0mprlBHOPP57nD+nGfoRCmoyWLYsrGxETQlBfz+L55wWsX9+6Vq4bMMe7A9mtZzsN\nKSdOnIgVK1agb9++CIfDmD9/flnO47AgXaB88gKADHKjVbXcOkYhc6Lxo7QsH7Wm6FjJZFJPM862\n5KckYLQcjOSWTqczYnipA44QopO8G5BlGclk0rKIDpA/GmHMmDSOPtqHhx9WcOWVSZ2QKJG7lVH3\n8stam/WuXVOIxbSPQ21tBHV1DAC6esp+TY3X1lyBDYCemJFIQI/RtSr8o0kQYtZ02gceiGD2bAmd\nOhVu5boBtyxdM3KNma8hJQDMmzevqDkVgsOCdMth6RrJiyYT2K1MVgrSpSQCaFZS9aGWCVbxtoUu\n+a2ImIZqybKsnz9NZzbH5jo5Hk26oB8zJ049MxH/8Y8qzjgjgpkzGVRXS7rjCmiWLMxztYNEArjy\nSh+WL+cgisCFFyb1iA+7dV2t6goYnXXGMDstBZggFGomLLMcli3NOZ1O48svWbz+ehXWrj2AdDp/\nXYNcKAVhFgKzpVuJOCxIFyi8pq6dLCpjZTJBEODz+RzFAjvVjfPBqCMzDINAIJBxDEqMpVryC4KQ\nEWKVyyI2a8RW4xr1YBqVUAz69yeYOlXGPfcIuO66BGpq/IhEfC3mmi+hw4zf/lbA8uUc0mktqmDJ\nkiDOPpvFT35SbOUvRg9Ro3IVlRfSacDvb+51ZpYmzHHEVL7geR733RfCFVfE0b49Z1nXwDxWNiex\n26FehcI4l1QqhWChcXStjMOKdAvZJxvpUr2SFsuhTrJ4PF6Qo8uNOZkrkgmCgGg0ikQioRNFOp0G\nwzCOrK98MIaAWY1rtojNoVbZiJgSLg0tc+vjIMsyLr44hjPO6ICdO9uhfXvgoosUnHiiajuhw0zE\nLMvi9dc5SJKWzMFxWo3bN97giiJdoy5uDIUDNDlBklhUVbGIRCKW4Ws04sGs6/773ww++ojDn/+c\nACFCC4vYmF2Xq8CMW3DbkdbQ0FCRxW6Aw4R0jREM1Bqwu58VweVykhUa2+tke7NlTJfeVLcNh8P6\nCxgOhzMK09D5mZfQhRAaXYo7jePNF2olSVJG7VwAGfJEoeRrnO/rr9fgmGNUfPUVh4kTFTz9NI8e\nPaQWhbrzacT0o6GqKmpr/di7lwfHaSFcHKe1Wi8U+aIoZBnYsoWBJGmOu2Awe03izB8Vc+e2w29+\n04RIpPme5ZMmrNKcqR5O5SSnshGdo9uo1ApjwGFCuhTGsBs7MEsSzR5+KaeTrNSOMeNxjNKGlW5L\nk0J8Pp8+32KW+8bjFhMCZobRSy/LMvx+v07INA1XFMWMugh2PxrG+Wrj+vDNNzxOO03FX//Ko7ER\nALRMNTtNHY0OKUo+fr8f99+fxnnn8YcqohF06KDi5z9vQCKRaRXnu1bm+VpFUaRSwI03CnjzTQ4H\nDwKXX+7D/fenccQRLedKiXjfPoItW2SsX0+wcyePGTOa2y2ZS2FmI2JaaIb+XZZlSJKUtYW7E/3e\nDUuXPgcNDQ0VWUsXOExIt1grlDrJaF+u2trarA+IU422UEeaUbelS2+7um2+5b4VudFY0VLpwTQq\nwWpcn8+Xsa3xo0Hnmk13NWe/0XE7diSIx4EZM7TeX9u3AzU19u+DMYqCjnvKKcC776bw5pscQiGC\nKVNkhMPNYWF2NGKrca3w2mscNmxgUVVFkEwy2LWLwVNP8fj//j/rYkC7dil4+WUVLMvgsccimDBB\ngaL49Zb05uuajYjp80cJlj5z9MNgdNZla81jJvZSREB4lm4bgVOCo9s2NjZCEARbEQmFyAVO5yRJ\nUot4W2OdBKoz2w3sz7bcN76ENJMMaLaIi62bC2Qu+e3M185y36qQjtlavOQSGQ89xKOhAVAUBmee\nqaBv3/z3wRxFYQ6F69ePoF8/SnwMAPsaMf34Ues217P2/fcMfD5g924GBw5oNTS+/z576vmHHzII\nh/34+msBggAMHKhi61YtTdrudTUTMf1QC4KQQcB0PCsiNqc5m1O5i3mWjO9RpbbqAX6gpGu0JAHo\nBWncPIb5ePlAX3ZqtVVXV2fE2wLQnWjmGgGFgL6EDMPoBV0CgYBOtoV4983n7FYigpEw6NKckgGd\nbyqVyphrly4cbr2Vw969PCIRretxrsObEzIKjaIwkxsd16ixUw0+13U95hgVy5ZxGDJExerVHL79\nlsUZZ2TWFqGOXp7nIQhB+P0cVJXBxRfLGdXb7M6VzpeuHgDtA0yNAKODzWgRG0E/gMYQPfqBpNmE\n+epN5JszoJFuO1q/s8JwWJCuE3nB3JaGLuGdHMtNR5pRt6X6aSqVQjqd1h9It/VVelyrAjJW2xl7\nWlGnEn1xaKaaMYSpFFlqQGYURSQSsYx3NVpugiChc2eNiJPJ7B8NaukDcDXqw1iLIhwO52xBZP7A\nnXoqhy++IHjtNR969lTx9desXlPCyho/5hgGK1cyGDxYhaJojrfevZ0bB1bFf7JZxEbtm66krIiY\nfizN/e/M9SbyRU4YLeWmpiYceeSRjs6vreCwIF2KfOFWVk4ymu7qxjGcbm+l21IvMfVqA5pGa3Q8\nFYt8IWDm+VOioqsBoz5sjJwwOjLdsMYpKHkZQ+XyzZUiH7nRvwUCgaIbbRqPST9o2ax8O3O95JIm\nnHuuCkJY3HNPBHfcweL+++OWMc09exJMmKBg0yats++xx6otnG65QLVmjuNafCjtShOUbM2ygiRJ\nLXwh2epNWDUiNWYVGi3dSmzVAxwmpJvL0rUKtzK+AKXWaK22t4q3pbotx2kdcenHwOjhp90nnHj2\nzcctJATM6pyMoUvUESlJkq7/Uevdaq52j2le8rtVSIdmbImiqFtodHXhVEYxw66jzO5cw2GN3G69\nVcQJJ0QwdSrB8cdnkiT96dlTI18nMFq3TtK5sxEx/SDTqAejn4CSplEvNsKKiI1pzoAmsT311FPY\nv3+/Kx/I1sBhQboUxq+pedmezUnmpuWaC/SBNHYANsbb0geI1nXIZimaw8FkWW4RhWAmi1KEgNFx\nzYH9xmtstIjtFtGhyFeDoVAYPzzGJb853tWYAWbno5Ftae4GtOdBxF138fjDH9phzZokOK7wzDoK\noybsRiagUeulzzRtm2WOJQaQoetSmImYPifG67t9+3asXbsWS5YsQceOHTFu3Di9kI0ZBw4cwM9/\n/nNs27YNvXv3xosvvmgZata7d2+9JKsgCFi/fn1R1yIXGFKKyOUyg1pE1GERCAT0LC2ax58NVNy3\nm1Koqqrjpc2BAwcQDAb1DsCBQEAfi4LO3xhvawfmZR79oS8goL1cPM8X3BHCCjS0DIBej9cOrOYK\nIIPUKIk4rcGQC3ZiY632scoAM5MwXU3Ra+wW2Rq1Zu3ecZg61Y8xYxT89reZoWP5ngPjDyUwWZYd\nl7bMBTuyCt3O6toCyPiwGaMejPHdDMPgvPPOw6JFi7Dv/2/vy8Oiqvf/X8Mw7IsbSwgVLoCoqaza\nk4oW2KJmZotLcl3IfLop6teF6knN65VcUq5leUsLs6vXSs0CNEnxVjIguKapiaKsLgw4DCDDLL8/\n+H1OZw5nZs45cwZhPK/n6XkCh3Pe58w578/783q/3u/3nTuoqKjA6NGjWW1asmQJevTogSVLluCD\nDz5AbW0t0tPT23wuNDQUxcXF7ZKccwinC7QWCdy7d49ytqSSzBqIgsGD4+Q/o9GI2tpazl+OVquF\nRqOBs7MzFWGa09uSKba2gpmBJuc0l/ziA3oSR6yIjq75JC8fPXrnS6MwQY+abV14mHy2Ttfq/EhR\ngRAahe0cZIFg3uNr12QYOdINx47ds9of15wjJvYqFArWnZEQ0BcIDw8PQfeYzVZCSRiNRhQVFcHf\n3x9nz57F8uXLceXKFavvbUREBI4dO4aAgABUV1cjISEBFy9ebPO50NBQFBUVoXv37rzt5guHcLpG\noxEqlYoa12OpuIEJenaZ67lqa2vRtWtXi+cgLzp5cMg2izxEZJsLgOJXxYA5p8glarO2faarHcjM\nMzFA7pVMJqOcIj26IfYSR0x3bpZssKa5FQom1+zi4tLGYfC5t2z3wtICsX69M379VY59+5otyuDo\noCcjyU7LWkTM5fsVsoPgCp1Oh4aGBuo7/7//+z8cOnQIt2/fRmxsLOLi4vDee+9Z3HV27doVtbW1\nlK3dunWjfqajV69e8PX1hVwux5w5c5CSkiLKNbDBIThdmUwGV1dXuLm5QaPRCNL9ifV58nBrtVq4\nu7vD09MTarWa4s6cnJyoiE5M7s+aBIyZ/CJ/Q4/azHGuZIGQycRtpGMpaiYJF3PFHJZ4TACiaG7Z\nQF8g6PeCRLoE5hQe5hwx4fu5LBDz5+uwe7cz9u6VY+xYPRQKwFxgSd/1mOPz+dxbNqkdXX0jBuj8\nLaE/srKycO7cOXzxxReIjo7GqVOnUFxcDA8PDyQmJqK6urrNcVavXm3yM1vugOC3337DQw89hNu3\nbyMxMREREREYPny4KNfDhENEugAoATdfvpWs0t7e3pz/pra2tk1ijpm4o/O2ZPtMz+aKuR2lS8BI\ngYNQsOkx6fbSqQmhECtqtrR9lslkJiPfbXW6YiTK6M8CfftMeGGSg+DyLBw+7IRp01yRmKiHr68R\ns2bpqOoz+vnITs7Dw4PXc2GJI6bvRMSU2gGmqgx3d3eo1WosWbIETk5O2LRpE2+ZWEREBPLy8hAY\nGIiqqiqMGjWKlV6gY+XKlfDy8sKiRYtsuRSzcIgR7AQkYrC3GoH+N8SB3L17Fy0tLfD29qb62xIH\nSyRKcrkc3t7e8PHxoV4uUqyhVquh0WioF5vUvluCXt86Op0kD/m+WGygNzwxGAxwcXGBl5eXyQj1\n+vp6qNVqNDQ0tBlRbg06nQ4ajcZk3LvQF5bQDWSRI9fu6uoKV1dXKnokI98bGxupXQaf77ylpQX1\n9fUUTSTUyZB7S5QrZGwSkQbKZDI0NDSgvr7e6rNw7pwTAgKMuHZNBk9PYMsWBSor/0o8NTc3Q6PR\nULpbvs8F/d6SMfU+Pj5wc3MzoU/u3bsHjUZDjZjn8yzQQb6rxsZGuLm5wd3dHXl5eRg/fjwmTpyI\nL7/8UpAud/z48cjMzAQAZGZmYsKECW0+09jYiPr6egCtifWffvoJAwcO5H0urnAIeoFAqOZUqNNl\nawFJHkgmb8vMEjO3+lw6gxFnaE8ezVzhBJHSAMK2o/ZIwBFb6BVw5irV+Gz1CZiRoj2y/GxqFWvP\ngpOTHOfPK5CUpMOOHQoMGmSAwWBEebkMAQE6u1TXmUvumXsWCJ1Fv7fmdkdMbXNTUxOWLl2Kmpoa\nZGdnw8/PT7Ddy5Ytw8svv4xt27ZRkjEAqKysREpKCrKyslBdXY2JEydStkydOhVJSUmCz2kNDkMv\nkC0b29bf2t81NDTwap5x9+5dKkp1d3enas2Z/WxtcTBMoTk9mUS2pLYUODDBrPrie1xzFUrEken1\netFla3TNrTVpIJu9lhKLxJmTKNpeMjChUrt33/WAVtsaAHh7y1BRIcfChY0IDW0Sfcuv17cW83BV\nf5ijJph6cpLfII7c2dkZBQUFSEtLw/z58zFlypROWwBhCQ7ndOvq6uDt7c35Ydbr9VCr1Zy2LsSZ\nkv62RK4iht7WGkg0AIByLnTHJpQfZhZO2MNmkpyjV93xzeqbs1nMSJ8egdKlSmyKCSGLqJg2X74M\nfPihAlqtETqdAcOG3cPkyU1wchJPamdJuibkWObKhnft2gWVSoWrV6+ipqYG27dvR8+ePQWdpzPA\nYegFetabvNhcQC87NAcS9RANMKEGuPa3tQXWtuVcts4kamW+MOZ60doKS0knS1VqXBybLWW2Qmzm\nQvtYc2z2sDksDFi9Wotr13RQKLSIiFBAoWhNBttqLzkGiW7FsJlEuSRpaDT+1dXO3d0dSqUSZWVl\nqKqqwhNPPIE9e/YgNjbWpnN2VDhMpEucTn19PRVJcIE13S29KxmR8jQ0NFCJMScnJyoxYy/ej2+G\nn42WAGDi0MiIb+JgxLKZXhbs6urK6WW1VqVGbCZlz2JqbgHTclhrBSrWqqmYNtuzNJhLwQfTXmaH\nMKYjFjO6ZYJZQKHX67F+/XoolUps3boVvXr1gsFgwJ9//omgoCBeiqLOBIdxugZD67RTjUZDvfBc\noVKp0KVLF5MHl6m3ZfK2hE8kCSTye3q0JlSqJKYEjH49bFVfYtgLmL5QthZ7MPls+laULluzVQpG\n78Nry+Jjjn8HOqZ0zRLnSvIFpKud2LQNCYguXryIBQsW4IUXXsC8efMEPeMHDx5Eamoq9Ho9Zs+e\njaVLl7J+7sSJExg2bBj27NlDJczuJxyGXiAQQwJGmom7urqanUtGOFDStcxSRp8r30pPDImZJCMO\ngdAf9HJkSwoEkuywZAPdCYjFrxJqgdhNChHIYsdHwG/OZrp6wNZ+CcReYgv9OwTQZoHmay+BpfaL\nfO2lF5KQ71Cr1VKOVqvVUpG0UHuB1oWNPiwAADZv3oycnBxs3boV/fr1E3QNer0ef//735Gbm4ue\nPXsiNjYW48ePb3M8vV6PpUuX4umnnxYkZbMHJKcLUJEq4bGInpY4JgLC28rl8jZyHOaDDFjnW4kz\nBv5KwIktASMRKBv9YclewgmaS3wBMJFqtRcnzLTXXIN1cwsdPSIXW1JlbpQ6m718Fg76/RCbWqFz\nt2xd4oTYy7wf5Jm+du0a5s2bh9GjRyM3N9em6ygsLESfPn3w6KOPAgBeffVVfP/9922c7ubNmzFp\n0iScOHFC8LnEhsM4XXoija/TNRqNJo1vSJKMUAaEauDL27KV3rK1ZiSfJVtQMSAkArVWKkxPfBFw\nmffFx2ZLAyzZ7CUvvbUG6+S+6vV6UW0mxyTPhyVHbm5hplepMRcOAFQEKmY5Mxfu1pK9lhwxiZQB\nUDuUbdu2Yffu3fj4448xZMgQm+2vqKhASEgI9XNwcDAKCgrafOb777/HkSNHcOLEiQ4jP3MYp0vA\njE4tgWx9DAYDXF1d4e7ubrKlBaz3t+ULIhKXyWTUlAj6GBO2bSiXbT4Bs1jA1giU7ohdXFzabEP1\nej00Gg2Av3qfCtmG0qkVWxKSbAsHaflJ/o1QC2x8Nl/pGpdWhtbsZXNsZFdFFn5SmWaL1I7AFmWC\ntYWDPjR01apVuH37NkpKSjBgwABkZ2eLNjady3WnpqYiPT2dCsQkekFk8Il0mbwtcWpsvK3YEQaX\nyixz22Zr+lZmdlhMJQU9AhV7G2qv6jpLXcZslYIRx0X4ZnvRFES1IrSqjnlseygTyGJFSoC9vLwA\ntHbuKi0tRUhICM6fP4+goCAcPXoU8fHxNp+zZ8+eKCsro34uKytDcHCwyWeKi4vx6quvAgDu3LmD\nnJwcKBQKjB8/3ubz2wKHcboEliJdEpmQZARxpqTumyRD9Ho95HK5XSVglhy5tW0zc1or4aTpzkVM\nTthaBGptG2pu4QBAceRic8LW+FWy46CXNtMVCGTiMFPh4eTkhObmZrvIwCyVHZujqugRpiXNM9nV\nia1vBv6S3JEuZrdv38bChQsRHByMb7/9lioiIt81F1hTJlRWVuLo0aOIjIyEq6sr6uvr8cMPP5h8\n5urVq9T/z5gxA+PGjbvvDhdwMKdLnCZbpMucAkznbUm0Sx548m+klyd9my+kGonPIEhL10ZeOjod\nQRJwxCYS2dkqA7M1ArW0cNDbSJLf0+eTCbnHBEITZXQFgrkeE2ThAECNkCELtC2Ol75I8BmnRBYO\naz08yH1uHdUuXhKOLrkjzZYOHDiADz/8EOnp6Rg9erTJdRA1hzVwUSaMGTMG+/fvR2pqKmpqatDQ\n0IB+/fpRY3vmzJkj2nWKDYfR6QKgElP19fUUd0RWeNKbkwxOZOuTwHQubKJygL0JDRuYEjAxH3i2\ncTnm9Jd8ZGtMKkGsaRbk2MyCDy6FBlykc+1FUxDdNF2Pa0tpsy3tF61Br9dT46jIc2+pmIPP90yi\nW/I91tXVYfHixXBzc8PGjRt59TJhIj8/HytXrsTBgwcBgBqvs2zZMrOfX7BgAZRKpeBzticcLtKl\nk+Z03tbX15d66AisbffZokv6lrm5uZnasjEdBNmC2tsBMBt/85GtMR2xvbpqAeYjULbokrnNZ95j\nZgRvz9JgczQFH2kgmyMWIwlnyW5L3C25x2ThMEelmOPgSXKZLBI///wzVq1ahffeew9jx461+Tq4\nKBMAYP/+/UhLS0NVVRV++uknm87ZnnAop0tgNBpx9+5dkymndGdLIlAh2302LpDpiI1GI/UA5nZw\nlgAAG4BJREFU0xN7tm5BuXLCdHCRrTG3oGJVwRG7+VRQcdnm0xN15HrE7qzFZwGyJLWj860AKEdG\nnkcxk3AAtwWI3GN6qTyXZCgAkwVIo9HgnXfeQUNDA3JyctCjRw9RroHrdzhhwgRMmDABv/zyC157\n7TVcunRJlPPbGw7ldAm1AFjX24q13SfRJXloiQSMZHOZLxw9wuT6cJkbESMU9IWDbBPl8tZ+uQaD\nwWwEz7d7GR/NrSUwI3iyAJHEjEKhMPnZFg6eWa0mdFw9myMmi31LSwuVe2DmDfhu8+l226JMsCYF\nI3y20WjE888/D39/f5w6dQopKSlYtGgRZ76WC7goE+gYPnw4dDodampq2mWwpK1wKE6XdPdvaGig\nmmWQSPN+bPcJmNs5OtdqyakxJU9ilQWTY5NIjm0BYor22bhLc06NzmXz7XNrDXSagtmP1lz/Az4y\nMCG9brnazexJy+W54JIMFXPasbljy+Wt3cAaGxuxYsUKXL9+HT4+Prhw4QKuXLmCW7ducZ6ozeWc\n4eHh+PnnnxEUFIS4uDjs2rXLJJFWUlKCXr16QSaT4eTJk3jppZdQUlIiyvntDYeKdF1dXaHT6eDs\n7Iz6+nqTBtoKhUJ0PSUfCZi57RyzEom8ZIRPFHuwIj0ishTJWeOH2SJ4Irezx+LGJVFmjUqxJAOz\nVxm2pQhU6Daf/r3YqyMYW+lxcXExFi9ejNdffx2bN2+mnDsZF8UF1qRgX3/9NdauXQuDwYB+/frB\nz88Pc+fObaNM+O6777Bjxw4oFAp4eXlh9+7dolx3e8ChIt2ZM2eiqqoKUVFR8PLywrlz57BmzRqq\njZzRaGyTKBASFdijCxhg+oLSIZSWYMIeERFxEC0tLWhpaQEgbvcyse1mOjV6BzNShs2XSrG33WwK\nD2KzTCYzKfARw+ky7dbpdPjggw9w8uRJbN26FY/+/34HfKHX6xEeHm4iBWNGsPn5+YiMjISvry8O\nHjyIFStWdBpVAlc4lNM1Go04fvw43nrrLZSXl2PEiBGoqKhA3759ERsbi6FDh6J3794A/po0wYe3\ntKcEjG27b237SZ/Ma+lls1SZJYbdTAqEzanxka1ZOrY9Ijkii2M6NaEyMGaUKLbdzc3NFOdMchV8\n6B9rx6ZHzhcuXMCCBQvwyiuv4M0337RpweMrBautrcXAgQNRXl4u+JwdEQ5FL8hkMmg0Gvztb3/D\n3LlzqUGRly5dQn5+Pv7973/jwoULcHV1RVRUFGJjYxEXF4cuXbqwbvHpkaXY7QsJLG33zW0/LcmT\n6I6YmRQSm6YwR6/YIlsjn7dWUWYL6E3L6ce21uiHfm3mdkqkCEfszmvk2CQCZRvAyaQluNoMtO3H\nYDAYsGnTJuTm5mLbtm0IDw+32X6uUjCCbdu24dlnn7X5vB0NDuV0gdZKlTFjxlA/y+VyREZGIjIy\nErNmzYLRaIRGo0FRURHy8/Pxn//8Bzdv3sTDDz+MmJgYxMfHo3///iZaWxLhuLi4iEolCOmPylUC\nRs5BKBB7TQDgwpPzsZmeYLJHVC5EBsbUaDNtJs6M7EhIEY5YYONX2cCUM3KxmdhNl/RduXIFqamp\nGDNmDA4fPixaIpTP83f06FFs374dv/32myjn7khwOKdrDTKZDN7e3hg1ahRGjRoFoPXBvH79OvLz\n8/Hdd9/hvffeg1qtRkNDA0JCQvDJJ5/Az8/PpOxRyHaZgE5TiFGEQH/ZiGPR6XRUhExvSC3UZkCc\nqQVsNgOgpGqkHFgmk9ksW6PbLaTM1prNZHEg2316KXZLS4sonLatjcu56MoB4Ndff8Xu3bvh4eGB\nM2fO4LPPPhOlMQ0dXKVgZ8+eRUpKCg4ePMhpYGxng0NxumIhLS0NX375JVJSUuDr64vCwkJcv34d\nPXr0QGxsLOLj4zF48GC4uLhQCgTAeumqPUtVmVQCc6qvuWSMuS0+89j2Kg0GzCecbJGtEdhTBkZf\nhAl3y3RqQvlhrtGtEDCr4RQKBU6fPo0NGzbgzp07aGpqwoULFzB37lxs2LBBtPNykYLduHEDo0eP\nxs6dOzF06FDRzt2RIDldFpw6dQq9e/emRvUArQ/qzZs3oVQqoVQqUVRUhKamJkRERFC0RGhoqEni\nix6lEYdLatXt4bRkMhkvx8LWpwFAmw5V9tLcComcuWpxCT1krwWOrf2iNZutLXhk8aBHt2I/K6QX\nCQBqTNHXX3+NL7/8Eps2baKi2+bmZty9exf+/v6cjmtNCnbx4kXMmDEDRUVF6Nq1K3x8fDBr1iyk\npaWZSMFmz56Nffv24eGHHwbQqigpLCwU6/I7BCSnawN0Oh3Onz+P/Px8KJVKXL58GZ6enoiOjkZc\nXBxiYmKgVqvR1NREbaPYEl5CwaUwgw+YURrpj0rnN4XQEmzn4eO0rIGpliCOmPDwtkrtmOcivLAt\ni5C5BQ9ovT9iT7eg33Ny7Js3b2LBggXo1asX/vnPf8Ld3V3QsblIwW7fvo3r169j//796Nq1KxYt\nWiTKdXVGPHCcrphwdnbGoEGDMGjQILzxxhtUz4fCwkIcOXIECxcuhEqlwtixYylaIjw8HDKZzCSp\nwVfTKrQPgzXQS5r1ej0VadGjS9KSUaje2R5NdegtDomthF6hn0+oBAwQjxdm2kxoA6KoIL8nfCvf\n6jQ20MvfPT094eTkhH379uFf//oX1q5di5EjR9p0LVzmlfn5+cHPzw9ZWVmCz+MokJyuiJDJZOjS\npQuSkpKwceNGJCUlYdWqVVCr1cjPz8fOnTtx7tw5yOVyDBo0iHLEPXr0gMFgoKqPLDkHektHsZul\nWNLFWpvzZq0gQqyeBubAlIExFwGmbM1cw2+2xYO+JRf7nlvibi1Vp3FJiDIXCldXV9TW1mLRokXw\n9fVFbm6uCYUmFHylYA86JKdrJ+zbt48qjQwICEDfvn0xffp0GI2tQzCLi4uhVCqRlpaGiooKBAYG\nUrrhxx57DDKZzERrSZyBXq8XvaMW38jZUkacTe8MgJLe2XuhMJdw4iNboyslyL+JzQsD1pUJ1kqx\nzbWRJBE/4eJJdHvo0CGsWbMGK1euxDPPPCPatYh5Tx4EdHqna43Av18wV4tOOoWNGDECI0aMAND6\nIpWXl0OpVCInJwerV6+GVqvFgAEDEBUVhYaGBmi1WsyYMQNyuZwaDsmnKs0c+GpuzV0Tm3Ngjr0h\nUjkhtAQTYhRQWJNT0Rc8EmmKxWkLVSawLR5sbSSJRvuTTz5Bnz59kJWVBblcjkOHDqFbt26CbWcD\n365gDzo6dSKNC4HfWaHVavHNN9/g3XffhU6nw4ABAwAA0dHRiI+PR3R0NNzd3XnLvwjE1NyyHZtN\nYsaW8BLCWYqVzDJnO131QLTPlkqE+ZTb2luZQCRspPnTsmXLkJ+fj9LSUvj7+yM2NhZffPEFPD09\nRTsvFykYwYoVK+Dt7S0l0joruBD4nRUuLi64dOkS3nnnHcycORMymQw1NTUoKChAfn4+PvroI6jV\naqqvRHx8PPr06QMAFpN0AEwcotilqpYSZXxpCSZnaW9emFkKS+4L4VAJLDUoN6fTtqfuFjAdn+Pl\n5YWmpiasWLECKpUKR44cgZ+fH65cuYLi4mLOLRi57CLnzZuHnJwcGAwGJCQkwNnZGbNmzWrTFay6\nuhqxsbFQq9VwcnJCRkYGLly4QE0OfpDQqSPdb7/9FocOHcJnn30GANi5cycKCgqwefPm+2xZ+4De\nV0KpVJrtK0EKDEhnKgAmHKUYjstacQaf47DJv0ikTLTI9opuhUT9pDKNqR8mUTzh5+Vyuej9bo1G\n0/E5zs7OKCwsxNKlS/Hmm29i2rRpgs7HZReZnZ2Njz76CNnZ2SgoKMD8+fMdriOYPdCpI10xnMXM\nmTORlZUFf39/nDt3TgSr2g9c+0oEBQXB2dkZJ0+exJEjR+Dp6Qm9Xg+NRkMdh05N8L2vJEIkfLUt\niTJmzwMSIWq1WipBRIYt2iqlAsSZr0YWLmYUTzhswmmTZjhCZWtstpPmOt7e3tBqtfjHP/6B33//\nHd988w1VYCAEXHaRBw4cQHJyMgAgPj4edXV1uHnzJgICAgSf90FAp3a6YhD4M2bMwFtvvYXp06eL\nbV67g62vhFKpxLRp0+Dr64tnnnkGr732GoxGIx577DFER0dj6NChCAwMbNNXgkuSzp68MGDKf9Jl\nYMyEFxstYY1ntfd2nyQo5XI5PDw8qEY+bLI1vppnNtvPnj2LhQsXYurUqUhPT7c5muYiA2P7THl5\nueR0raBTO92YmBj8+eefKC0tRVBQEP773/9i165dvI4xfPhwlJaW2sfADgB3d3esX78ezz//PMWL\nNjc349SpU1AqldToFdJXIi4uDkOGDIFcLrfYhpE4FXvwwvQtM5tDpKslSFMfLq0YiSO2Z/tFS86c\ni2zNmuaZuRDpdDqsW7cO//vf/5CZmYm+ffuKch18ikaE/N2DjE7tdJ2dnfHRRx9hzJgx0Ov1FIHf\n3igrK8P06dNx69YtyGQyvP7665g3b16728EGUjFHIJPJ4ObmhmHDhmHYsGEATPtKHDt2DBs2bEBj\nYyMiIiKoJB3pK6FWqyknQBwC4VptfeHoqge+MjA2WoLOs5JhmwQKhUIw72wOQjqCWUsusk0+VqlU\nCAkJweXLl5GamoqxY8fip59+ElX/zGUXyfxMeXk5evbsKZoNjopOnUgTC6WlpRg3bpxgTre6uhrV\n1dUYPHgwNBoNoqOjsX///k6tomD2lbh06RLUajWqqqqQlpaGV155Bd7e3iZJOqEyKsC+MjCgNbtP\noltnZ2dR5F8E7aVMILzx5MmToVQqoVAo8MILL+DZZ5/FU089hS5duoh2Ti4yMHoiTalUIjU1VUqk\ncUCnjnQ7CgIDAxEYGAgA8PLyQr9+/VBZWdmpnS69r0RycjISEhLg6+uLhQsXory8HG+88QZUKhVC\nQ0OpaDgiIgJOTk68xs7bWwZG56o9PT3bOHNr8i9ryUVmMktM29lUFdevXwcALFiwAI8//jhOnjyJ\nr776CmFhYaI6XXO7SLoM7Nlnn0V2djb69OkDT09PfPHFF6Kd35EhRbqwPdJlHmvkyJE4f/68Q2kQ\n8/LyMGLECJMts8FgQElJCRUNs/WV8PPzMyksoCfpSNtIvi0puUBoJzNrc+notpPJIvZKxNFHtgPA\njh07sHPnTmRkZCA2Ntbmc6hUKrzyyiu4fv06Hn30UezZs4fVcXdmhU9HxAPvdCdPnoxjx46hpqYG\n/v7+eP/99zFjxgxBx9JoNEhISMC7776LCRMm8Prbe/fuYeTIkVRk8/zzz2PNmjWC7LhfYPaVKCgo\nQGVlJQIDAxETE4O4uDgMGjQIMpkMN27cQFBQEACYRJS2lAYTiE1VMJup01teKhQKk2hYDF6bGd1W\nV1dj/vz56NevH1atWsV53Lk1LFmyBD169MCSJUvwwQcfoLa2lhoWSccvv/wCLy8vTJ8+XXK6IuCB\nd7pioaWlBWPHjsUzzzyD1NRUQcdobGyEh4cHdDodnnjiCaxfvx5PPPGEyJa2L+h9JZRKJY4ePYqy\nsjL07dsXs2fPRnR0NB555BETpyZUg8vWVUvs7T5dIkf6SbA1Uxeieab3wSA0y7fffostW7ZQz4KY\n1xMREYFjx44hICAA1dXVSEhIwMWLF1k/K+Zu8EGHxOmKAKPRiFmzZiEyMlKwwwVAlWeSMl6xG5Pc\nD8hkMoSEhCAkJARyuRy7du3Cxo0bERYWhsLCQqxbtw4lJSXw9fWlouGYmBi4uLiYLQ1mS3aJ0bjH\nEsxxt8x5aWTh4KN5Zo7PcXFxQU1NDRYuXAh/f3/k5ubC29tb1OsBYFLIEBAQgJs3b4p+DgltIUW6\nIuDXX3/FiBEjqJaMALBmzRo8/fTTvI5jMBgQFRWFkpISzJ07F2vXrrWHufcNGo0GWq22zWJiNBpN\n+kqcOHGC6itBRiGFhYWZaHGBv6JK4rTEbnlJbBOqTLA0pofJawOgSoSzsrKwbt06rF69GomJiTZd\nT2JiIqqrq9v8fvXq1UhOTkZtbS31u27dukGlUrEeR4p0xYPkdDsg7t69izFjxiA9PR0JCQm8/16v\n1yMmJgbBwcH44YcfxDewHcClr4RGo0FtbS1VEGDrtGMm6I3RyTwxW8HU4JJ+GOvWrYOnpyeKi4vh\n6+uLLVu22H0SbkREBPLy8hAYGIiqqiqMGjVKohfaARK90AHh6+uL5557DkVFRYKcbkZGBiIjI1Ff\nXy++ce0Etr4S9fX1KCoqwm+//YYVK1agrKwMY8eORVRUFOLi4jBgwACq4syWsULWKuJsATm/Vqul\nelXIZDJ0794dubm5qKioQEVFBf744w9kZmaaFLaIjfHjxyMzMxNLly5FZmYm7+SvBGEQr/5Rgk24\nc+cO6urqAABNTU04fPgwhgwZwvs45eXlyM7OxuzZs9uUaHZmyGQy+Pj4YPTo0bh16xZCQ0Nx/vx5\nrFmzBsHBwdi7dy9efPFFTJgwAcuXL0dOTg5qa2sph9nc3Iz6+nrU19ejsbGRaq7OvEctLS3UYuXt\n7W2X0ecajQZyuRyenp5obm7G4sWLcfr0aezZsweXL1+GSqXC1q1bERoayvscKpUKiYmJCAsLQ1JS\nEvVM0VFWVoZRo0bhwIEDWLVqFfz9/XHkyBEsW7YMAFBZWYnnnnuO+vzkyZPx+OOP4/LlywgJCZH0\nuDZCohc6CM6dO4fk5GQYDAYYDAa89tprWLx4Me/jvPTSS3j77behVquxfv36TksvWIJGo6EiRDqY\nfSWUSqVJX4nY2FhERUXB1dWVlWMlv/Pw8BBdd8vsM+zk5ESNa5o/fz6mTJkiCn3BRQbmiBWUnQmS\n03Ug/Pjjj8jJycHHH3+MvLw8bNiwQZDTffTRR+Hj4wO5XA6FQoHCwkI7WNs+MBqNqK6upnTDRUVF\nJn0l4uLioFKp0NLSgqioKOrvmNphMbhhImNrbm7G6tWrcfnyZXz66aei9ivgIwMjmDBhAt566y08\n+eSTotkhwTwkp+tAePvtt/HVV1/B2dkZ9+7dg1qtxosvvogdO3bwOk5oaCiKi4sdQrLGBtJX4siR\nI/jkk09w584dJCQkoG/fvoiLi0NsbCx8fHzaVKTxTdLRS5BJkcbp06exaNEizJgxA7Nnzxa1wxkA\ndO3alVIkGI1GdOvWzUShwISjVlB2ZEhO10Fx7NgxwfRCaGgoioqK0L17dztY1nGQnJwMd3d3rF27\nFgaDAYWFhcjPz0dBQYFJX4m4uDj069ePGlCp0+kAmLaNZCbp6ONz3NzcoNPpsH79eiiVSnz66afo\n3bu3YLvFkoHZUkEpQTgk9YIDQ+iWWCaT4amnnoJcLsecOXOQkpIismUdA59//rkJd5uUlISkpCQA\nrVHqlStXqAkcZ8+ehVwux+DBg036ShgMBiopRwohCFfs4uICd3d3/PHHH0hNTcXEiRNx8OBBmws3\nDh8+bPbfCK1AZGD+/v6sn2tpacGLL76IadOmSQ63nSFFuhLaoKqqCg899BBu376NxMREbN68GcOH\nD+d9nLq6OsyePRvnz5+HTCbD9u3bMXToUDtYbH+w9ZWoqKhAYGAglaTT6/W4efMmnn76adTV1SEm\nJgZ9+/bFnTt3sHjxYkyaNInqN2EvLFmyBN27d8fSpUuRnp6Ourq6Nok0o9GI5ORkdO/eHRs3brSr\nPRLaQnK6Eixi5cqV8PLyEjQyOzk5GSNHjsTMmTOh0+nQ0NAAX19fO1h5f0D6SuTl5eHDDz9ESUkJ\nRowYgZ49e+KRRx5Bbm4uIiMj4efnhxMnTqC4uBhXr16luobZAyqVCi+//DJu3Lhh0jmssrISKSkp\nyMrKEq2CUoIwSE5XggkaGxuh1+vh7e2NhoYGJCUlYfny5dS2myvu3r2LIUOG4OrVq3aytONg+fLl\nuHbtGjIyMuDp6YkzZ87gq6++QmJiIsaNG0d9jnQm4wsuLRgdoUvdgwLJ6UowwbVr1/DCCy8AaM3y\nT506FWlpabyPc/r0acyZMweRkZE4c+YMoqOjkZGRQTX1cSTo9XrRG+zQwbUFoyN2qXNESE5Xgl1Q\nVFSEYcOG4fjx44iNjUVqaip8fHzw/vvv8zrOpUuX8Oqrr1I/X716FatWreowM+jaA3y1t42NjRg5\nciQyMzMRGRnZjpZK4AKpDFiCXRAcHIzg4GBqwsGkSZNw8uRJ3scJDw/HqVOncOrUKRQXF8PDw4OK\nxB8UcG3BaDAYMHjwYAQEBGDUqFGSw+2gkCRjEuyCwMBAamJtWFgYcnNz0b9/f5uOmZubi969eyMk\nJEQkKzsOLGlv6bA0ncLJyQmnT5+mutTl5eUJapgkwb6QnK4Eu2Hz5s2YOnUqtFotevfubXOjlN27\nd2PKlCkiWdexIIb2lsDWLnUS7AuJXpBgNwwaNAgnTpzAmTNnsHfvXpvkYlqtFj/88ANeeuklQX+/\nZs0a9O/fHwMHDsSUKVPQ3Nws2Jb2BmnBCMBsC0axutRJsD8kpyuhUyAnJwfR0dHw8/Pj/belpaX4\n7LPPcPLkSZw7dw56vR67d++2g5X2wbJly3D48GGEhYWZbcFYWVmJ0aNHY/DgwYiPj8e4ceOkBjYd\nFBK9IKFTYNeuXZg8ebKgv/Xx8YFCoUBjYyPkcjkaGxtF7ewlBFzHnwOtdEFNTQ3Cw8NNemkEBQUh\nKysLAPDYY48JSlRKaH9Ika6EDo+Ghgbk5uZi4sSJgv6+W7duWLRoER5++GEEBQWhS5cueOqpp0S2\nkh/S09ORmJiIy5cv48knn2TV3RKQSSBizn6TcP8gOV0JHR6enp64c+eO4Im4JSUl2LRpE0pLS1FZ\nWQmNRoOvv/5asD0ZGRkYOHAgBgwYgIyMDEHHOHDgAJKTkwG0lkvv37+f9XOOOgnkQYbkdCU4PIqK\nivD444+je/fucHZ2xsSJE3H8+HFBx/r999/x+eefUwnCH3/8ESUlJbyPw1V7u2DBAqxbt070vrsS\n7h+kb1KCwyMiIgJKpRJNTU0wGo1UIxohuHjxIuLj4+Hm5ga5XI6RI0di7969rJ9NTEzEwIED2/x3\n4MABk8+Z097++OOP8Pf3x5AhQ6Qo14EgJdIkODwGDRqE6dOnIyYmBk5OToiKisLrr78u6FgDBgzA\nO++8A5VKBTc3N2RlZSEuLo71s7Zqb48fP44DBw4gOzubmgQyffp03pNAJHQsSL0XJEjgie3bt2PL\nli3w9PRE//794erqyrsvLZe+t3TYMglEQseCRC9IkMATM2fORFFREY4dO4YuXbogPDyc9zG4aG+Z\nkNQLjgEp0pUggSdu3boFf39/3LhxA2PGjEFBQQF8fHzut1kSOgkkTleCBJ6YNGkSampqoFAosGXL\nFsnhSuAFKdKVIEGChHaExOlKkCBBQjtCcroSJEiQ0I74f7vlau6CODNmAAAAAElFTkSuQmCC\n",
       "text": [
        "<matplotlib.figure.Figure at 0x7fbd40928550>"
       ]
      }
     ],
     "prompt_number": 92
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
     "prompt_number": 93
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
     "prompt_number": 94
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
     "prompt_number": 95
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
       "metadata": {},
       "output_type": "pyout",
       "png": "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\nMzE6NDUrMDE6MDDVFmMsAAAAAElFTkSuQmCC\n",
       "prompt_number": 96,
       "text": [
        "<IPython.core.display.Image at 0x7fbd408b0d50>"
       ]
      }
     ],
     "prompt_number": 96
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
     "prompt_number": 97
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 98,
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
     "prompt_number": 98
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
     "prompt_number": 99
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
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 100,
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
     "prompt_number": 100
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
     "prompt_number": 101
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "John.show3(topos=True,pattern=True)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 102
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