#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Setup script for pylayers 
"""
import numpy
from setuptools import setup,find_packages
import os

setup(name='pylayers' ,
      version='0.5',
      description='Python LocAlization mobilitY Environement aware Radio Simulator',
      author='UGUEN Bernard, AMIOT Nicolas, LAARAIEDH Mohamed, MHEDHBI Meriem',
      author_email='bernard.uguen@univ-rennes1.fr, nicolas.amiot@univ-rennes1.fr, mohamed.laaraeidh@gmail.com',
      url='https://github.com/pylayers/pylayers',
      include_dirs = [numpy.get_include()],
      install_requires=[
        'numpy>=1.9.0',
        'scipy>=0.17.1',
        'networkx>=1.10',
        'matplotlib>=1.5.0',
        'SimPy==2.2',
        'PIL>=1.1.7',
        'bitstring>=3.1.3',
        'descartes>=1.0.1',
        'osmapi>=0.3.1',
        'GDAL>=2.0.0',
        'h5py>=2.5.0',
        'seaborn>=0.7.0',
        'basemap>=1.0.7',
        'protobuf==2.6.1',
        'mayavi>=4.4.0'




                        ],
      packages=find_packages()
)


# detect if project and source are already locate in .pylayers
# if not create both
# if so just update the source path


# detect user's home
home = os.path.expanduser('~')
# detect directory from which setup is launched
source = os.getcwd()
# project directory name. If no BASENAME env is already set, it is created @home
project = 'pylayers_project'

# check fresh install of pylayers
if not os.path.isfile(os.path.join(home,'.pylayers')):
    with open(os.path.join(home,'.pylayers'),'a') as f:
        f.write('source\n')
        f.write(source)
        basen_env = os.getenv('BASENAME')
        # test if env var BASENAME already set 
        #if not create a pyproject directory
        if basen_env == None:
            if not os.path.isdir(os.path.join(home,project)):
                os.mkdir(os.path.join(home,project))
        else:
            project=basen_env
        f.write('\nproject\n')
        f.write(os.path.join(home,project))

# if pylayers has already been installer, a .pylayers exists. 
# The idea here is to maintain the project path and to update the source path.
else: 
    with open(os.path.join(home,'.pylayers'),'r') as f:
        lines = f.readlines()
        # line corresponding to the source's path
    with open(os.path.join(home,'.pylayers'),'w') as f:
        for il,l in enumerate(lines):
            # line corresponding to the source's path
            if il != 1:
                f.write(l)
            else :
                f.write(source+"\n")

