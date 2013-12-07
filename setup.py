#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Setup script for pylayers 
"""
import numpy

from setuptools import setup,find_packages

setup(name='pylayers' ,
      version='0.1',
      description='Python LocAlization mobilitY Environement aware Radio Simulator',
      author='UGUEN Bernard, AMIOT Nicolas, and LAARAIEDH Mohamed',
      author_email='bernard.uguen@univ-rennes1.fr, nicolas.amiot@univ-rennes1.fr, mohamed.laaraeidh@univ-rennes1.fr',
      url='https://github.com/pylayers/pylayers',
      include_dirs = [numpy.get_include()],
      install_requires=[
        'numpy>=1.6.1',
        'scipy>=0.10.1',
        'networkx>=1.7',
        'matplotlib>=1.1.0',
        'shapely>=1.2.14',
        'descartes>=1.0',
        'SimPy>=2.2',
        'PIL>=1.1.5',
        'bitstring>=3.0.2',
        'shapely>=1.2.14',
        'descartes>=1.0',
        #'pyinterval>=1.0',
                        ],
      packages=find_packages()
)


