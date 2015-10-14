import pdb
import numpy as np
from xml.etree.ElementTree import Element
from xml.etree import ElementTree

class ProbeBox(Element):
	def __init__(self, Name='port_ut1', Type='wv', Weight='1'):
"""

type:     0 for voltage probing => 'vp'
          1 for current probing => 'cp'
          2 for E-field probing => 'Efp'
          3 for H-field probing => 'Hfp'
          10 for waveguide voltage mode matching => 'wv'
         11 for waveguide current mode matching => 'wc'

"""
	Element.__init__(self,'ProbeBox')
	self.attrib['Weight']=Weight
	if Type=='vp'
		self.attrib['Type']='0'
	if Type=='cp'
		self.attrib['Type']='1'
	if Type=='Efp'
		self.attrib['Type']='2'
	if Type=='Hfp'
		self.attrib['Type']='3'
	if Type=='wv'
		self.attrib['Type']='10'
	if Type=='wc'
		self.attrib['Type']='11'
    self.append(Attributes('Attributes',a))



class Attributes(Element):
	def __init__(self, name, a):
	Element.__init__(self, name, ModeFunctionX=str(a[0]), ModeFunctionY=str(a[1]), ModeFunctionZ=str(a[2]))
