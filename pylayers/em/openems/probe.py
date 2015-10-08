import pdb
from xml.etree.ElementTree import Element
from xml.etree import ElementTree

class ProbeBox(Element):
    def __init__(self,Name,Type=0,Weight=-1):
        Element.__init__(self,'ProbeBox',
                         Name=Name,
                         Type=Type,
                         Wieght=Weight)

