from xml.etree import ElementTree
from xml.etree.ElementTree import Element
from xml.etree.ElementTree import SubElement


openEMS=Element('openEMS')
openEMS.append(Element('Metal',Name="gnd"))
openEMS.append(Element('Metal',Name="gnd2"))
FDTD=SubElement(openEMS,'FDTD',NumberOfTimeSteps="30000")
output_file = open( 'openEMS.xml', 'w' )
output_file.write( '<?xml version="1.0"?>' )
output_file.write( ElementTree.tostring( openEMS ) )
output_file.close()
