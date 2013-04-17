import numpy as np
import scipy.stats as sp
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import string
import struct

class Marker:
  def __init__(self, x, y, z):
    self.x=0.0
    self.y=0.0
    self.z=0.0

  def __repr__(self):
    return str("[x = "+str(self.x) +" y = " + str(self.y)+" z = "+ str(self.z)+"]")


class ParameterGroup:
  def __init__(self, nom, description, parameter):
    self.name = nom
    self.description = description
    self.parameter = parameter

  def __repr__(self):
    return self.name, " ", self.description, " ", self.parameter



class Parameter:
  def __init__(self, name, datatype, dim, data, description):
    self.name = name
    self.datatype = datatype
    self.dim = dim
    self.data = data
    self.description = description
 
    def __repr__(self):
      return self.name, " ", self.description, " ", self.dim


def getFloat(str):
    return struct.unpack('f',str[0:4])[0], str[4:]

def getNumber(str, length):
    sum = 0
        #sum = (sum << 8) + ord(str[i]) for big indian
        sum = sum + ord(str[i])*(2**(8*i))
    return sum, str[length:]

FullFileName = '07_01.c3d'


Markers=[];
VideoFrameRate=0;
AnalogSignals=[];
AnalogFrameRate=0;
Event=[];
ParameterGroups=[];
CameraInfo=[];
ResidualError=[];

print "*********************"
print "**** Opening File ***"
print "*********************"

print "FileName = ", FullFileName
fid=open(FullFileName,'rb'); # native format (PC-intel)
content = fid.read();
content_memory = content



NrecordFirstParameterblock, content = getNumber(content,1)     # Reading record number of parameter section


key, content = getNumber(content,1)

#if key~=80,
#h=errordlg(['File: ',FileName,' does not comply to the C3D format'],'application error');
#uiwait(h)
#fclose(fid)
#return
#end

if key!=80:
	print 'File: ',FullFileName,' does not comply to the C3D format'
	fid.close()

#fseek(fid,512*(NrecordFirstParameterblock-1)+3,'bof'); % jump to processortype - field
#proctype=fread(fid,1,'int8')-83;                       % proctype: 1(INTEL-PC); 2(DEC-VAX); 3(MIPS-SUN/SGI)
content = content[512*(NrecordFirstParameterblock-1)+1:]
proctype,content =getNumber(content,1)
proctype = proctype-83                      # proctype: 1(INTEL-PC); 2(DEC-VAX); 3(MIPS-SUN/SGI)

#print "*************************"
#print "**** Processor coding ***"
#print "*************************"

if proctype==1:
	print "Intel-PC"
elif proctype==2:
	print "DEC-VAX"
elif proctype==3:
	print "MIPS-SUN/SGI"
else:
	print "unknown processor type"


#if proctype==2,
#    fclose(fid);
#    fid=fopen(FullFileName,'r','d'); % DEC VAX D floating point and VAX ordering
#end
print "***********************"
print "**** Reading Header ***"
print "***********************"

# ###############################################
# ##                                           ##
# ##    read header                            ##
# ##                                           ##
# ###############################################

#%NrecordFirstParameterblock=fread(fid,1,'int8');     % Reading record number of parameter section
#%key1=fread(fid,1,'int8');                           % key = 80;

content = content_memory

#
#fseek(fid,2,'bof');
content = content[2:]


Nmarkers, content=getNumber(content, 2)
NanalogSamplesPerVideoFrame, content = getNumber(content, 2)

StartFrame,  content = getNumber(content, 2)
EndFrame,  content = getNumber(content, 2)

MaxInterpolationGap,  content = getNumber(content, 2)

Scale, content = getFloat(content)

NrecordDataBlock,  content = getNumber(content, 2)

NanalogFramesPerVideoFrame,  content = getNumber(content, 2)


if NanalogFramesPerVideoFrame > 0:
	NanalogChannels=NanalogSamplesPerVideoFrame/NanalogFramesPerVideoFrame
else:
	NanalogChannels=0

VideoFrameRate, content = getFloat(content)

AnalogFrameRate=VideoFrameRate*NanalogFramesPerVideoFrame

print "NanalogFramesPerVideoFrame= ", NanalogFramesPerVideoFrame
print "AnalogFrameRate= ", AnalogFrameRate
print "VideoFrameRate= ", VideoFrameRate
print "Scale= ", Scale
print "Nmarkers= ", Nmarkers
print "StartFrame= ", StartFrame
print "EndFrame= ", EndFrame

#print "***********************"
#print "**** Reading Events ..."
#print "***********************"

content = content_memory
content = content[298:] #bizarre .. ce devrait être 150 selon la doc

EventIndicator,  content = getNumber(content, 2)


EventTime=[]
EventValue=[]
EventName=[]
#
print "EventIndicator = ", EventIndicator
if EventIndicator==12345:
	Nevents,  content = getNumber(content, 2)
	print "Nevents= " ,Nevents
	content = content[2:]
	if Nevents>0:
		for i in range(Nevents):
			letime, content = getFloat(content)
			EventTime.append(letime)
		content = content_memory
		content = content[188*2:]
		for i in range(Nevents):
			lavalue, content = getNumber(content, 1)
			EventValue.append(lavalue)
		content = content_memory
		content = content[198*2:]
		for i in range(Nevents):
			lenom = content[0:4]
			content = content[4:]
			EventName.append(lenom)


print "***************************"
print "**** Reading Parameters ..."
print "***************************"

content = content_memory
content = content[512*(NrecordFirstParameterblock-1):] 

ParameterGroups = []
ParameterNumberIndex = []
##
dat1, content = getNumber(content, 1)
key2, content = getNumber(content, 1)
#
NparameterRecords, content = getNumber(content, 1)
print "NparameterRecords=", NparameterRecords
proctype,content =getNumber(content,1)
proctype = proctype-83                      # proctype: 1(INTEL-PC); 2(DEC-VAX); 3(MIPS-SUN/SGI)
#
for i in range(NparameterRecords):
	leparam = ParameterGroup(None, None, [])
	ParameterGroups.append(leparam)
	ParameterNumberIndex.append(0)
Ncharacters, content = getNumber(content, 1)
if Ncharacters>=128:
	Ncharacters = -(2**8)+(Ncharacters)
GroupNumber, content = getNumber(content, 1)
if GroupNumber>=128:
	GroupNumber = -(2**8)+(GroupNumber)
	print "GroupNumber = ", GroupNumber

while Ncharacters > 0:
	if GroupNumber<0:
		GroupNumber=abs(GroupNumber)
		GroupName = content[0:Ncharacters]
		content = content[Ncharacters:]
		print "Group Number = ", GroupNumber
		ParameterGroups[GroupNumber].name = GroupName
		print "ParameterGroupName =", GroupName
		offset, content = getNumber(content, 2)
		deschars, content = getNumber(content, 1)
		GroupDescription = content[0:deschars]
		content = content[deschars:]
		ParameterGroups[GroupNumber].description = GroupDescription
        
		ParameterNumberIndex[GroupNumber]=0
		ParameterNumber = 0 
		content = content[offset-3-deschars:]
	else:
		
		ParameterNumberIndex[GroupNumber]=ParameterNumberIndex[GroupNumber]+1
		ParameterNumber=ParameterNumberIndex[GroupNumber]
		print "ParameterNumber=", ParameterNumber
		ParameterGroups[GroupNumber].parameter.append(Parameter(None, None, [], [], None))
		ParameterName = content[0:Ncharacters]
		content = content[Ncharacters:]
		#print "ParameterName = ",ParameterName 
		if len(ParameterName)>0:
			ParameterGroups[GroupNumber].parameter[ParameterNumber-1].name=ParameterName
		offset, content = getNumber(content, 2)
		filepos = len(content_memory)-len(content)
		nextrec = filepos+offset-2
		
		type, content=getNumber(content, 1)
		if type>=128:
			type = -(2**8)+type
		ParameterGroups[GroupNumber].parameter[ParameterNumber-1].type=type
		
		dimnum, content=getNumber(content, 1)
		if dimnum == 0:
			datalength = abs(type)
		else:
			mult=1
			dimension=[]
			for j in range (dimnum):
				ladim, content = getNumber(content, 1)
			  	dimension.append(ladim)
			 	mult=mult*dimension[j]
			  	ParameterGroups[GroupNumber].parameter[ParameterNumber-1].dim.append(dimension[j])
			datalength = abs(type)*mult

#print "ParameterNumber = ", ParameterNumber, " Group Number = ", GroupNumber

	if type==-1:
		data = ""
		wordlength=dimension[0]
		if dimnum==2 and datalength>0:
			for j in range(dimension[1]):
				data=string.rstrip(content[0:wordlength])
	    			content = content[wordlength:]
	    			ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data.append(data)
		elif dimnum==1 and datalength>0:
	    		data=content[0:wordlength]
	    		ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data.append(data) # ???
		if string.rstrip(ParameterName) == "LABELS" and string.rstrip(GroupName) == "POINT":
	  		print "POINT = ", ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data
		elif string.rstrip(ParameterName) == "LABEL_PREFIXES" and string.rstrip(GroupName) == "SUBJECTS":
	  		print "SUBJECTS = ", ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data
		else:
	  #print  ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data
	  		pass


	elif type == 1:
		data = []
		Nparameters=datalength/abs(type)
		print "Nparameters=", Nparameters
		for i in range(Nparameters):
			ladata,content = getNumber(content, 1)
  			data.append(ladata)
		ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data=data
#print ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data

#print "type boolean"
	elif type == 2 and datalength>0:
		data = []
		Nparameters=datalength/abs(type)
		for i in range(Nparameters):
 			ladata,content = getNumber(content, 2)
  			data.append(ladata)
#ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data=data
		if dimnum>1:
  #???? print "arg je comprends pas"
  			ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data=data
  #???ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data=reshape(data,dimension)
		else:
  			ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data=data
#print ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data

#pass
#print "type integer"
	elif type == 4 and datalength>0:
		data = []
		Nparameters=datalength/abs(type)
		for i in range(Nparameters):
  			ladata,content = getFloat(content)  
  			data.append(ladata)
		if dimnum>1:
  			ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data=data
  #print "arg je comprends pas"
  #???ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data=reshape(data,dimension)
		else:
  			ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data=data
#print ParameterGroups[GroupNumber].parameter[ParameterNumber-1].data
	else:
#print "error"
		pass
	deschars, content= getNumber(content, 1)
	if deschars>0:
		description = content[0:deschars]
		content = content[deschars:]
		ParameterGroups[GroupNumber].parameter[ParameterNumber-1].description=description

	content = content_memory
	content = content[nextrec:] 

Ncharacters,content = getNumber(content, 1)
if Ncharacters>=128:
	Ncharacters = -(2**8)+(Ncharacters)
GroupNumber,content = getNumber(content, 1)
if GroupNumber>=128:
	GroupNumber = -(2**8)+(GroupNumber)
print "GroupNumber = ", GroupNumber

#
#
### ###############################################
### ##                                           ##
### ##    read data block                        ##
### ##                                           ##
### ###############################################
###  Get the coordinate and analog data
##
#
#content = content_memory
#content = content[(NrecordDataBlock-1)*512:] 
#
#NvideoFrames = EndFrame - StartFrame + 1
#print "NVideoFrames = ", NvideoFrames
#
#for i in range (NvideoFrames):
#Markers.append([])
#ResidualError.append([])
#CameraInfo.append([])
#for j in range (Nmarkers):
#Markers[i].append(Marker(0.0,0.0,0.0))
#ResidualError[i].append(0)
#CameraInfo[i].append(0)
#
##print Markers
##
##if Scale < 0
##    for i=1:NvideoFrames
##        for j=1:Nmarkers
##            Markers(i,j,1:3)=fread(fid,3,'float32')'; 
##            a=fix(fread(fid,1,'float32'));  
##            highbyte=fix(a/256);
##            lowbyte=a-highbyte*256; 
##            CameraInfo(i,j)=highbyte; 
##            ResidualError(i,j)=lowbyte*abs(Scale); 
##        end
##        waitbar(i/NvideoFrames)
##        for j=1:NanalogFramesPerVideoFrame,
##            AnalogSignals(j+NanalogFramesPerVideoFrame*(i-1),1:NanalogChannels)=...
##                fread(fid,NanalogChannels,'int16')';
##        end
##    end
#
#print "***************************"
#print "**** Reading DataBlock ...."
#print "***************************"
#ptr_read = 0
#if Scale < 0.0:
#for i in range (NvideoFrames):
#print "*",
#for j in range (Nmarkers):
##x, content = getFloat(content)
##y, content = getFloat(content)
##z, content = getFloat(content)
#x= struct.unpack('f',content[ptr_read:ptr_read+4])[0]
#ptr_read+=4
#y= struct.unpack('f',content[ptr_read:ptr_read+4])[0]
#ptr_read+=4
#z= struct.unpack('f',content[ptr_read:ptr_read+4])[0]
#ptr_read+=4
#Markers[i][j].x = x*Scale
#Markers[i][j].y = y*Scale
#Markers[i][j].z = z*Scale
##a, content  = getFloat(content)
#a= struct.unpack('f',content[ptr_read:ptr_read+4])[0]
#ptr_read+=4
#a = int(a)
#highbyte = int(a/256)
#lowbyte=a-highbyte*256
#CameraInfo[i][j] = highbyte
#ResidualError[i][j] = lowbyte*abs(Scale)
##if i< 2:
#  #print Markers[i][j]
#ptr_read+=NanalogFramesPerVideoFrame*NanalogChannels*2
##for j in range (NanalogFramesPerVideoFrame):
##  for k in range(NanalogChannels):
##    val, content = getNumber(content, 2)
##    AnalogSignals[j+NanalogFramesPerVideoFrame*(i)][k]=val #??? i-1
##else
##    for i=1:NvideoFrames
##        for j=1:Nmarkers
##            Markers(i,j,1:3)=fread(fid,3,'int16')'.*Scale;
##            ResidualError(i,j)=fread(fid,1,'int8');
##            CameraInfo(i,j)=fread(fid,1,'int8');
##        end
##        waitbar(i/NvideoFrames)
##        for j=1:NanalogFramesPerVideoFrame,
##            AnalogSignals(j+NanalogFramesPerVideoFrame*(i-1),1:NanalogChannels)=...
##                fread(fid,NanalogChannels,'int16')';
##        end
##    end
##end
#
#else:
#for i in range(NvideoFrames):
##print "**** Frame ", i, " *****"
#print "*",
#for j in range(Nmarkers):
##x, content = getNumber(content,2)
#x = ord(content[ptr_read]) + ord(content[ptr_read+1])*(2**8)
#ptr_read+=2
#if x > 32768:
#  x=-(2**16)+(x)
##y, content = getNumber(content,2)
#y = ord(content[ptr_read]) + ord(content[ptr_read+1])*(2**8)
#ptr_read+=2
#if y > 32768:
#  y=-(2**16)+(y)
##z, content = getNumber(content,2)
#z = ord(content[ptr_read]) + ord(content[ptr_read+1])*(2**8)
#ptr_read+=2
#if z > 32768:
#  z=-(2**16)+(z)
#Markers[i][j].x = x*Scale
#Markers[i][j].y = y*Scale
#Markers[i][j].z = z*Scale
##if i< 2:
##  print Markers[i][j]
#ResidualError[i][j], content = getNumber(content, 1)
#CameraInfo[i][j], content = getNumber(content, 1)
#ptr_read += NanalogFramesPerVideoFrame*NanalogChannels*2
##for j in range (NanalogFramesPerVideoFrame):
##  for k in range(NanalogChannels):
##    val, content = getNumber(content, 2)
#  #AnalogSignals(j+NanalogFramesPerVideoFrame*(i-1),1:NanalogChannels)=val
#print   
#
#empties=[]
#ipos=[]
#curvesX=[]
#curvesY=[]
#curvesZ=[]
#
#
#
#for i in range(Nmarkers):
#emptyname = 'em_'+str(i)
#em=GetOrCreateEmpty(emptyname)
#empties.append(em)
#eraseIPO(emptyname)
#
#try:
#lipo = Ipo.Get("Ipo"+str(emptyname)) 
#except:
#lipo = Ipo.New('Object', "Ipo"+str(emptyname))
#
##empties[i].clearIpo()
#empties[i].setLocation(0.0, 0.0, 0.0)
##lipo = Ipo.New('Object','ipo'+emptyname)
#ipos.append(lipo)
#empties[i].setIpo(ipos[i])
#curvesX.append(GetOrCreateCurve(ipos[i],'LocX'))
#curvesY.append(GetOrCreateCurve(ipos[i],'LocY'))
#curvesZ.append(GetOrCreateCurve(ipos[i],'LocZ'))
#
#print "**** Processing Frames ..."
#print "**************************"
#
#index=1
#while EndFrame > Blender.Get("curframe"):
#print "*",
#for i in range(Nmarkers):
#curvesX[i].addBezier((Blender.Get("curframe"), -Markers[index][i].x))
#curvesY[i].addBezier((Blender.Get("curframe"), -Markers[index][i].y))
#curvesZ[i].addBezier((Blender.Get("curframe"), -Markers[index][i].z))
#
#index = index+1
#num_frame = Blender.Get("curframe")
#Blender.Set("curframe", num_frame+1)
#print
#
#for i in range(Nmarkers):
#curvesX[i].Recalc()
#curvesY[i].Recalc()
#curvesZ[i].Recalc()
#
#context=scn.getRenderingContext() 
#context.startFrame(StartFrame)
#context.endFrame(EndFrame)
#context.framesPerSec(int(VideoFrameRate))
#Blender.Redraw()
#Blender.Set("curframe",StartFrame)
#print "Done."
#
#
##h = c3d.Header(open(file_name))
##r = c3d.Reader(open(file_name))
##frm=[]
#var = r.read_frames()
#for pnt, ang in var:
#	frm.append((pnt, ang))
#	x = pnt[:,0]
#	y = pnt[:,1]
#	z = pnt[:,2]
#		
#		
#		
#		 
#		
#	
#
