# -*- coding: latin1 -*-
import os
import shutil
import pkgutil
#class Project(object)
#       """
#       Création d'une arborescence de projet
#       """
#       def __init__(self):
#       def
currentdir = os.getcwd()
dir1 = pkgutil.get_loader('pylayers').filename
os.chdir(dir1)
os.chdir('..')
pylayersdir = os.getcwd()
print "pylayers is in : " + pylayersdir
try:
    pulsraydir = os.environ['PULSRAY']
    print "PULSRAY  : ", pulsraydir
except:
    raise EnvironmentError('Please set the PULSRAY environment variable')
try:
    figuredir = os.environ['FIGURDIR']
except:
    figuredir = os.environ['HOME'] + "/Pyproject/figures"
#print currentdir
print "FIGURDIR : ", figuredir
try:
    basename = os.environ['BASENAME']
except:
    basename = os.environ['HOME'] + "/Pyproject"

#
# Dictionnary which associate PULSRAY environment variable with sub direrories
# of the project 
#
pstruc = {}

pstruc['DIRSTRUC'] ='struc'
pstruc['DIRSTRUC2'] = 'struc'
pstruc['DIRSLAB'] = 'ini'
pstruc['DIRSLAB2'] = 'ini'
pstruc['DIRMAT'] = 'ini'
pstruc['DIRMAT2'] = 'ini'
pstruc['DIRANT'] = 'ant'
pstruc['DIRTRA'] = 'output'
pstruc['DIRLCH'] = 'output'
pstruc['DIRTUD'] = 'output'
pstruc['DIRGEOM'] = 'geom'
pstruc['DIRTRA'] = 'output'
pstruc['DIRCIR'] = 'output'
pstruc['DIRMES'] = 'meas'

print "BASENAME : ", basename

# if basename directory does not exit it is created 
try:
    os.chdir(basename)
except:
    print "Create directory " + basename
    os.mkdir(basename)

fd = open(basename+'/project.conf','w')
fd.close()
for nm in pstruc.keys():
    dirname =  basename + '/'+pstruc[nm] 
    try:
        os.chdir(dirname)
        os.chdir('..')
    except:
        print "create "+ dirname
        os.mkdir(dirname)
        os.chdir('..')

    if nm == 'DIRSTRUC':
        strdir = dirname
    if nm == 'DIRGEOM':
        geomdir = dirname
    if nm == 'DIRLCH':
        lchdir = dirname
    if nm == 'DIRTUD':
        tuddir = dirname
    if nm == 'DIRSLAB':
        slabdir = dirname
    if nm == 'DIRMA':
        matdir = dirname
    if nm == 'DIRTRA':
        tradir = dirname
    if nm == 'DIRMES':
        mesdir = dirname
    if nm == 'DIRCIR':
        mesdir = dirname

    fd = open(basename+'/project.conf','a')
    fd.write(nm+' '+dirname +'\n')
    fd.close()
#
# copy files from /data/ini in project directory 
#
filelist = os.listdir(pylayersdir+'/data/ini')
for fi in filelist:
    if os.path.isfile(basename+'/ini/'+fi):
        print fi + '  already exists'
    else:
        shutil.copy(pylayersdir+'/data/ini/'+fi,basename+'/ini/'+fi)

#
# copy files from /data/struc in project directory 
#
filelist = os.listdir(pylayersdir+'/data/struc')
for fi in filelist:
    if os.path.isfile(basename+'/struc/'+fi):
        print fi + '  already exists'
    else:
        shutil.copy(pylayersdir+'/data/struc/'+fi,basename+'/struc/'+fi)
#
# copy files from /data/ant in project directory 
#
filelist = os.listdir(pylayersdir+'/data/ant')
for fi in filelist:
    if os.path.isfile(basename+'/ant/'+fi):
        print fi + '  already exists'
    else:
        shutil.copy(pylayersdir+'/data/ant/'+fi,basename+'/ant/'+fi)
os.chdir(currentdir)
