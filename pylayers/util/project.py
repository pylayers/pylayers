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
try:
    pylayersdir = os.environ['PYLAYERS']
    print pylayersdir
except:
    raise NameError('PYLAYERS environemnt variable need to be defined')

try:
    puslarydir = os.environ['PULSRAY']
except:
    pulsraydir = pylayersdir

try:
    basename = os.environ['BASENAME']
except:
    basename = os.environ['HOME'] + "/PyLayers/Project"
    try:
        os.chdir(dirname)
    except:
        os.mkdir(dirname)

try:
    figuredir = os.environ['FIGURDIR']
except:
    figuredir = os.environ['HOME'] + "/Pyproject/figures"
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
        if os.path.getsize(pylayersdir+'/data/ini/'+fi) != os.path.getsize(basename+'/ini/'+fi):
            print ('Would you like to replace your ' +fi +' by the default configuration file ?')
            A=raw_input()
            if A == 'y':
                shutil.move(basename+'/ini/'+fi,basename+'/ini/'+fi+'.old')
                shutil.copy(pylayersdir+'/data/ini/'+fi,basename+'/ini/'+fi)
        else: 
            print fi + '  already exists'
    else:
        shutil.copy(pylayersdir+'/data/ini/'+fi,basename+'/ini/'+fi)


#
# copy files from /data/struc in project directory 
#
filelist = os.listdir(pylayersdir+'/data/struc')
for fi in filelist:
    if os.path.isfile(basename+'/struc/'+fi):
        if os.path.getsize(pylayersdir+'/data/struc/'+fi) != os.path.getsize(basename+'/struc/'+fi):
            print ('Would you like to replace your ' +fi +' by the default configuration file ?')
            A=raw_input()
            if A == 'y':
                shutil.move(basename+'/struc/'+fi,basename+'/struc/'+fi+'.old')
                shutil.copy(pylayersdir+'/data/struc/'+fi,basename+'/struc/'+fi)
        else: 
            print fi + '  already exists'
    else:
        shutil.copy(pylayersdir+'/data/struc/'+fi,basename+'/struc/'+fi)
#
# copy files from /data/ant in project directory 
#
filelist = os.listdir(pylayersdir+'/data/ant')
for fi in filelist:
    if os.path.isfile(basename+'/ant/'+fi):
        if os.path.getsize(pylayersdir+'/data/ant/'+fi) != os.path.getsize(basename+'/ant/'+fi):
            print ('Would you like to replace your ' +fi +' by the default configuration file ?')
            A=raw_input()
            if A == 'y':
                shutil.move(basename+'/ant/'+fi,basename+'/ant/'+fi+'.old')
                shutil.copy(pylayersdir+'/data/ant/'+fi,basename+'/ant/'+fi)
        else:
            print fi + '  already exists'
    else:
        shutil.copy(pylayersdir+'/data/ant/'+fi,basename+'/ant/'+fi)
#
# copy files from /data/output in project directory 
#
filelist = os.listdir(pylayersdir+'/data/output')
for fi in filelist:
    if os.path.isfile(basename+'/output/'+fi):
        if os.path.getsize(pylayersdir+'/data/output/'+fi) != os.path.getsize(basename+'/output/'+fi):
            print ('Would you like to replace your ' +fi +' by the default configuration file ?')
            A=raw_input()
            if A == 'y':
                shutil.move(basename+'/output/'+fi,basename+'/output/'+fi+'.old')
                shutil.copy(pylayersdir+'/data/output/'+fi,basename+'/output/'+fi)
        else :
            print fi + '  already exists'
    else:
        shutil.copy(pylayersdir+'/data/output/'+fi,basename+'/output/'+fi)
os.chdir(currentdir)
