# -*- coding: latin1 -*-
import os
import sys
import shutil
import pkgutil
#class Project(object)
#       """
#       Création d'une arborescence de projet 
#       """
#       def __init__(self):
#       def     
<<<<<<< HEAD


currentdir=os.getcwd()
pylayersdir=currentdir.split('pylayers')[0] + '/pylayers/'

try:
    basename = os.environ['BASENAME']
    print "BASENAME  : ", basename
except:
    raise EnvironmentError('Please position an environement variable $BASENAME where your all your pylayers project will be hosted')


=======
currentdir = os.getcwd()
>>>>>>> master_bernard/master
try:
    pylayersdir = os.environ['PYLAYERS']
    print pylayersdir
except:
<<<<<<< HEAD
    raise EnvironmentError('Please set the PULSRAY environment variable')

try:
    os.path.isdir(basename +'/figures')
except:
    os.mkdir(basename+'/figures')



=======
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
>>>>>>> master_bernard/master
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
dirlist=['ini','struc','ant','output']

for dl in dirlist:


    filelist = os.listdir(pylayersdir+'/data/' + dl)
    for fi in filelist:
        if os.path.isfile(basename+'/' + dl +'/' +fi):
            if os.path.getsize(pylayersdir+'/data/' + dl +'/'+fi) != os.path.getsize(basename + '/' +dl +'/' + fi):
                print ('Would you like to replace your ' +fi +' by the default configuration file ?')
                A=raw_input()
                if A == 'y':
                    shutil.move(basename+'/' + dl +'/'+fi,basename+'/' +dl +'/'+fi+'.old')
                    shutil.copy(pylayersdir+'/data/' + dl +'/'+fi,basename+'/' + dl +'/'+fi)
        else:
            shutil.copy(pylayersdir+'/data/' + dl + '/'+fi,basename+'/' + dl +'/'+fi)


os.chdir(currentdir)
