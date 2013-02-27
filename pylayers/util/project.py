# -*- coding: latin1 -*-
import os
import sys
import shutil
import pkgutil
import pdb
#class Project(object)
#       """
#       Création d'une arborescence de projet 
#       """
#       def __init__(self):
#       def     

currentdir = os.getcwd()

try:
    pylayersdir = os.environ['PYLAYERS']
except:
    pylayersdir = currentdir.split('pylayers')[0] + '/pylayers/'


if len(pylayersdir) == 1:
    raise EnvironmentError('Please verify that pylayers sources are into the "pylayers/" directory')

try:
    basename = os.environ['BASENAME']
except:
    raise EnvironmentError('Please position an environement variable $BASENAME where your all your pylayers project will be hosted')

try:
    mesdir = os.environ['MESDIR']
except:
    mesdir = basename + '/meas'

try:
    os.path.isdir(basename +'/figures')
except:
    os.mkdir(basename+'/figures')



# Dictionnary which associate PULSRAY environment variable with sub direrories
# of the project 
#
pstruc = {}
pstruc['DIRSIMUL'] ='ini'
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
pstruc['DIRNETSAVE'] = 'netsave'


# if basename directory does not exit it is created 
try:
    os.chdir(basename)
except:
    print "Create directory " + basename
    os.mkdir(basename)

#
# write file project.conf
#
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


    fd = open(basename+'/project.conf','a')
    fd.write(nm+' '+dirname +'\n')
    fd.close()
#
# copy files from /data/ini in project directory 
#


dirlist=['ini','struc','ant','output','geom']

for dl in dirlist:
    filelist = os.listdir(pylayersdir+'/data/' + dl)
    for fi in filelist:
        if os.path.isfile(basename+'/' + dl +'/' +fi):
            pass
        else:
            shutil.copy(pylayersdir+'/data/' + dl + '/'+fi,basename+'/' + dl +'/'+fi)


os.chdir(currentdir)
