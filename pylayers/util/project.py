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
    pylayersdir = currentdir.split('pylayers')[0] + 'pylayers/'

if pylayersdir[-1] == '/':
    pylayersdir = pylayersdir[:-1]

if len(pylayersdir) == 1:
    raise EnvironmentError('Please verify that pylayers sources are into the "pylayers/" directory')

try:
    basename = os.environ['BASENAME']
except:
    raise EnvironmentError('Please position an environement variable $BASENAME where your pylayers project will be hosted')

try:
    mesdir = os.environ['MESDIR']
except:
    mesdir = basename + '/meas'

try:
    datadir = os.environ['DATADIR']
except:
    datadir = basename + '/meas'

try:
    os.path.isdir(basename +'/figures')
except:
    os.mkdir(basename+'/figures')


# Dictionnary which associate PULSRAY environment variable with sub direrories
# of the project 
#
pstruc = {}
pstruc['DIRSIMUL'] ='ini'
pstruc['DIRSTRUC'] ='struc/str'
pstruc['DIRWRL'] ='struc/wrl'
pstruc['DIRINI'] ='struc/ini'
pstruc['DIROSM'] ='struc/osm'
pstruc['DIRSTRUC2'] = 'struc/str'
pstruc['DIRFUR'] = 'struc/furnitures'
pstruc['DIRIMAGE'] = 'struc/images'
pstruc['DIRPICKLE'] = 'struc/gpickle'
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
#for nm in pstruc.keys():
for nm,nv in pstruc.items():
    dirname =  basename + '/'+pstruc[nm] 
    sp = nv.split('/')  
    if len(sp)>1:
        if not os.path.isdir(basename + '/'+sp[0]):     
            os.mkdir(basename + '/'+sp[0])
            os.mkdir(basename + '/'+nv)
            print "create ",basename + '/'+nv
        else:
            if not os.path.isdir(basename + '/'+nv):
                os.mkdir(basename + '/'+nv)
                print "create ",basename + '/'+nv
    else :
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
            print "create ",dirname
        


#    try:
#        os.chdir(dirname)
#        os.chdir('..')
#    except:
#        pdb.set_trace()
#        sp = nv.split('/')  
#        if len(sp)>1:
#            try:
#                os.chdir(basename + '/'+sp[0])
#                os.chdir('..')
#            except:
#                os.mkdir(basename + '/'+sp[0])
#                os.chdir(basename + '/'+sp[0])
#                os.mkdir(basename + '/'+sp[1])
#                os.chdir('..')
#        else:
#            print "create "+ dirname
#            os.mkdir(dirname)
#            os.chdir('..')


    if nm == 'DIRANT':
        antdir = dirname 
    if nm == 'DIRSTRUC':
        strdir = dirname
    if nm == 'DIRFUR':
        furdir = dirname
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

if basename<>pylayersdir+'/data':
    dirlist=['ini','struc','struc/furnitures','struc/osm','struc/str','struc/wrl','struc/images','struc/ini','ant','output','geom']
    for dl in dirlist:
        filelist = os.listdir(pylayersdir+'/data/' + dl)
        for fi in filelist:
            if not os.path.isdir(basename+'/'+dl+'/'+fi):
                if os.path.isfile(basename+'/' + dl +'/' +fi): # file already exists
                    pass
                else:
                    shutil.copy(pylayersdir+'/data/' + dl + '/'+fi,basename+'/' + dl +'/'+fi)
            

os.chdir(currentdir)
