# -*- coding: latin1 -*-
#class Project(object)
#       """
#       Création d'une arborescence de projet
#       """
#       def __init__(self):
#       def

#-*- coding:Utf-8 -*-
"""
"""
import numpy as np
import os
import sys
import shutil
import pkgutil
import pdb

class PyLayers(object):
    """ Generic PyLayers Meta Class
    """

    def help(self,letter='az',typ='mt'):
        """ generic help

        Parameters
        ----------

        txt : string
            'mb' | 'mt'

        """

        members = self.__dict__.keys()
        lmeth = np.sort(dir(self))

        if typ=='mb':
            print np.sort(self.__dict__.keys())
        if typ=='mt':
            for s in lmeth:
                if s not in members:
                    if s[0]!='_':
                        if len(letter)>1:
                            if (s[0]>=letter[0])&(s[0]<letter[1]):
                                try:
                                    doc = eval('self.'+s+'.__doc__').split('\n')
                                    print s+': '+ doc[0]
                                except:
                                    pass
                        else:
                            if (s[0]==letter[0]):
                                try:
                                    doc = eval('self.'+s+'.__doc__').split('\n')
                                    print s+': '+ doc[0]
                                except:
                                    pass


currentdir = os.getcwd()


try:
    pylayersdir = os.environ['PYLAYERS']
except:
    pylayersdir = currentdir.split('pylayers')[0] + 'pylayers'


if pylayersdir[-1] == '/' or '\\':
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
    mesdir = os.path.join(basename ,'meas')

try:
    datadir = os.environ['DATADIR']
except:
    datadir = os.path.join(basename, 'meas')

try:
    os.path.isdir(os.path.join(basename ,'figures'))
except:
    os.mkdir(os.path.join(basename,'figures'))


# Dictionnary which associate PULSRAY environment variable with sub direrories
# of the project 
#
pstruc = {}
pstruc['DIRSIMUL'] ='ini'
pstruc['DIRSTRUC'] =os.path.join('struc','str')
pstruc['DIRWRL'] =os.path.join('struc','wrl')
pstruc['DIRINI'] =os.path.join('struc','ini')
pstruc['DIROSM'] =os.path.join('struc','osm')
pstruc['DIRSTRUC2'] = os.path.join('struc','str')
pstruc['DIRFUR'] = os.path.join('struc','furnitures')
pstruc['DIRIMAGE'] = os.path.join('struc','images')
pstruc['DIRPICKLE'] = os.path.join('struc','gpickle')
pstruc['DIRSLAB'] = 'ini'
pstruc['DIRSLAB2'] = 'ini'
pstruc['DIRMAT'] = 'ini'
pstruc['DIRMAT2'] = 'ini'
pstruc['DIRANT'] = 'ant'
pstruc['DIRTRA'] = 'output'
pstruc['DIRLCH'] = 'output'
pstruc['DIRTUD'] = 'output'
pstruc['DIRTx'] = os.path.join('output','Tx001')
pstruc['DIRGEOM'] = 'geom'
pstruc['DIRTRA'] = 'output'
pstruc['DIRCIR'] = 'output'
pstruc['DIRMES'] = 'meas'
pstruc['DIRNETSAVE'] = 'netsave'
pstruc['DIRSIG'] = os.path.join('output','sig')
pstruc['DIRR2D'] = os.path.join('output','r2d')
pstruc['DIRR3D'] = os.path.join('output','r3d')
pstruc['DIRCT'] = os.path.join('output','Ct')
pstruc['DIRH'] = os.path.join('output','H')
pstruc['DIRLNK'] = 'output'
pstruc['DIRBODY'] = 'body'
pstruc['DIRC3D'] = os.path.join('body','c3d')
pstruc['DIROOSM'] = os.path.join('gis','osm')
pstruc['DIRWEAR'] = os.path.join('body','wear')

# if basename directory does not exit it is created
try:
    os.chdir(basename)
except:
    print "Create directory " + basename
    os.mkdir(basename)

#
# write file project.conf
#
fd = open(os.path.join(basename,'project.conf'),'w')
fd.close()
#for nm in pstruc.keys():
for nm,nv in pstruc.items():
    dirname =  os.path.join(basename , pstruc[nm])
    spl = nv.split('/') # never again a variable called sp
    if len(spl)>1:
        if not os.path.isdir(os.path.join(basename ,spl[0])):
            os.mkdir(os.path.join(basename ,spl[0]))
            os.mkdir(os.path.join(basename,nv))
            print "create ",os.path.join(basename ,nv)
        else:
            if not os.path.isdir(os.path.join(basename ,nv)):
                os.mkdir(os.path.join(basename ,nv))
                print "create ",os.path.join(basename ,nv)
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
    if nm == 'DIROOSM':
        osmdir = dirname


    fd = open(os.path.join(basename,'project.conf'),'a')
    fd.write(nm+' '+dirname +'\n')
    fd.close()

#
# copy files from /data/ini in project directory
#

if basename != os.path.join(pylayersdir,'data'):
    if not 'win' in sys.platform:
        dirlist=['ini','struc','struc/furnitures'
        ,'struc/osm','struc/str','struc/wrl'
        ,'struc/images','struc/ini'
        ,'ant','output/Tx001','output'
        ,'geom','output/sig','output/r2d'
        ,'output/r3d','body','body/c3d','body/wear']
    else :
        dirlist=['ini','struc',os.path.join('struc','furnitures')
        ,os.path.join('struc','osm')
        ,os.path.join('struc','str'),os.path.join('struc','wrl')
        ,os.path.join('struc','images')
        ,os.path.join('struc','ini')
        ,'ant',os.path.join('output','Tx001'),'output'
        ,'geom',os.path.join('output','sig')
        ,os.path.join('output','r2d')
        ,os.path.join('output','r3d'),'body'
        ,os.path.join('body','c3d')
        ,os.path.join('body','wear')]
    for dl in dirlist:
        filelist = os.listdir(os.path.join(pylayersdir,'data', dl))
        for fi in filelist:
            if not os.path.isdir(os.path.join(basename,dl,fi)):
                if os.path.isfile(os.path.join(basename,dl,fi)): # file already exists
                    pass
                else:
                    shutil.copy(
                        os.path.join(pylayersdir,'data',dl,fi),
                        os.path.join(basename,dl,fi))


os.chdir(currentdir)
