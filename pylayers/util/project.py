# -*- coding: latin1 -*-
import os
#class Project(object)
#       """
#       Création d'une arborescence de projet
#       """
#       def __init__(self):
#       def
currentdir = os.getcwd()
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

print "BASENAME : ", basename
try:
    os.chdir(basename)
    for nm in pstruc.keys():
        dirname =  basename + '/'+pstruc[nm] 
        try:
            os.chdir(dirname)
            os.chdir('..')
        except:
            print "create "+ dirname
            os.mkdir(dirname)
            os.chdir('..')

        if nm = ['DIRSTRUC']:
            strdir = dirname 
    geomdir = basename + '/geom'
    try:
        os.chdir(geomdir)
        #print "geom dir already exist"
        os.chdir('..')
    except:
        os.mkdir(geomdir)
        print "create geom dir"
        os.chdir('..')

    lchdir = basename + '/launch'
    try:
        os.chdir(lchdir)
        #print "launch dir already exist"
        os.chdir('..')
    except:
        print "create launch dir"
        os.mkdir(lchdir)
        os.chdir('..')

    tuddir = basename + '/tud'
    try:
        os.chdir(tuddir)
        #print "tud dir already exist"
        os.chdir('..')
    except:
        print "create tud dir"
        os.mkdir(tuddir)
        os.chdir('..')

    slabdir = basename + '/slab'
    try:
        os.chdir(slabdir)
        #print "slab dir already exist"
        os.chdir('..')
    except:
        print "create slab dir"
        os.mkdir(slabdir)
        os.chdir('..')

    matdir = basename + '/mat'
    try:
        os.chdir(matdir)
        #print "mat dir already exist"
        os.chdir('..')
    except:
        print "create mat dir"
        os.mkdir(matdir)
        os.chdir('..')

    tradir = basename + '/trace'
    try:
        os.chdir(geomdir)
        #print "trace dir already exist"
        os.chdir('..')
    except:
        print "create trace dir"
        os.mkdir(tradir)
        os.chdir('..')

    geomdir = basename + 'geom'
    try:
        os.chdir(geomdir)
        #print "geom dir already exist"
        os.chdir('..')
    except:
        print "create geom dir"
        os.mkdir(geomdir)
        os.chdir('..')

    antdir = basename + '/ant'
    try:
        os.chdir(antdir)
        #print "ant dir already exist"
        os.chdir('..')
    except:
        print "create ant dir"
        os.mkdir(antdir)
        os.chdir('..')

    tuddir = basename + '/tud'
    try:
        os.chdir(tuddir)
        #print "tud dir already exist"
        os.chdir('..')
    except:
        print "create tud dir"
        os.mkdir(tuddir)
        os.chdir('..')

    mesdir = basename + '/measures'
    try:
        os.chdir(mesdir)
        #print "tud dir already exist"
        os.chdir('..')
    except:
        print "create mes dir"
        os.mkdir(mesdir)
        os.chdir('..')
    simuldir = basename + '/simul'
    try:
        os.chdir(simuldir)
        #print "tud dir already exist"
        os.chdir('..')
    except:
        print "create simul dir"
        os.mkdir(simuldir)
        os.chdir('..')
    wavedir = basename + '/wave'
    try:
        os.chdir(wavedir)
        #print "wavedir already exist"
        os.chdir('..')
    except:
        print "create wavedir"
        os.mkdir(wavedir)
        os.chdir('..')
    cirdir = basename + '/cir'
    try:
        os.chdir(cirdir)
        #print "wavedir already exist"
        os.chdir('..')
    except:
        print "create cirdir"
        os.mkdir(cirdir)
        os.chdir('..')


except:
    print "directory " + basename + " does'n exist"
    print "Create directory "
    os.mkdir(basename)

os.chdir(currentdir)
