# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import scipy as sp
from pylayers.util import geomutil as geo
from pylayers.util import pyutil as pyu
import re

_filename = 'P-0007-BEAULIEU_BAT11D-x-2-E1.ifc'
filename = pyu.getlong(_filename,os.path.join('struc','bat11'))

fd  = open(filename,'r')
lig = fd.readlines()
dic1 = {}
tentry = []
for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        entry = t[1].split('(')[0]
        dic1[Id] = li.split(entry)[1]
        if entry not in tentry:
            tentry.append(entry)

#
#
#
for k in dic1:
    lv  = dic1[k]
    while lv.find('#')!=-1:
        m = re.search('(#[0-9]*)',lv)
        cle = m.group(0)
        lv  = lv.replace(m.group(0),dic1[cle])
    dic1[k]=lv

dico = {}
for entry in tentry:
    cle = entry.replace('IFC','')
    dico[cle]={}
    for li in lig:
        li = li.replace(';\r\n','')
        if li.find('=')!=-1:
            t = li.split('=')
            v  = t[1]
            Id = t[0]
            if v.find(entry)!=-1:
                dico[cle][Id] = dic1[Id].replace('(','').replace(')','').split(',')

#
# Reshaping
#
for k in dico['CARTESIANPOINT']:
    p = map(float,dico['CARTESIANPOINT'][k])
    dico['CARTESIANPOINT'][k] = p

