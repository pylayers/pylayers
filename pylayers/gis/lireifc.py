# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import scipy as sp
from pylayers.util import geomutil as geo
from pylayers.util import pyutil as pyu

_filename = 'P-0007-BEAULIEU_BAT11D-x-2-E1.ifc'
filename = pyu.getlong(_filename,os.path.join('struc','bat11'))

fd  = open(filename,'r')
lig = fd.readlines()
dico = {}
dico['CARTESIANPOINT']={}   # OK
dico['VECTOR']={}           # OK
dico['DIRECTION']={}        # OK
dico['LINE']={}             # OK
dico['POLYLINE']={}         # OK
dico['POLYLOOP']={}         # OK
dico['AXIS2PLACEMENT2D']={}
dico['AXIS2PLACEMENT3D']={}
dico['LOCALPLACEMENT']={}
dico['OPENINGELEMENT']={}   #
dico['SPACE']={}
dico['QUANTITYAREA']={}
dico['QUANTITYLENGTH']={}
dico['RELSPACEBOUNDARY']={}
dico['STYLEDITEM']={}
dico['CURVESTYLE']={}
dico['COLOURRGB']={}
dico['ELEMENTQUANTITY']={}
tentry = []
for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        entry = t[1].split('(')[0]
        if entry not in tentry:
            tentry.append(entry)
        if v.find('IFCCARTESIANPOINT')!=-1:
            spt = v.split('POINT')[1]
            sspt = spt.replace('(','').replace(')','').split(',')
            tp = []
            for item in sspt:
                tp.append(float(item))
            pt   = np.array(tp)
            dico['CARTESIANPOINT'][Id] = pt 

for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        if v.find('IFCDIRECTION')!=-1:
            spt  = v.split('DIRECTION')[1]
            sspt = spt.replace('(','').replace(')','').split(',')
            print sspt
            if len(sspt)==3:
                vec  = np.array([float(sspt[0]),float(sspt[1]),float(sspt[2])])
            else:
                vec  = np.array([float(sspt[0]),float(sspt[1])])
            dico['DIRECTION'][Id] = vec

for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        if v.find('IFCVECTOR')!=-1:
            spt = v.split('VECTOR')[1]
            sspt = spt.replace('(','').replace(')','').split(',')
            vec  = dico['DIRECTION'][sspt[0]]*float(sspt[1])
            dico['VECTOR'][Id] = vec

for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        if v.find('IFCLINE')!=-1:
            spt = v.split('LINE')[1]
            sspt = spt.replace('(','').replace(')','').split(',')
            idcp = sspt[0]
            idv = sspt[1]
            dico['LINE'][Id] = (dico['CARTESIANPOINT'][idcp],dico['VECTOR'][idv])

for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        if v.find('IFCPOLYLINE')!=-1:
            spt = v.split('POLYLINE')[1]
            sspt = spt.replace('(','').replace(')','').split(',')
            idta = sspt[0]
            idhe = sspt[1]
            dico['POLYLINE'][Id] = (dico['CARTESIANPOINT'][idta],dico['CARTESIANPOINT'][idhe])

for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        if v.find('IFCPOLYLOOP')!=-1:
            spt = v.split('POLYLOOP')[1]
            sspt = spt.replace('(','').replace(')','').split(',')
            tp = []
            for item in sspt:
                tp.append(dico['CARTESIANPOINT'][item])
            dico['POLYLOOP'][Id] = tp

for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        if v.find('IFCAXIS2PLACEMENT2D')!=-1:
            spt = v.split('PLACEMENT2D')[1]
            sspt = spt.replace('(','').replace(')','').split(',')
            tp   = (dico['CARTESIANPOINT'][sspt[0]],dico['DIRECTION'][sspt[1]])
            dico['AXIS2PLACEMENT2D'][Id] = tp

for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        if v.find('IFCAXIS2PLACEMENT3D')!=-1:
            spt = v.split('PLACEMENT3D')[1]
            sspt = spt.replace('(','').replace(')','').split(',')
            if sspt[0][0]=='#':
                v1 = dico['CARTESIANPOINT'][sspt[0]]
            else:
                v1 = sspt[0]
            if sspt[1][0]=='#':
                v2 =dico['DIRECTION'][sspt[1]]
            else:
                v2= sspt[1]
            if sspt[2][0]=='#':
                v3 =dico['DIRECTION'][sspt[2]]
            else:
                v3 =sspt[2]
            tp   = (v1,v2,v3)
            dico['AXIS2PLACEMENT3D'][Id] = tp

for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        if v.find('IFCLOCALPLACEMENT')!=-1:
            spt = v.split('LOCALPLACEMENT')[1]
            sspt = spt.replace('(','').replace(')','').split(',')

for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        if v.find('IFCOPENINGELEMENT')!=-1:
            spt = v.split('OPENINGELEMENT')[1]
            sspt = spt.replace('(','').replace(')','').split(',')
            dp = {}
            dp['ID']=sspt[0].replace('\'','')
            dp['name']=sspt[2].replace('\'','')
            dico['OPENINGELEMENT'][Id] = dp 

for li in lig:
    li = li.replace(';\r\n','')
    if li.find('=')!=-1:
        t = li.split('=')
        Id = t[0]
        v  = t[1]
        if v.find('IFCQUANTITYAREA')!=-1:
            dico['QUANTITYAREA'][Id] = v.split('AREA')[1]
        if v.find('IFCQUANTITYLENGTH')!=-1:
            dico['QUANTITYLENGTH'][Id] = v.split('LENGTH')[1]
        if v.find('IFCRELSPACEBOUNDARY')!=-1:
            dico['RELSPACEBOUNDARY'][Id] = v.split('BOUNDARY')[1]
        if v.find('IFCSPACE')!=-1:
            dico['SPACE'][Id] = v.split('SPACE')[1]
        if v.find('IFCCOLOURRGB')!=-1:
            dico['COLOURRGB'][Id] = v.split('RGB')[1]
        if v.find('IFCCURVESTYLE')!=-1:
            dico['CURVESTYLE'][Id] = v.split('CURVESTYLE')[1]
        if v.find('IFCSTYLEDITEM')!=-1:
            dico['STYLEDITEM'][Id] = v.split('STYLEDITEM')[1]
        if v.find('IFCELEMENTQUANTITY')!=-1:
            dico['ELEMENTQUANTITY'][Id] = v.split('ELEMENTQUANTITY')[1]


