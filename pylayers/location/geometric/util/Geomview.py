# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of RGPA.

#Foobar is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#Foobar is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
#authors :
#Nicolas AMIOT          : nicolas.amiot@univ-rennes1.fr
#Bernard UGUEN          : bernard.uguen@univ-rennes1.fr
#Mohamed LAARAIEDH      : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
import numpy as np
import scipy as sp
import os
import pdb


def cloud(p, name="cloud", display=False, color='r', dice=2, R=0.5, access='new'):
    """
    cloud(p,filename,display,color) : display cloud of points p

    p            : cloud of points array(Npx3)
    display      : Boolean to switch display on/off
    color        : 'r','b','g','k'
    dice         : sphere sampling (2)
    R            : sphere radius   (0.5)
    access       : 'new' create a new file append mode neither

    """
    sh = np.shape(p)
    if len(sh) == 1:
        p = p.reshape((1, len(p)))
    Np = np.shape(p)[0]
    filename = './geom/' + name + '.list'
    if access == 'new':
        fd = open(filename, "w")
        fd.write("LIST\n")
    else:
        fd = open(filename, "a")
    sdice = " " + str(dice) + " " + str(dice) + " "
    radius = " " + str(R) + " "
    if color == 'r':
        col = " 1 0 0 "
    elif color == 'b':
        col = " 0 0 1 "
    elif color == 'm':
        col = " 1 0 1 "
    elif color == 'y':
        col = " 1 1 0 "
    elif color == 'c':
        col = " 0 1 1 "
    elif color == 'g':
        col = " 0 1 0 "
    elif color == 'k':
        col = " 0 0 0 "
    else:
        col = color
    for k in range(Np):
        try:
            c1 = str(p[k, 0]) + " " + str(p[k, 1]) + " " + str(p[k, 2])
        except:
            c1 = str(p[k, 0]) + " " + str(p[k, 1]) + " " + str(0)
        chaine = "{appearance {-edge patchdice" + sdice + "material { diffuse " + col + "}}{SPHERE" + radius + c1 + " }}\n"
        fd.write(chaine)

    fd.close()
    if display:
        chaine = "geomview  -nopanel  -b 1 1 1 " + filename + " 2>/dev/null &"
        os.system(chaine)

    return(filename)

if __name__ == "__main__":
    p = 10 * sp.randn(3000, 3)
    cloud(p, display=True)
