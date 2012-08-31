#!/usr/bin/python
# -*- coding: latin1 -*-
import numpy as np
from math import sqrt
#from scipy.linalg import *


class Vect(object):
        """
        Vector class
        ------------------------------------
        """

        def __init__(self,x=0,y=0,z=0):
                self.x = x
                self.y = y
                self.z = z

        def __add__(self,v):
                w = Vect()
                w.x = self.x + v.x
                w.y = self.y + v.y
                w.z = self.z + v.z
                return(w)

        def __sub__(self,v):
                w = Vect()
                w.x = self.x - v.x
                w.y = self.y - v.y
                w.z = self.z - v.z
                return(w)

        def __repr__(self):
                return ("("+str(self.x)+',' +
                           str(self.y)+',' +
                           str(self.z)+ ")")

        
        def __mul__(self,v):
                if isinstance(v,Vect):
                        ps = self.x*v.x+self.y*v.y+self.z*v.z
                        return(ps)
                else:
                        w = Vect()
                        w.x = v*self.x
                        w.y = v*self.y
                        w.z = v*self.z
                return(w)
        
        def __div__(self,a):
                w = Vect()
                try:
                        w.x = self.x/a
                        w.y = self.y/a
                        w.z = self.z/a
                except:
                        print("error : Vect division i")
                return(w)

        def __pow__(self,v):
                w   = Vect()
                w.x = self.y*v.z-v.y*self.z
                w.y = self.z*v.x-v.z*self.x
                w.z = self.x*v.y-v.x*self.y
                return(w)

        def show(self):
                print "(",self.x,',',self.y,',',self.z,")"

        def n2(self):
                norm = sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
                return(norm)

        def normalize(self):
                w  = Vect()
                n2 = self.n2()
                w  = self/n2
                return(w)

        def col(self):
                """
                Transfer a vector in a (3x1) numpy matrix
                """
                c = np.matrix(zeros([3,1]))
                c[0] = self.x
                c[1] = self.y
                c[2] = self.z
                return(c)
        
        def lgn(self):
                """
                Transfer a vector in a (1x3) numpy matrix
                """
                l = np.matrix(zeros([1,3]))
                l[0,0] = self.x
                l[0,1] = self.y
                l[0,2] = self.z
                return(l)

        def afficheMethodesEtAttributs(self):
             print
             print "Method List"
             for methode in dir(self):
                   print methode

class Point(Vect):

        def __repr__(self):
                return ("["+str(self.x)+',' +
                           str(self.y)+',' +
                           str(self.z)+ "]")

if (__name__ == "__main__"):
        v1 = Vect(1,0,.0)
        v2 = Vect(0,1,.0)
        v3 = Vect(1,1,1)

        v4 = v1*v2
        v5 = v1**v2

