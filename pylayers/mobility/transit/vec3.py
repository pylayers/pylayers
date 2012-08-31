#####################################################################
# vec3 - 3-dimensional vector
#
# Copyright (C) 2002, Matthias Baas (baas@ira.uka.de)
#
# You may distribute under the terms of the BSD license, as
# specified in the file license.txt.
####################################################################

import types, math, copy

# vec3
class vec3:
    """Three-dimensional vector.

    This class can be used to represent points, vectors, normals
    or even colors. The usual vector operations are available.
    """

    def __init__(self, *args):
        """Constructor.

        There are several possibilities how to initialize a vector:

        v = vec3()       -> v = <0,0,0>
        v = vec3(a)      -> v = <a,a,a>
        v = vec3(x,y)    -> v = <x,y,0>
        v = vec3(x,y,z)  -> v = <x,y,z>

        Note that specifying just one value sets all three components to
        that value (except when that single value is a another vec3, then
        that vector is copied).

        Additionally you can wrap those values in a list or a tuple or
        specify them as a string:

        v = vec3([1,2,3]) -> v = <1,2,3>
        v = vec3("4,5")   -> v = <4,5,0>        
        """
        
        if len(args)==0:
            self.x, self.y, self.z = (0.0, 0.0, 0.0)

        elif len(args)==1:
            T = type(args[0])
            # scalar
            if T==types.FloatType or T==types.IntType or T==types.LongType:
                self.x, self.y, self.z = (args[0], args[0], args[0])
            # vec3  
            elif isinstance(args[0], vec3):
                self.x, self.y, self.z = args[0]
            # Tuple/List
            elif T==types.TupleType or T==types.ListType:
                if len(args[0])==0:
                    self.x = self.y = self.z = 0.0
                elif len(args[0])==1:
                    self.x = self.y = self.z = args[0][0]
                elif len(args[0])==2:
                    self.x, self.y = args[0]
                    self.z         = 0.0
                elif len(args[0])==3:
                    self.x, self.y, self.z = args[0]
                else:
                    raise TypeError, "vec3() takes at most 3 arguments"
            # String
            elif T==types.StringType:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                f=map(lambda x: float(x), s)
                dummy = vec3(f)
                self.x, self.y, self.z = dummy
            # error
            else:
                raise TypeError,"vec3() arg can't be converted to vec3"

        elif len(args)==2:
            self.x, self.y, self.z = (args[0], args[1], 0.0)
            
        elif len(args)==3:
            self.x, self.y, self.z = args

        else:
            raise TypeError, "vec3() takes at most 3 arguments"


    def __repr__(self):
        return 'vec3('+`self.x`+', '+`self.y`+', '+`self.z`+')'

    def __str__(self):
        fmt="%1.4f"
        return '('+fmt%self.x+', '+fmt%self.y+', '+fmt%self.z+')'


    def __eq__(self, other):
        """== operator

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> c=vec3(-0.3, 0.75, 0.5)
        >>> print a==b
        0
        >>> print b==c
        1
        >>> print a==None
        0
        """
        if isinstance(other, vec3):
            return self.x==other.x and self.y==other.y and self.z==other.z
        else:
            return 0

    def __ne__(self, other):
        """!= operator

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> c=vec3(-0.3, 0.75, 0.5)
        >>> print a!=b
        1
        >>> print b!=c
        0
        >>> print a!=None
        1
        """
        if isinstance(other, vec3):
            return self.x!=other.x or self.y!=other.y or self.z!=other.z
        else:
            return 1


    def __add__(self, other):
        """Vector addition.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> print a+b
        (0.7000, 1.2500, -1.3000)
        """
        if isinstance(other, vec3):
            return vec3(self.x+other.x, self.y+other.y, self.z+other.z)
        else:
            raise TypeError, "unsupported operand type for +"

    def __sub__(self, other):
        """Vector subtraction.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> print a-b
        (1.3000, -0.2500, -2.3000)
        """
        if isinstance(other, vec3):
            return vec3(self.x-other.x, self.y-other.y, self.z-other.z)
        else:
            raise TypeError, "unsupported operand type for -"

    def __mul__(self, other):
        """Multiplication with a scalar or dot product.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> print a*2.0
        (2.0000, 1.0000, -3.6000)
        >>> print 2.0*a
        (2.0000, 1.0000, -3.6000)
        >>> print a*b
        -0.825
        """

        T = type(other)
        # vec3*scalar
        if T==types.FloatType or T==types.IntType or T==types.LongType:
            return vec3(self.x*other, self.y*other, self.z*other)
        # vec3*vec3
        if isinstance(other, vec3):
            return self.x*other.x + self.y*other.y + self.z*other.z
        # unsupported
        else:
            # Try to delegate the operation to the other operand
            if getattr(other,"__rmul__",None)!=None:
                return other.__rmul__(self)
            else:
                raise TypeError, "unsupported operand type for *"

    __rmul__ = __mul__

    def __div__(self, other):
        """Division by scalar

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print a/2.0
        (0.5000, 0.2500, -0.9000)
        """
        T = type(other)
        # vec3/scalar
        if T==types.FloatType or T==types.IntType or T==types.LongType:
            return vec3(self.x/other, self.y/other, self.z/other)
        # unsupported
        else:
            raise TypeError, "unsupported operand type for /"

    def __mod__(self, other):
        """Modulo (component wise)

        >>> a=vec3(3.0, 2.5, -1.8)
        >>> print a%2.0
        (1.0000, 0.5000, 0.2000)
        """
        T = type(other)
        # vec3%scalar
        if T==types.FloatType or T==types.IntType or T==types.LongType:
            return vec3(self.x%other, self.y%other, self.z%other)
        # unsupported
        else:
            raise TypeError, "unsupported operand type for %"

    def __iadd__(self, other):
        """Inline vector addition.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> a+=b
        >>> print a
        (0.7000, 1.2500, -1.3000)
        """
        if isinstance(other, vec3):
            self.x+=other.x
            self.y+=other.y
            self.z+=other.z
            return self
        else:
            raise TypeError, "unsupported operand type for +="

    def __isub__(self, other):
        """Inline vector subtraction.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> a-=b
        >>> print a
        (1.3000, -0.2500, -2.3000)
        """
        if isinstance(other, vec3):
            self.x-=other.x
            self.y-=other.y
            self.z-=other.z
            return self
        else:
            raise TypeError, "unsupported operand type for -="

    def __imul__(self, other):
        """Inline multiplication (only with scalar)

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> a*=2.0
        >>> print a
        (2.0000, 1.0000, -3.6000)
        """
        T = type(other)
        # vec3*=scalar
        if T==types.FloatType or T==types.IntType or T==types.LongType:
            self.x*=other
            self.y*=other
            self.z*=other
            return self
        else:
            raise TypeError, "unsupported operand type for *="

    def __idiv__(self, other):
        """Inline division with scalar

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> a/=2.0
        >>> print a
        (0.5000, 0.2500, -0.9000)
        """
        T = type(other)
        # vec3/=scalar
        if T==types.FloatType or T==types.IntType or T==types.LongType:
            self.x/=other
            self.y/=other
            self.z/=other
            return self
        else:
            raise TypeError, "unsupported operand type for /="

    def __imod__(self, other):
        """Inline modulo

        >>> a=vec3(3.0, 2.5, -1.8)
        >>> a%=2.0
        >>> print a
        (1.0000, 0.5000, 0.2000)
        """
        T = type(other)
        # vec3%=scalar
        if T==types.FloatType or T==types.IntType or T==types.LongType:
            self.x%=other
            self.y%=other
            self.z%=other
            return self
        else:
            raise TypeError, "unsupported operand type for %="

    def __neg__(self):
        """Negation

        >>> a=vec3(3.0, 2.5, -1.8)
        >>> print -a
        (-3.0000, -2.5000, 1.8000)
        """
        return vec3(-self.x, -self.y, -self.z)

    def __pos__(self):
        """
        >>> a=vec3(3.0, 2.5, -1.8)
        >>> print +a
        (3.0000, 2.5000, -1.8000)
        """
        return vec3(+self.x, +self.y, +self.z)

    def __abs__(self):
        """Return the length of the vector.

        abs(v) is equivalent to v.length().

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print abs(a)
        2.11896201004
        """
        return math.sqrt(self*self)


    def __len__(self):
        """Length of the sequence (always 3)"""
        return 3

    def __getitem__(self, key):
        """Return a component by index (0-based)

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print a[0]
        1.0
        >>> print a[1]
        0.5
        >>> print a[2]
        -1.8
        """
        T=type(key)
        if T!=types.IntType and T!=types.LongType:
            raise TypeError, "index must be integer"

        if   key==0: return self.x
        elif key==1: return self.y
        elif key==2: return self.z
        else:
            raise IndexError,"index out of range"

    def __setitem__(self, key, value):
        """Set a component by index (0-based)

        >>> a=vec3()
        >>> a[0]=1.5; a[1]=0.7; a[2]=-0.3
        >>> print a
        (1.5000, 0.7000, -0.3000)
        """
        T=type(key)
        if T!=types.IntType and T!=types.LongType:
            raise TypeError, "index must be integer"

        if   key==0: self.x = value
        elif key==1: self.y = value
        elif key==2: self.z = value
        else:
            raise IndexError,"index out of range"

    def cross(self, other):
        """Cross product.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> c=a.cross(b)
        >>> print c
        (1.6000, 0.0400, 0.9000)
        """
        
        if isinstance(other, vec3):
            return vec3(self.y*other.z-self.z*other.y,
                        self.z*other.x-self.x*other.z,
                        self.x*other.y-self.y*other.x)
        else:
            raise TypeError, "unsupported operand type for cross()"
        

    def length(self):
        """Return the length of the vector.

        v.length() is equivalent to abs(v).

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print a.length()
        2.11896201004
        """

        return math.sqrt(self*self)

    def normalize(self):
        """Return normalized vector.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print a.normalize()
        (0.4719, 0.2360, -0.8495)
        """

        nlen = 1.0/math.sqrt(self*self)
        return vec3(self.x*nlen, self.y*nlen, self.z*nlen)

    def angle(self, other):
        """Return angle (in radians) between self and other.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> print a.angle(b)
        1.99306755584
        """
        
        if isinstance(other, vec3):
            return math.acos((self*other) / (abs(self)*abs(other)))
        else:
            raise TypeError, "unsupported operand type for angle()"

    def reflect(self, N):
        """Return the reflection vector.

        N is the surface normal which has to be of unit length.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print a.reflect(vec3(1,0,1))
        (2.6000, 0.5000, -0.2000)
        """

        return self - 2.0*(self*N)*N

    def refract(self, N, eta):
        """Return the transmitted vector.

        N is the surface normal which has to be of unit length.
        eta is the relative index of refraction. If the returned
        vector is zero then there is no transmitted light because
        of total internal reflection.
        
        >>> a=vec3(1.0, -1.5, 0.8)
        >>> print a.refract(vec3(0,1,0), 1.33)
        (1.3300, -1.7920, 1.0640)
        """

        dot = self*N
        k   = 1.0 - eta*eta*(1.0 - dot*dot)
        if k<0:
            return vec3(0.0,0.0,0.0)
        else:
            return eta*self - (eta*dot + math.sqrt(k))*N

    def ortho(self):
        """Returns an orthogonal vector.

        Returns a vector that is orthogonal to self (where
        self*self.ortho()==0).

        >>> a=vec3(1.0, -1.5, 0.8)
        >>> print round(a*a.ortho(),8)
        0.0
        """

        x=abs(self.x)
        y=abs(self.y)
        z=abs(self.z)
        # Is z the smallest element? Then use x and y
        if z<=x and z<=y:
            return vec3(-self.y, self.x, 0.0)
        # Is y smallest element? Then use x and z
        elif y<=x and y<=z:
            return vec3(-self.z, 0.0, self.x)
        # x is smallest
        else:
            return vec3(0.0, -self.z, self.y)
        
	
    def ang0(self):
	C    = vec3(1.0,0.0,0.0)
	try:
		return self.angle(C)
	except: 
		return 0.0

    def vquarter(self):
	a0 	= copy.copy(self.normalize())
	a1 	= copy.copy(self.normalize())
	a0.x	= math.cos(self.ang0()+math.pi/6.0)*self.length()
	a1.x	= math.cos(self.ang0()-math.pi/6.0)*self.length()
	a0.y	= math.sin(self.ang0()+math.pi/6.0)*self.length()
	a1.y	= math.sin(self.ang0()-math.pi/6.0)*self.length()
	return [a0,self,a1]

######################################################################
def _test():
    import doctest, vec3
    failed, total = doctest.testmod(vec3)
    print "%d/%d failed" % (failed, total)

if __name__=="__main__":

    _test()

#    a = vec3(1,2,3.03)
#    b = vec3("-2,0.5,1E10")

#    print a.angle(b)
