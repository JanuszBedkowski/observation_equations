from sympy import *

def rotX(a):
    return Matrix([[1, 0, 0, 0],
		[0, cos(a), -sin(a), 0],
		[0, sin(a), cos(a), 0],
		[0, 0, 0, 1]])
def rotY(a):
    return Matrix([[cos(a), 0, sin(a), 0],
		[0, 1, 0, 0],
		[-sin(a), 0, cos(a), 0],
		[0, 0, 0, 1]])
def rotZ(a):
    return Matrix([[cos(a), -sin(a), 0, 0],
		[sin(a), cos(a), 0, 0],
		[0, 0, 1, 0],
		[0, 0, 0, 1]])
def transXYZ(tx, ty, tz):
    return Matrix([[1, 0, 0, tx],
		[0, 1, 0, ty],
		[0, 0, 1, tz],
		[0, 0, 0, 1]])

def matrix44FromTaitBryan(tx, ty, tz, om, fi, ka):
    return transXYZ(tx, ty, tz) * rotX(om) * rotY(fi) * rotZ(ka)

def matrix44FromTaitBryanZYX(tx, ty, tz, om, fi, ka):
    return transXYZ(tx, ty, tz) * rotZ(ka) * rotY(fi) * rotX(om) 

def taitBryanFromMatrix44Case1(m):
    omfika = Matrix([0, 0, 0])
    omfika[0] = atan2(-m[1,2], m[2,2])#om
    omfika[1] = asin(m[0,2])          #fi
    omfika[2] = atan2(-m[0,1], m[0,0])#ka
    return omfika
   
def taitBryanFromMatrix44Case2(m):
    omfika = Matrix([0, 0, 0])
    omfika[0] = -atan2(m[1,0], m[1,1])#om
    omfika[1] = -pi / 2.0             #fi
    omfika[2] = 0                     #ka
    return omfika

def taitBryanFromMatrix44Case3(m):
    omfika = Matrix([0, 0, 0])
    omfika[0] = atan2(m[1,0], m[1,1]) #om
    omfika[1] = pi / 2.0              #fi
    omfika[2] = 0                     #ka
    return omfika
