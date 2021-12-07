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
def transXYZ(x, y, z):
    return Matrix([[1, 0, 0, x],
		[0, 1, 0, y],
		[0, 0, 1, z],
		[0, 0, 0, 1]])

def matrix44FromTaitBryan(px, py, pz, om, fi, ka):
    return transXYZ(px, py, pz) * rotX(om) * rotY(fi) * rotZ(ka)

def matrix44FromTaitBryanZYX(px, py, pz, om, fi, ka):
    return transXYZ(px, py, pz) * rotZ(ka) * rotY(fi) * rotX(om) 

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
