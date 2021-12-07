from sympy import *

def matrix44FromRodrigues_utheta(px, py, pz, ux, uy, uz, theta):
    c = cos(theta)
    s = sin(theta)
    c1 = 1. - c
    rrt = Matrix([[ux*ux, ux*uy, ux*uz], [ux*uy, uy*uy, uy*uz], [ux*uz, uy*uz, uz*uz]])
    r_x = Matrix([[0, -uz,  uy], [uz,  0, -ux], [-uy, ux, 0]])
    c_eye =  Matrix([[c, 0, 0], [0, c ,0], [0, 0, c]])
    c1_rrt = c1 * rrt
    s_r_x = s * r_x
    R = c_eye + c1_rrt + s_r_x
    return Matrix([[R[0,0],R[0,1],R[0,2],px],[R[1,0],R[1,1],R[1,2],py],[R[2,0],R[2,1],R[2,2],pz],[0,0,0,1]])


