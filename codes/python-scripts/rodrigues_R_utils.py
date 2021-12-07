from sympy import *

def matrix44FromRodrigues(px, py, pz, sx, sy, sz):
    norm = sqrt( (sx)*(sx) + (sy)*(sy) + (sz)*(sz) )
    ux=sx/norm
    uy=sy/norm
    uz=sz/norm
    theta = norm
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

def rodriguesFromMatrix44(m):
    rodrigues = Matrix([0, 0, 0])
    r_x = m[2,1] - m[1,2]
    r_y = m[0,2] - m[2,0]
    r_z = m[1,0] - m[0,1]
    s = sqrt((r_x*r_x + r_y*r_y + r_z*r_z)*0.25)
    c = (m[0, 0] + m[1, 1] + m[2, 2] - 1)*0.5
    theta = acos(c)
    vth = 1/(2*s)
    vth = vth * theta
    r_x = r_x * vth
    r_y = r_y * vth
    r_z = r_z * vth
    rodrigues[0] = r_x
    rodrigues[1] = r_y
    rodrigues[2] = r_z
    return rodrigues
