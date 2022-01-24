from sympy import *

def matrix44FromQuaternion(tx, ty, tz, q0, q1, q2, q3):
    q11 = q1*q1;
    q22 = q2*q2;
    q33 = q3*q3;
    q03 = q0*q3;
    q13 = q1*q3;
    q23 = q2*q3;
    q02 = q0*q2;
    q12 = q1*q2;
    q01 = q0*q1;
    mat0 = 1 - 2 * (q22 + q33);
    mat1 = 2.0*(q12+q03);
    mat2 = 2.0*(q13-q02);
    mat4 = 2.0*(q12-q03);
    mat5 = 1 - 2 * (q11 + q33);
    mat6 = 2.0*(q23+q01);
    mat8 = 2.0*(q13+q02);
    mat9 = 2.0*(q23-q01);
    mat10 = 1 - 2 * (q11 + q22);
    return Matrix([[mat0,mat4,mat8,tx],[mat1,mat5,mat9,ty],[mat2,mat6,mat10,tz],[0,0,0,1]])

def quaternionFromMatrix44(m):
    q = Matrix([1, 0, 0, 0])
    T = 1 + m[0,0] + m[1,1] + m[2,2]
    S = sqrt(T) * 2
    X = ( m[1,2] - m[2,1] ) / S
    Y = ( m[2,0] - m[0,2] ) / S
    Z = ( m[0,1] - m[1,0] ) / S
    W = 0.25 * S
    q[0] = W
    q[1] = -X
    q[2] = -Y
    q[3] = -Z
    sq = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3])
    q[0] = q[0]/sq
    q[1] = q[1]/sq
    q[2] = q[2]/sq
    q[3] = q[3]/sq
    return q
