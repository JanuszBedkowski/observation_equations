from sympy import *
from sympy import acos

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

#from Michal Pelka "Automation of the multi-sensor system calibration for mobile robotic applications" PhD

def quat_len(q):
    return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3])

def quat_norm(q):
    d = quat_len(q)
    return [q[0]/d, q[1]/d, q[2]/d, q[3]/d]
    
def quat_conj(q):
    return [q[0], -q[1],-q[2],-q[3]]

def quat_inv(q):
    l = quat_len (q)
    return [q[0]/l, -q[1]/l,-q[2]/l,-q[3]/l]

def quat_mul(q_a, q_b):
    [s_a, x_a, y_a, z_a] = q_a
    [s_b, x_b, y_b, z_b] = q_b
    q = [0,0,0,0]
    q[0] = s_a*s_b - x_a*x_b - y_a*y_b - z_a*z_b
    q[1] = s_a*x_b + x_a*s_b + y_a*z_b - z_a*y_b
    q[2] = s_a*y_b - x_a*z_b + y_a*s_b + z_a*x_b
    q[3] = s_a*z_b + x_a*y_b - y_a*x_b + z_a*s_b
    return Matrix(q)

def quat_log_map(q):
    phi = acos(q[0])
    sin_theta = sin(phi)
    u = [q[1]/sin_theta,q[2]/sin_theta,q[3]/sin_theta]
    return Matrix([0,phi*u[0],phi*u[1],phi*u[2]])

def quat_exp_map(q):
    phi = quat_len(q)
    sin_phi = sin(phi)
    u = quat_norm(q)
    return Matrix([cos(phi), sin_phi*u[1],sin_phi*u[2],sin_phi*u[3]])

def quat_slerp (q0, q1, t):
    from_q0_to_q1 = quat_mul(q1, quat_conj(q0))
    tan_from_q0_to_q1 = quat_log_map(from_q0_to_q1)
    temp = quat_exp_map(t*tan_from_q0_to_q1)
    return quat_mul(temp, q0)
