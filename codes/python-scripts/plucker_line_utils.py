from sympy import *

def plucker_line_K(fx, fy, cx, cy):
    return Matrix([[fy,0,0],[0,fx,0],[-fy*cx,-fx*cy,fx*fy]])

def plucker_line_motion_matrix_cw(Rt_wc):
    a=Rt_wc[0,3]
    b=Rt_wc[1,3]
    c=Rt_wc[2,3]
    rt=Matrix([[Rt_wc[0,0], Rt_wc[1,0], Rt_wc[2,0]],
               [Rt_wc[0,1], Rt_wc[1,1], Rt_wc[2,1]],
               [Rt_wc[0,2], Rt_wc[1,2], Rt_wc[2,2]]])
    tx=Matrix([[0,-c,b],[c,0,-a],[-b,a,0]])
    rttx=rt*tx
    return Matrix([[rt[0,0], rt[0,1], rt[0,2], -rttx[0,0], -rttx[0,1], -rttx[0,2]],
                   [rt[1,0], rt[1,1], rt[1,2], -rttx[1,0], -rttx[1,1], -rttx[1,2]],
                   [rt[2,0], rt[2,1], rt[2,2], -rttx[2,0], -rttx[2,1], -rttx[2,2]],
                   [0,0,0, rt[0,0], rt[0,1], rt[0,2]],
                   [0,0,0, rt[1,0], rt[1,1], rt[1,2]],
                   [0,0,0, rt[2,0], rt[2,1], rt[2,2]]])

def plucker_line_motion_matrix_wc(Rt_wc):
    a=Rt_wc[0,3]
    b=Rt_wc[1,3]
    c=Rt_wc[2,3]
    r=Matrix([[Rt_wc[0,0], Rt_wc[0,1], Rt_wc[0,2]],
              [Rt_wc[1,0], Rt_wc[1,1], Rt_wc[1,2]],
              [Rt_wc[2,0], Rt_wc[2,1], Rt_wc[2,2]]])
    tx=Matrix([[0,-c,b],[c,0,-a],[-b,a,0]])
    txr=tx*r
    return Matrix([[r[0,0], r[0,1], r[0,2], txr[0,0], txr[0,1], txr[0,2]],
                   [r[1,0], r[1,1], r[1,2], txr[1,0], txr[1,1], txr[1,2]],
                   [r[2,0], r[2,1], r[2,2], txr[2,0], txr[2,1], txr[2,2]],
                   [0,0,0, r[0,0], r[0,1], r[0,2]],
                   [0,0,0, r[1,0], r[1,1], r[1,2]],
                   [0,0,0, r[2,0], r[2,1], r[2,2]]])
