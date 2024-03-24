#https://docs.opencv.org/4.8.0/db/d58/group__calib3d__fisheye.html
from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *

fx,fy,cx,cy = symbols('fx fy cx cy');
px, py, pz = symbols('px py pz')
om, fi, ka = symbols('om fi ka')
tie_px, tie_py, tie_pz = symbols('tie_px tie_py tie_pz'); 
u_kp, v_kp = symbols('u_kp v_kp')
k1,k2,k3,k4 = symbols("k1 k2 k3 k4")
alpha = symbols("alpha")

position_symbols = [px, py, pz]
orientation_symbols = [om, fi, ka]
tie_point_symbols = [tie_px, tie_py, tie_pz]
intrinsics_symbols = [k1, k2, k3, k4]
all_symbols = position_symbols + orientation_symbols

RT_wc = matrix44FromTaitBryan(px, py, pz, om, fi, ka)
R_cw=RT_wc[:-1,:-1].transpose()
T_wc=Matrix([px, py, pz]).vec()
T_cw=-R_cw*T_wc

point_global=Matrix([tie_px, tie_py, tie_pz]).vec()
point_local = R_cw * point_global + T_cw;
x = point_local[0];
y = point_local[1];
z = point_local[2];

a = x/z;
b = y/z;
r = sqrt(a*a + b*b);
Theta = atan(r);

Theta_d = Theta * (1 + k1 * Theta * Theta + k2 * Theta * Theta * Theta * Theta + k3 * Theta * Theta * Theta * Theta * Theta * Theta + k4 * Theta * Theta * Theta * Theta * Theta * Theta * Theta * Theta);

x_prim = (Theta_d/r) * a;
y_prim = (Theta_d/r) * b;

u = fx * (x_prim + alpha * y_prim) + cx;
v = fy * y_prim + cy;

u_delta = u_kp - u;
v_delta = v_kp - v;

uv = Matrix([u, v]).vec()
obs_eq = Matrix([u_delta, v_delta]).vec()
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(obs_eq)
print(obs_eq_jacobian)

with open("fisheye_camera_calibRT_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("#ifndef _fisheye_camera_calibRT_tait_bryan_wc_jacobian_\n");
    f_cpp.write("#define _fisheye_camera_calibRT_tait_bryan_wc_jacobian_\n\n");
    f_cpp.write("inline void observation_equation_fisheye_camera_tait_bryan_wc(Eigen::Matrix<double, 2, 1> &delta, double fx, double fy, double cx, double cy, double px, double py, double pz, double om, double fi, double ka, double tie_px, double tie_py, double tie_pz, double u_kp, double v_kp, double k1, double k2, double k3, double k4, double alpha)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(obs_eq[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(obs_eq[1,0])))
    f_cpp.write("}\n\n")
    f_cpp.write("inline void projection_fisheye_camera_tait_bryan_wc(double &u, double &v, double fx, double fy, double cx, double cy, double px, double py, double pz, double om, double fi, double ka, double tie_px, double tie_py, double tie_pz, double k1, double k2, double k3, double k4, double alpha)\n")
    f_cpp.write("{")
    f_cpp.write("u = %s;\n"%(ccode(uv[0,0])))
    f_cpp.write("v = %s;\n"%(ccode(uv[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n\n")
    f_cpp.write("inline void observation_equation_fisheye_camera_tait_bryan_wc_jacobian(Eigen::Matrix<double, 2, 6> &j, double fx, double fy, double px, double py, double pz, double om, double fi, double ka, double tie_px, double tie_py, double tie_pz, double k1, double k2, double k3, double k4, double alpha)\n")
    f_cpp.write("{")
    for i in range (2):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_jacobian[i,j])))
    f_cpp.write("}\n\n")
    f_cpp.write("#endif\n");
  

