from sympy import *
import sys
sys.path.insert(1, '..')
from rodrigues_R_utils import *

fx,fy,cx,cy = symbols('fx fy cx cy');
px, py, pz = symbols('px py pz')
sx, sy, sz = symbols('sx sy sz')
tie_px, tie_py, tie_pz = symbols('tie_px tie_py tie_pz'); 
u_kp, v_kp = symbols('u_kp v_kp')
k1,k2,k3,p1,p2 = symbols("k1 k2 k3 p1 p2")

position_symbols = [px, py, pz]
rodrigues_symbols = [sx, sy, sz]
tie_point_symbols = [tie_px, tie_py, tie_pz]
intrinsics_symbols = [k1,k2,k3,p1,p2]
all_symbols = intrinsics_symbols + position_symbols + rodrigues_symbols + tie_point_symbols

RT_cw = matrix44FromRodrigues(px, py, pz, sx, sy, sz)
R_cw=RT_cw[:-1,:-1]
T_cw=Matrix([px, py, pz]).vec()

point_global=Matrix([tie_px, tie_py, tie_pz]).vec()
point_local = R_cw * point_global + T_cw;
s = 1.0 / point_local[2];

u1 = s * point_local[0]
v1 = s * point_local[1]
r2 = u1 * u1 + v1 * v1;
u2 = u1 * (1 + k1*r2 + k2*r2*r2 + k3 * r2*r2*r2) + 2*p1*u1*v1 + p2*(r2 + 2 * u1 * u1)
v2 = v1 * (1 + k1*r2 + k2*r2*r2 + k3 * r2*r2*r2) + p1*(r2 + 2*v1*v1) + 2*(p2 * u1 * v1);

u3 = fx * u2 + cx
v3 = fy * v2 + cy
u_delta = u_kp - u3;
v_delta = v_kp - v3;

uv = Matrix([u3, v3]).vec()
obs_eq = Matrix([u_delta, v_delta]).vec()
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(obs_eq)
print(obs_eq_jacobian)

with open("perspective_camera_with_intrinsics_rodrigues_cw_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void observation_equation_perspective_camera_with_intrinsics_rodrigues_cw(Eigen::Matrix<double, 2, 1> &delta, double fx, double fy, double cx, double cy, double px, double py, double pz, double sx, double sy, double sz, double tie_px, double tie_py, double tie_pz, double u_kp, double v_kp, double k1, double k2, double k3, double p1, double p2)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(obs_eq[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(obs_eq[1,0])))
    f_cpp.write("}")
    f_cpp.write("inline void projection_perspective_camera_with_intrinsics_rodrigues_cw(double &u, double &v, double fx, double fy, double cx, double cy, double px, double py, double pz, double sx, double sy, double sz, double tie_px, double tie_py, double tie_pz, double k1, double k2, double k3, double p1, double p2)\n")
    f_cpp.write("{")
    f_cpp.write("u = %s;\n"%(ccode(uv[0,0])))
    f_cpp.write("v = %s;\n"%(ccode(uv[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_perspective_camera_with_intrinsics_rodrigues_cw_jacobian(Eigen::Matrix<double, 2, 14> &j, double fx, double fy, double px, double py, double pz, double sx, double sy, double sz, double tie_px, double tie_py, double tie_pz, double k1, double k2, double k3, double p1, double p2)\n")
    f_cpp.write("{")
    for i in range (2):
        for j in range (14):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_jacobian[i,j])))
    f_cpp.write("}")
    
  

