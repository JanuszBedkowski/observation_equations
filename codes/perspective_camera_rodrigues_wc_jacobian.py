from sympy import *
from rodrigues_R_utils import *

fx,fy,cx,cy = symbols('fx fy cx cy');
px, py, pz = symbols('px py pz')
sx, sy, sz = symbols('sx sy sz')
tie_px, tie_py, tie_pz = symbols('tie_px tie_py tie_pz'); 
u_kp, v_kp = symbols('u_kp v_kp')

position_symbols = [px, py, pz]
rodrigues_symbols = [sx, sy, sz]
tie_point_symbols = [tie_px, tie_py, tie_pz]
all_symbols = position_symbols + rodrigues_symbols + tie_point_symbols

RT_wc = matrix44FromRodrigues(px, py, pz, sx, sy, sz)
R_cw=RT_wc[:-1,:-1].transpose()
T_wc=Matrix([px, py, pz]).vec()
T_cw=-R_cw*T_wc

point_global=Matrix([tie_px, tie_py, tie_pz]).vec()
point_local = R_cw * point_global + T_cw;
s = 1.0 / point_local[2];
u = s * fx * point_local[0] + cx;
v = s * fy * point_local[1] + cy;
u_delta = u_kp - u;
v_delta = v_kp - v;

obs_eq = Matrix([u_delta, v_delta]).vec()
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(obs_eq)
print(obs_eq_jacobian)

with open("perspective_camera_rodrigues_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void observation_equation_perspective_camera_rodrigues_wc(Eigen::Matrix<double, 2, 1> &delta, double fx, double fy, double cx, double cy, double px, double py, double pz, double sx, double sy, double sz, double tie_px, double tie_py, double tie_pz, double u_kp, double v_kp)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(obs_eq[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(obs_eq[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_perspective_camera_rodrigues_wc_jacobian(Eigen::Matrix<double, 2, 9, Eigen::RowMajor> &j, double fx, double fy, double cx, double cy, double px, double py, double pz, double sx, double sy, double sz, double tie_px, double tie_py, double tie_pz)\n")
    f_cpp.write("{")
    for i in range (2):
        for j in range (9):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_jacobian[i,j])))
    f_cpp.write("}")
    
  

