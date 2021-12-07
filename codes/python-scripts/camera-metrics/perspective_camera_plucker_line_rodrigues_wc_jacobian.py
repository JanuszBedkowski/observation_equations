from sympy import *
import sys
sys.path.insert(1, '..')
from rodrigues_R_utils import *
from plucker_line_utils import *

fx,fy,cx,cy = symbols('fx fy cx cy');
px, py, pz = symbols('px py pz')
sx, sy, sz = symbols('sx sy sz')
mx_w, my_w, mz_w, lx_w, ly_w, lz_w = symbols('mx_w my_w mz_w lx_w ly_w lz_w')
u_s,v_s = symbols('u_s v_s')
u_e,v_e = symbols('u_e v_e')

position_symbols = [px, py, pz]
rodrigues_symbols = [sx, sy, sz]
plucker_line_symbols = [mx_w, my_w, mz_w]
all_symbols = position_symbols + rodrigues_symbols + plucker_line_symbols

RT_wc = matrix44FromRodrigues(px, py, pz, sx, sy, sz)
K=plucker_line_K(fx, fy, cx, cy)
mm_cw = plucker_line_motion_matrix_cw(RT_wc)
l_w = Matrix([[mx_w],[my_w],[mz_w],[lx_w],[ly_w],[lz_w]])
l_c = mm_cw*l_w
m_c=Matrix([[l_c[0,0]],[l_c[1,0]],[l_c[2,0]]])
l=K*m_c

x_s=Matrix([[u_s,v_s,1]])
x_e=Matrix([[u_e,v_e,1]])
xs_delta = Matrix([[0]]) - (x_s*l)/sqrt(l[0,0]*l[0,0] + l[1,0]*l[1,0])
xe_delta = Matrix([[0]]) - (x_e*l)/sqrt(l[0,0]*l[0,0] + l[1,0]*l[1,0])

obs_eq = Matrix([xs_delta, xe_delta]).vec()
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(obs_eq)
print(obs_eq_jacobian)

with open("perspective_camera_plucker_line_rodrigues_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void observation_equation_perspective_camera_plucker_line_rodrigues_wc(Eigen::Matrix<double, 2, 1> &delta, double fx, double fy, double cx, double cy, double px, double py, double pz, double sx, double sy, double sz, double mx_w, double my_w, double mz_w, double lx_w, double ly_w, double lz_w, double u_s, double v_s, double u_e, double v_e)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(obs_eq[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(obs_eq[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_perspective_camera_plucker_line_rodrigues_wc_jacobian(Eigen::Matrix<double, 2, 6, Eigen::RowMajor> &j, double fx, double fy, double cx, double cy, double px, double py, double pz, double sx, double sy, double sz, double mx_w, double my_w, double mz_w, double lx_w, double ly_w, double lz_w, double u_s, double v_s, double u_e, double v_e)\n")
    f_cpp.write("{")
    for i in range (2):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_jacobian[i,j])))
    f_cpp.write("}")
