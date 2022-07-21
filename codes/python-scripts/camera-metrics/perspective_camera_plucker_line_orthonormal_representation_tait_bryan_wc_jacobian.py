from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *
from plucker_line_utils import *

fx,fy,cx,cy = symbols('fx fy cx cy');
tx, ty, tz = symbols('tx ty tz')
om, fi, ka = symbols('om fi ka')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
#mx_w, my_w, mz_w, lx_w, ly_w, lz_w = symbols('mx_w my_w mz_w lx_w ly_w lz_w')
pl_om, pl_fi, pl_ka, pl_w = symbols('pl_om pl_fi pl_ka pl_w')
u_s,v_s = symbols('u_s v_s')
u_e,v_e = symbols('u_e v_e')

position_symbols = [tx, ty, tz]
orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
plucker_line_symbols = [pl_om, pl_fi, pl_ka, pl_w]
all_symbols = position_symbols + orientation_symbols + plucker_line_symbols

Rt_wc = matrix44FromTaitBryan(tx, ty, tz, om, fi, ka)
#Rt_wc = matrix44FromRodrigues(tx, ty, tz, sx, sy, sz)
#Rt_wc = matrix44FromQuaternion(tx, ty, tz, q0, q1, q2, q3)
R_pl = matrix44FromTaitBryan(0, 0, 0, pl_om, pl_fi, pl_ka) 

K=plucker_line_K(fx, fy, cx, cy)
mm_cw = plucker_line_motion_matrix_cw(Rt_wc)
#l_w = Matrix([[mx_w],[my_w],[mz_w],[lx_w],[ly_w],[lz_w]])
l_w = Matrix([[R_pl[0,0]*cos(pl_w)],[R_pl[1,0]*cos(pl_w)],[R_pl[2,0]*cos(pl_w)],[R_pl[0,1]*sin(pl_w)],[R_pl[1,1]*sin(pl_w)],[R_pl[2,1]*sin(pl_w)]])
l_c = mm_cw*l_w
m_c=Matrix([[l_c[0,0]],[l_c[1,0]],[l_c[2,0]]])
l=K*m_c

x_s=Matrix([[u_s,v_s,1]])
x_e=Matrix([[u_e,v_e,1]])

target_value = Matrix([0,0]).vec()
model_function = Matrix([(x_s*l)/sqrt(l[0,0]*l[0,0] + l[1,0]*l[1,0]),(x_e*l)/sqrt(l[0,0]*l[0,0] + l[1,0]*l[1,0])]).vec()
obs_eq = target_value - model_function
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(obs_eq)
print(obs_eq_jacobian)


with open("perspective_camera_plucker_line_orthonormal_representation_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void observation_equation_perspective_camera_plucker_line_orthonormal_representation_tait_bryan_wc(Eigen::Matrix<double, 2, 1> &delta, double fx, double fy, double cx, double cy, double tx, double ty, double tz, double om, double fi, double ka, double pl_om, double pl_fi, double pl_ka, double pl_w, double u_s, double v_s, double u_e, double v_e)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(obs_eq[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(obs_eq[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_perspective_camera_plucker_line_orthonormal_representation_tait_bryan_wc_jacobian(Eigen::Matrix<double, 2, 10, Eigen::RowMajor> &j, double fx, double fy, double cx, double cy, double tx, double ty, double tz, double om, double fi, double ka, double pl_om, double pl_fi, double pl_ka, double pl_w, double u_s, double v_s, double u_e, double v_e)\n")
    f_cpp.write("{")
    for i in range (2):
        for j in range (10):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_jacobian[i,j])))
    f_cpp.write("}")
