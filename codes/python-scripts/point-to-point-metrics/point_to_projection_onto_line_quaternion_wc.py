from sympy import *
import sys
sys.path.insert(1, '..')
from quaternion_R_utils import *

x_src_l, y_src_l, z_src_l = symbols('x_src_l y_src_l z_src_l')
x_trg_g, y_trg_g, z_trg_g = symbols('x_trg_g y_trg_g z_trg_g')
x_trg_ln, y_trg_ln, z_trg_ln = symbols('x_trg_ln y_trg_ln z_trg_ln')
px, py, pz = symbols('px py pz')
q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

position_symbols = [px, py, pz]
quaternion_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + quaternion_symbols

point_source_local = Matrix([x_src_l, y_src_l, z_src_l,1]).vec()
RT_wc = matrix44FromQuaternion(px, py, pz, q0, q1, q2, q3)[:-1,:]
point_source_global = RT_wc * point_source_local
point_on_line_target_global = Matrix([x_trg_g, y_trg_g, z_trg_g])

a=point_source_global-point_on_line_target_global
b=Matrix([x_trg_ln, y_trg_ln, z_trg_ln]).vec()

p_proj = point_on_line_target_global + (a.dot(b)/b.dot(b))*b

delta = Matrix([0,0,0]).vec()-(point_source_global - p_proj)
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_projection_onto_line_quaternion_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_projection_onto_line_quaternion_wc(Eigen::Matrix<double, 3, 1> &delta, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double x_trg_ln, double y_trg_ln, double z_trg_ln)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta.coeffRef(2,0) = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_projection_onto_line_quaternion_wc_jacobian(Eigen::Matrix<double, 3, 7> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double x_trg_ln, double y_trg_ln, double z_trg_ln)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (7):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")

