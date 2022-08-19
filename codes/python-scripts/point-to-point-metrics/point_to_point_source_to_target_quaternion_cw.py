from sympy import *
import sys
sys.path.insert(1, '..')
from quaternion_R_utils import *

x_t, y_t, z_t = symbols('x_t y_t z_t')
x_s, y_s, z_s = symbols('x_s y_s z_s')
px, py, pz = symbols('px py pz')
q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

position_symbols = [px, py, pz]
quaternion_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + quaternion_symbols

point_t = Matrix([x_t, y_t, z_t]).vec()
point_s = Matrix([x_s, y_s, z_s, 1]).vec()

Rt_wc=matrix44FromQuaternion(px, py, pz, q0, q1, q2, q3)
R_cw=Rt_wc[:-1,:-1].transpose()
t_wc=Matrix([px, py, pz]).vec()
t_cw=-R_cw*t_wc
Rt_cw=Matrix.hstack(R_cw, t_cw)
Rt_cw=Matrix.vstack(Rt_cw, Matrix([[0,0,0,1]]))


transformed_point_s = (Rt_cw * point_s)[:-1,:]

delta=Matrix([0,0,0]).vec()-(transformed_point_s-point_t)
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_point_source_to_target_quaternion_cw_jacobian.h",'w') as f_cpp: 
    f_cpp.write("#ifndef _point_to_point_source_to_target_quaternion_cw_jacobian_h_")
    f_cpp.write("\n")
    f_cpp.write("#define _point_to_point_source_to_target_quaternion_cw_jacobian_h_")
    f_cpp.write("\n")  
    f_cpp.write("inline void point_to_point_source_to_target_quaternion_cw(double &delta_x, double &delta_y, double &delta_z, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_target_quaternion_cw_jacobian(Eigen::Matrix<double, 3, 7, Eigen::RowMajor> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_s, double y_s, double z_s)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (7):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("#endif")


