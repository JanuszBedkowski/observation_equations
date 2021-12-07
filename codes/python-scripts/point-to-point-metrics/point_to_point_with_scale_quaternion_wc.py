from sympy import *
import sys
sys.path.insert(1, '..')
from quaternion_R_utils import *

x_1, y_1, z_1 = symbols('x_1 y_1 z_1')
s_1 = symbols('s_1')
px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
x_2, y_2, z_2 = symbols('x_2 y_2 z_2')
s_2 = symbols('s_2')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')

position_symbols_1 = [px_1, py_1, pz_1]
quaternion_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
position_symbols_2 = [px_2, py_2, pz_2]
quaternion_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
all_symbols = position_symbols_1 + quaternion_symbols_1 + [s_1] + position_symbols_2 + quaternion_symbols_2 + [s_2]

point_1 = Matrix([x_1, y_1, z_1, 1]).vec()
point_2 = Matrix([x_2, y_2, z_2, 1]).vec()
m_scale_1 = Matrix([[s_1, 0, 0, 0], [0, s_1, 0, 0], [0, 0, s_1, 0], [0,0,0,1]])
m_scale_2 = Matrix([[s_2, 0, 0, 0], [0, s_2, 0, 0], [0, 0, s_2, 0], [0,0,0,1]])

transformed_point_1 = ( matrix44FromQuaternion(px_1, py_1, pz_1, q0_1, q1_1, q2_1, q3_1) * m_scale_1 * point_1)[:-1,:]
transformed_point_2 = ( matrix44FromQuaternion(px_2, py_2, pz_2, q0_2, q1_2, q2_2, q3_2) * m_scale_2 * point_2)[:-1,:]

delta=transformed_point_1-transformed_point_2
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_point_with_scale_quaternion_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_with_scale_quaternion_wc(double &delta_x, double &delta_y, double &delta_z, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, double s_1, double s_2)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_with_scale_quaternion_wc_jacobian(Eigen::Matrix<double, 3, 16, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, double s_1, double s_2)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (16):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")



