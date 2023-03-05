from sympy import *
import sys
sys.path.insert(1, '..')
from quaternion_R_utils import *

x_t, y_t, z_t = symbols('x_t y_t z_t')
x_s, y_s, z_s = symbols('x_s y_s z_s')
px_0, py_0, pz_0 = symbols('px_0 py_0 pz_0')
q0_0, q1_0, q2_0, q3_0 = symbols('q0_0 q1_0 q2_0 q3_0')
px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
t0, t1, t = symbols('t0 t1 t') 
all_symbols = [t]

point_t = Matrix([x_t, y_t, z_t]).vec()
point_s = Matrix([x_s, y_s, z_s, 1]).vec()

alpha = (t - t0) / (t1 - t0)
px = (1.0 - alpha) * px_0 + alpha * px_1
py = (1.0 - alpha) * py_0 + alpha * py_1
pz = (1.0 - alpha) * pz_0 + alpha * pz_1

q0 = Matrix([q0_0, q1_0, q2_0, q3_0])
q1 = Matrix([q0_1, q1_1, q2_1, q3_1])
q = quat_slerp(q0, q1, alpha).vec()

transformed_point_s = (matrix44FromQuaternion(px, py, pz, q[0], q[1], q[2], q[3]) * point_s)[:-1,:]

delta=Matrix([0,0,0]).vec()-(transformed_point_s-point_t)
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("slerp_point_to_point_source_to_target_quaternion_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void slerp_point_to_point_source_to_target_quaternion_wc(double &delta_x, double &delta_y, double &delta_z, double px_0, double py_0, double pz_0, double q0_0, double q1_0, double q2_0, double q3_0, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double t0, double t1, double t, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void slerp_point_to_point_source_to_target_quaternion_wc_jacobian(Eigen::Matrix<double, 3, 1> &j, double px_0, double py_0, double pz_0, double q0_0, double q1_0, double q2_0, double q3_0, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double t0, double t1, double t, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (1):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void slerp_point_to_point_source_to_target_quaternion_wc_interpolate(double &px, double &py, double &pz, double &q0, double &q1, double &q2, double &q3, double px_0, double py_0, double pz_0, double q0_0, double q1_0, double q2_0, double q3_0, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double t0, double t1, double t)\n")
    f_cpp.write("{")
    f_cpp.write("px = %s;\n"%(ccode(px)))
    f_cpp.write("py = %s;\n"%(ccode(py)))
    f_cpp.write("pz = %s;\n"%(ccode(pz)))
    f_cpp.write("q0 = %s;\n"%(ccode(q[0])))
    f_cpp.write("q1 = %s;\n"%(ccode(q[1])))
    f_cpp.write("q2 = %s;\n"%(ccode(q[2])))
    f_cpp.write("q3 = %s;\n"%(ccode(q[3])))
    f_cpp.write("}")
    f_cpp.write("\n")

