from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_L, y_L, z_L = symbols('x_L y_L z_L')
x_s, y_s, z_s = symbols('x_s y_s z_s')
px, py, pz = symbols('px py pz')
#om, fi, ka = symbols('om fi ka')
#sx, sy, sz = symbols('sx sy sz')
q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

position_symbols = [px, py, pz]
#orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
orientation_symbols = [q0, q1, q2, q3]
landmark_symbols = [x_L, y_L, z_L]

all_symbols = position_symbols + orientation_symbols + landmark_symbols

point_Landmark = Matrix([x_L, y_L, z_L]).vec()
point_source = Matrix([x_s, y_s, z_s, 1]).vec()
#transformed_point_source = (matrix44FromTaitBryan(px, py, pz, om, fi, ka) * point_source)[:-1,:]
#transformed_point_source = (matrix44FromRodrigues(px, py, pz, sx, sy, sz) * point_source)[:-1,:]
transformed_point_source = (matrix44FromQuaternion(px, py, pz, q0, q1, q2, q3) * point_source)[:-1,:]

target_value = Matrix([0,0,0]).vec()
model_function = transformed_point_source-point_Landmark
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_point_source_to_landmark_quaternion_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_source_to_landmark_quaternion_wc(double &delta_x, double &delta_y, double &delta_z, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_s, double y_s, double z_s, double x_L, double y_L, double z_L)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_landmark_quaternion_wc_jacobian(Eigen::Matrix<double, 3, 10, Eigen::RowMajor> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_s, double y_s, double z_s)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (10):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")



