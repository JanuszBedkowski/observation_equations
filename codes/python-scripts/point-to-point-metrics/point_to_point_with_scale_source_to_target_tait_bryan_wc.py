from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x, y, z = symbols('x y z')
x_target, y_target, z_target = symbols('x_target y_target z_target')
s = symbols('s')
tx, ty, tz = symbols('tx ty tz')
om, fi, ka = symbols('om fi ka')

position_symbols = [tx, ty, tz]
orientation_symbols = [om, fi, ka]

all_symbols = position_symbols + orientation_symbols + [s]

point = Matrix([x, y, z, 1]).vec()
point_target = Matrix([x_target, y_target, z_target]).vec()

m_scale = Matrix([[s, 0, 0, 0], [0, s, 0, 0], [0, 0, s, 0], [0,0,0,1]])

transformed_point = ( matrix44FromTaitBryan(tx, ty, tz, om, fi, ka) * m_scale * point)[:-1,:]

target_value = Matrix([0,0,0]).vec()
model_function = transformed_point - point_target
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_point_with_scale_source_to_target_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_with_scale_source_to_target_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, [[maybe_unused]] double tx, [[maybe_unused]] double ty, [[maybe_unused]] double tz, [[maybe_unused]] double om, [[maybe_unused]] double fi, [[maybe_unused]] double ka, [[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] double x_target, [[maybe_unused]] double y_target, [[maybe_unused]] double z_target, [[maybe_unused]] double s)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_with_scale_source_to_target_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 7, Eigen::RowMajor> &j, [[maybe_unused]] double tx, [[maybe_unused]] double ty, [[maybe_unused]] double tz, [[maybe_unused]] double om, [[maybe_unused]] double fi, [[maybe_unused]] double ka, [[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] double x_target, [[maybe_unused]] double y_target, [[maybe_unused]] double z_target, [[maybe_unused]] double s)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (7):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")



