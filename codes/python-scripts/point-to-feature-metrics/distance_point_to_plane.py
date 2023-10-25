from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x, y, z = symbols('x y z')
a, b, c, d = symbols('a b c d')

position_symbols = [x, y, z]
all_symbols = position_symbols

point = Matrix([x, y, z, 1]).vec()
point_source_global = point

target_value = Matrix([0]).vec()
model_function = Matrix([a*point_source_global[0] + b * point_source_global[1] + c * point_source_global[2] + d]).vec()
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("distance_point_to_plane_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void delta_distance_point_to_plane(Eigen::Matrix<double, 1, 1> &delta, double x, double y, double z, double a, double b, double c, double d)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void delta_distance_point_to_plane_jacobian(Eigen::Matrix<double, 1, 3> &j, double x, double y, double z, double a, double b, double c, double d)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (3):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")

