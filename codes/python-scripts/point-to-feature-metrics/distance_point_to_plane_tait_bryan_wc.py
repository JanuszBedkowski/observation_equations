from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x, y, z = symbols('x y z')
a, b, c, d = symbols('a b c d')
px, py, pz = symbols('px py pz')
om, fi, ka = symbols('om fi ka')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

position_symbols = [px, py, pz]
orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + orientation_symbols

point_source_local = Matrix([x, y, z, 1]).vec()
RT_wc = matrix44FromTaitBryan(px, py, pz, om, fi, ka)[:-1,:]
#RT_wc = matrix44FromRodrigues(px, py, pz, sx, sy, sz)[:-1,:]
#RT_wc = matrix44FromQuaternion(px, py, pz, q0, q1, q2, q3)[:-1,:]
point_source_global = RT_wc * point_source_local

target_value = Matrix([0]).vec()
model_function = Matrix([a*point_source_global[0] + b * point_source_global[1] + c * point_source_global[2] + d]).vec()
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("distance_point_to_plane_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void delta_distance_point_to_plane_tait_bryan_wc(Eigen::Matrix<double, 1, 1> &delta, double px, double py, double pz, double om, double fi, double ka, double x, double y, double z, double a, double b, double c, double d)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void delta_distance_point_to_plane_tait_bryan_wc_jacobian(Eigen::Matrix<double, 1, 6> &j, double px, double py, double pz, double om, double fi, double ka, double x, double y, double z, double a, double b, double c, double d)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")

