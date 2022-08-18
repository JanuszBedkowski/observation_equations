from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x, y, z = symbols('x y z')
a, b, c, d = symbols('a b c d')
tx, ty, tz = symbols('tx ty tz')
om, fi, ka = symbols('om fi ka')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

position_symbols = [tx, ty, tz]
orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + orientation_symbols

point_source_local = Matrix([x, y, z, 1]).vec()
#Rt_wc = matrix44FromTaitBryan(tx, ty, tz, om, fi, ka)[:-1,:]
#Rt_wc = matrix44FromRodrigues(tx, ty, tz, sx, sy, sz)[:-1,:]
#Rt_wc = matrix44FromQuaternion(tx, ty, tz, q0, q1, q2, q3)[:-1,:]
RT_cw = matrix44FromTaitBryan(tx, ty, tz, om, fi, ka)
R_wc=RT_cw[:-1,:-1].transpose()
T_cw=Matrix([tx, ty, tz]).vec()
T_wc=-R_wc*T_cw
Rt_wc=Matrix.hstack(R_wc, T_wc)


point_source_global = Rt_wc * point_source_local

target_value = Matrix([0]).vec()
model_function = Matrix([a*point_source_global[0] + b * point_source_global[1] + c * point_source_global[2] + d]).vec()
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("distance_point_to_plane_tait_bryan_cw_jacobian.h",'w') as f_cpp:
    f_cpp.write("#ifndef _distance_point_to_plane_tait_bryan_cw_jacobian_h_")
    f_cpp.write("\n")
    f_cpp.write("#define _distance_point_to_plane_tait_bryan_cw_jacobian_h_")
    f_cpp.write("\n")     
    f_cpp.write("inline void delta_distance_point_to_plane_tait_bryan_cw(Eigen::Matrix<double, 1, 1> &delta, double tx, double ty, double tz, double om, double fi, double ka, double x, double y, double z, double a, double b, double c, double d)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void delta_distance_point_to_plane_tait_bryan_cw_jacobian(Eigen::Matrix<double, 1, 6> &j, double tx, double ty, double tz, double om, double fi, double ka, double x, double y, double z, double a, double b, double c, double d)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("#endif")

