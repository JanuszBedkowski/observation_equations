from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_src_l, y_src_l, z_src_l = symbols('x_src_l y_src_l z_src_l')
x_trg_g, y_trg_g, z_trg_g = symbols('x_trg_g y_trg_g z_trg_g')
a, b, c = symbols('a b c')
tx, ty, tz = symbols('tx ty tz')
om, fi, ka = symbols('om fi ka')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

position_symbols = [tx, ty, tz]
orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + orientation_symbols

point_source_local=Matrix([x_src_l, y_src_l, z_src_l, 1]).vec()
Rt_wc = matrix44FromTaitBryan(tx, ty, tz, om, fi, ka)[:-1,:]
#Rt_wc = matrix44FromRodrigues(tx, ty, tz, sx, sy, sz)[:-1,:]
#Rt_wc = matrix44FromQuaternion(tx, ty, tz, q0, q1, q2, q3)[:-1,:]
point_source_global = Rt_wc * point_source_local
point_source_global_1=Matrix([point_source_global[0], point_source_global[1], point_source_global[2], 1]).vec()
point_on_plane_target_global = Matrix([x_trg_g, y_trg_g, z_trg_g])
v_pl=Matrix([a, b, c]).vec()
d = -v_pl.dot(point_on_plane_target_global)

p_proj = point_source_global - (point_source_global_1.dot([[a,b,c,d]])) * v_pl

target_value = Matrix([0,0,0]).vec()
model_function = point_source_global - p_proj
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

delta_x = Matrix([delta[0]]).vec()
delta_hessian_x=delta_x.jacobian(all_symbols).jacobian(all_symbols)
delta_y = Matrix([delta[1]]).vec()
delta_hessian_y=delta_y.jacobian(all_symbols).jacobian(all_symbols)
delta_z = Matrix([delta[2]]).vec()
delta_hessian_z=delta_z.jacobian(all_symbols).jacobian(all_symbols)
delta_hessian = delta_hessian_x + delta_hessian_y + delta_hessian_z;

#print(delta)
#print(delta_jacobian)
#print(delta_hessian)


with open("point_to_projection_onto_plane_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_projection_onto_plane_tait_bryan_wc(Eigen::Matrix<double, 3, 1> &delta, double tx, double ty, double tz, double om, double fi, double ka, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double a, double b, double c)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta.coeffRef(2,0) = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_projection_onto_plane_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 6> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double a, double b, double c)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void point_to_projection_onto_plane_tait_bryan_wc_hessian(Eigen::Matrix<double, 6, 6> &h, double tx, double ty, double tz, double om, double fi, double ka, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double a, double b, double c)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (6):
            f_cpp.write("h.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_hessian[i,j])))
    f_cpp.write("}")

