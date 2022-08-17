from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_src_l, y_src_l, z_src_l = symbols('x_src_l y_src_l z_src_l')
x_trg_g, y_trg_g, z_trg_g = symbols('x_trg_g y_trg_g z_trg_g')
a, b, c = symbols('a b c')
tx_cw, ty_cw, tz_cw = symbols('tx_cw ty_cw tz_cw')
om_cw, fi_cw, ka_cw = symbols('om_cw fi_cw ka_cw')
#sx_cw, sy_cw, sz_cw = symbols('sx sy sz')
#q0_cw, q1_cw, q2_cw, q3_cw = symbols('q0_cw q1_cw q2_cw q3_cw')

position_symbols = [tx_cw, ty_cw, tz_cw]
orientation_symbols = [om_cw, fi_cw, ka_cw]
#orientation_symbols = [sx_cw, sy_cw, sz_cw]
#orientation_symbols = [q0_cw, q1_cw, q2_cw, q3_cw]
all_symbols = position_symbols + orientation_symbols

point_source_local=Matrix([x_src_l, y_src_l, z_src_l, 1]).vec()

#Rt_wc = matrix44FromTaitBryan(tx_cw, ty_cw, tz_cw, om_cw, fi_cw, ka_cw)[:-1,:]
RT_cw = matrix44FromTaitBryan(tx_cw, ty_cw, tz_cw, om_cw, fi_cw, ka_cw)
R_wc=RT_cw[:-1,:-1].transpose()
T_cw=Matrix([tx_cw, ty_cw, tz_cw]).vec()
T_wc=-R_wc*T_cw
Rt_wc=Matrix.hstack(R_wc, T_wc)
#Rt_wc=Matrix.vstack(Rt_wc, Matrix([[0,0,0,1]]))

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

print(delta)
print(delta_jacobian)

with open("point_to_projection_onto_plane_tait_bryan_cw_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_projection_onto_plane_tait_bryan_cw(Eigen::Matrix<double, 3, 1> &delta, double tx_cw, double ty_cw, double tz_cw, double om_cw, double fi_cw, double ka_cw, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double a, double b, double c)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta.coeffRef(2,0) = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_projection_onto_plane_tait_bryan_cw_jacobian(Eigen::Matrix<double, 3, 6> &j, double tx_cw, double ty_cw, double tz_cw, double om_cw, double fi_cw, double ka_cw, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double a, double b, double c)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")

