from sympy import *
import sys
sys.path.insert(1, '..')
from rodrigues_R_utils import *

x_src_l, y_src_l, z_src_l = symbols('x_src_l y_src_l z_src_l')
x_trg_g, y_trg_g, z_trg_g = symbols('x_trg_g y_trg_g z_trg_g')
a, b, c = symbols('a b c')
px, py, pz = symbols('px py pz')
sx, sy, sz = symbols('sx sy sz')

position_symbols = [px, py, pz]
rodrigues_symbols = [sx, sy, sz]
all_symbols = position_symbols + rodrigues_symbols

point_source_local = Matrix([x_src_l, y_src_l, z_src_l, 1]).vec()
#RT_wc = matrix44FromRodrigues(px, py, pz, sx, sy, sz)[:-1,:]
RT_cw = matrix44FromRodrigues(px, py, pz, sx, sy, sz)
R_wc=RT_cw[:-1,:-1].transpose()
T_cw=Matrix([px, py, pz]).vec()
T_wc=-R_wc*T_cw
RT_wc=Matrix.hstack(R_wc, T_wc)

point_source_global = RT_wc * point_source_local
point_source_global_1=Matrix([point_source_global[0], point_source_global[1], point_source_global[2], 1]).vec()
point_on_plane_target_global = Matrix([x_trg_g, y_trg_g, z_trg_g])
v_pl=Matrix([a, b, c]).vec()
d = -v_pl.dot(point_on_plane_target_global)

p_proj = point_source_global - (point_source_global_1.dot([[a,b,c,d]])) * v_pl

delta = Matrix([0,0,0]).vec()-(point_source_global - p_proj)
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_projection_onto_plane_rodrigues_cw_jacobian.h",'w') as f_cpp:  
    f_cpp.write("#ifndef _point_to_projection_onto_plane_rodrigues_cw_jacobian_h_")
    f_cpp.write("\n")
    f_cpp.write("#define _point_to_projection_onto_plane_rodrigues_cw_jacobian_h_")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_projection_onto_plane_rodrigues_cw(Eigen::Matrix<double, 3, 1> &delta, double px, double py, double pz, double sx, double sy, double sz, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double a, double b, double c)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta.coeffRef(2,0) = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_projection_onto_plane_rodrigues_cw_jacobian(Eigen::Matrix<double, 3, 6> &j, double px, double py, double pz, double sx, double sy, double sz, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double a, double b, double c)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("#endif")

