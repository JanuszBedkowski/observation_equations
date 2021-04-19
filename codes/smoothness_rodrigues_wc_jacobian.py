from sympy import *
from rodrigues_R_utils import *

px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
sx_2, sy_2, sz_2 = symbols('sx_2 sy_2 sz_2')
px_3, py_3, pz_3 = symbols('px_3 py_3 pz_3')
sx_3, sy_3, sz_3 = symbols('sx_3 sy_3 sz_3')

position_symbols_1 = [px_1, py_1, pz_1]
rodrigues_symbols_1 = [sx_1, sy_1, sz_1]
position_symbols_2 = [px_2, py_2, pz_2]
rodrigues_symbols_2 = [sx_2, sy_2, sz_2]
position_symbols_3 = [px_3, py_3, pz_3]
rodrigues_symbols_3 = [sx_3, sy_3, sz_3]
all_symbols = position_symbols_1 + rodrigues_symbols_1 + position_symbols_2 + rodrigues_symbols_2 + position_symbols_3 + rodrigues_symbols_3

RT_wc_1 = matrix44FromRodrigues(px_1, py_1, pz_1, sx_1, sy_1, sz_1)
RT_wc_2 = matrix44FromRodrigues(px_2, py_2, pz_2, sx_2, sy_2, sz_2)
RT_wc_3 = matrix44FromRodrigues(px_3, py_3, pz_3, sx_3, sy_3, sz_3)

R_cw_1=RT_wc_1[:-1,:-1].transpose()
T_wc_1=Matrix([px_1, py_1, pz_1]).vec()
T_cw_1=-R_cw_1*T_wc_1
RT_cw_1=Matrix.hstack(R_cw_1, T_cw_1)
RT_cw_1=Matrix.vstack(RT_cw_1, Matrix([[0,0,0,1]]))

R_cw_2=RT_wc_2[:-1,:-1].transpose()
T_wc_2=Matrix([px_2, py_2, pz_2]).vec()
T_cw_2=-R_cw_2*T_wc_2
RT_cw_2=Matrix.hstack(R_cw_2, T_cw_2)
RT_cw_2=Matrix.vstack(RT_cw_2, Matrix([[0,0,0,1]]))

relative_pose12 = RT_cw_1 * RT_wc_2
relative_pose23 = RT_cw_2 * RT_wc_3

rodrigues12=rodriguesFromMatrix44(relative_pose12)
rodrigues23=rodriguesFromMatrix44(relative_pose23)
t12 = relative_pose12[0 : 3, 3]
t23 = relative_pose23[0 : 3, 3]

rel_pose12=t12.col_join(rodrigues12)
rel_pose23=t23.col_join(rodrigues23)

measurement = Matrix([0, 0, 0, 0, 0, 0])

delta = measurement-(rel_pose12-rel_pose23)
delta_jacobian = delta.jacobian(all_symbols)

with open("smoothness_rodrigues_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void smoothness_obs_eq_rodrigues_wc(Eigen::Matrix<double, 6, 1> &delta, double px_1, double py_1, double pz_1, double sx_1, double sy_1, double sz_1, double px_2, double py_2, double pz_2, double sx_2, double sy_2, double sz_2, double px_3, double py_3, double pz_3, double sx_3, double sy_3, double sz_3)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void smoothness_obs_eq_rodrigues_wc_jacobian(Eigen::Matrix<double, 6, 18, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double sx_1, double sy_1, double sz_1, double px_2, double py_2, double pz_2, double sx_2, double sy_2, double sz_2, double px_3, double py_3, double pz_3, double sx_3, double sy_3, double sz_3)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (18):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    
