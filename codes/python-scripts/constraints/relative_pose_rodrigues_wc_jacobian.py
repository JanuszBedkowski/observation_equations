from sympy import *
import sys
sys.path.insert(1, '..')
from rodrigues_R_utils import *

px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
sx_2, sy_2, sz_2 = symbols('sx_2 sy_2 sz_2')
px_m, py_m, pz_m = symbols('px_m py_m pz_m')
sx_m, sy_m, sz_m = symbols('sx_m sy_m sz_m')

position_symbols_1 = [px_1, py_1, pz_1]
rodrigues_symbols_1 = [sx_1, sy_1, sz_1]
position_symbols_2 = [px_2, py_2, pz_2]
rodrigues_symbols_2 = [sx_2, sy_2, sz_2]
all_symbols = position_symbols_1 + rodrigues_symbols_1 + position_symbols_2 + rodrigues_symbols_2

RT_wc_1 = matrix44FromRodrigues(px_1, py_1, pz_1, sx_1, sy_1, sz_1)
RT_wc_2 = matrix44FromRodrigues(px_2, py_2, pz_2, sx_2, sy_2, sz_2)

R_cw_1=RT_wc_1[:-1,:-1].transpose()
T_wc_1=Matrix([px_1, py_1, pz_1]).vec()
T_cw_1=-R_cw_1*T_wc_1
RT_cw_1=Matrix.hstack(R_cw_1, T_cw_1)
RT_cw_1=Matrix.vstack(RT_cw_1, Matrix([[0,0,0,1]]))

relative_pose = RT_cw_1 * RT_wc_2

rodrigues = rodriguesFromMatrix44(relative_pose)
t = relative_pose[0 : 3, 3]

observation = t.col_join(rodrigues)
measurement = Matrix([px_m, py_m, pz_m, sx_m, sy_m, sz_m])
delta = (measurement-observation)
delta_jacobian = delta.jacobian(all_symbols)

with open("relative_pose_rodrigues_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void relative_pose_obs_eq_rodrigues_wc(Eigen::Matrix<double, 6, 1> &delta, double px_1, double py_1, double pz_1, double sx_1, double sy_1, double sz_1, double px_2, double py_2, double pz_2, double sx_2, double sy_2, double sz_2, double px_m, double py_m, double pz_m, double sx_m, double sy_m, double sz_m)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_obs_eq_rodrigues_wc_jacobian(Eigen::Matrix<double, 6, 12, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double sx_1, double sy_1, double sz_1, double px_2, double py_2, double pz_2, double sx_2, double sy_2, double sz_2)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void relative_pose_rodrigues_wc(Eigen::Matrix<double, 6, 1> &relative_pose, double px_1, double py_1, double pz_1, double sx_1, double sy_1, double sz_1, double px_2, double py_2, double pz_2, double sx_2, double sy_2, double sz_2)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("relative_pose.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(observation[i])))
    f_cpp.write("}")
    f_cpp.write("\n")

