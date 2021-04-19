from sympy import *
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
om_1, fi_1, ka_1 = symbols('om_1 fi_1 ka_1')
#sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
#q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
om_2, fi_2, ka_2 = symbols('om_2 fi_2 ka_2')
#sx_2, sy_2, sz_2 = symbols('sx_2 sy_2 sz_2')
#q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')
px_m, py_m, pz_m = symbols('px_m py_m pz_m')
om_m, fi_m, ka_m = symbols('om_m fi_m ka_m')
#sx_m, sy_m, sz_m = symbols('sx_m sy_m sz_m')
#q0_m, q1_m, q2_m, q3_m = symbols('q0_m q1_m q2_m q3_m')

position_symbols_1 = [px_1, py_1, pz_1]
orientation_symbols_1 = [om_1, fi_1, ka_1]
#orientation_symbols_1 = [sx_1, sy_1, sz_1]
#orientation_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
position_symbols_2 = [px_2, py_2, pz_2]
orientation_symbols_2 = [om_2, fi_2, ka_2]
#orientation_symbols_2 = [sx_2, sy_2, sz_2]
#orientation_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
all_symbols = position_symbols_1 + orientation_symbols_1 + position_symbols_2 + orientation_symbols_2

RT_wc_1 = matrix44FromTaitBryan(px_1, py_1, pz_1, om_1, fi_1, ka_1)
#RT_wc_1 = matrix44FromRodrigues(px_1, py_1, pz_1, sx_1, sy_1, sz_1)
#RT_wc_1 = matrix44FromQuaternion(px_1, py_1, pz_1, q0_1, q1_1, q2_1, q3_1)
RT_wc_2 = matrix44FromTaitBryan(px_2, py_2, pz_2, om_2, fi_2, ka_2)
#RT_wc_2 = matrix44FromRodrigues(px_2, py_2, pz_2, sx_2, sy_2, sz_2)
#RT_wc_2 = matrix44FromQuaternion(px_2, py_2, pz_2, q0_2, q1_2, q2_2, q3_2)

R_cw_1=RT_wc_1[:-1,:-1].transpose()
T_wc_1=Matrix([px_1, py_1, pz_1]).vec()
T_cw_1=-R_cw_1*T_wc_1
RT_cw_1=Matrix.hstack(R_cw_1, T_cw_1)
RT_cw_1=Matrix.vstack(RT_cw_1, Matrix([[0,0,0,1]]))

relative_pose = RT_cw_1 * RT_wc_2

parametrization = taitBryanFromMatrix44Case1(relative_pose)
#parametrization = rodriguesFromMatrix44(relative_pose)
#parametrization = quaternionFromMatrix44(relative_pose)
t = relative_pose[0 : 3, 3]

target_value=Matrix([px_m, py_m, pz_m, om_m, fi_m, ka_m])
#target_value=Matrix([px_m, py_m, pz_m, sx_m, sy_m, sz_m])
#target_value=Matrix([px_m, py_m, pz_m, q0_m, q1_m, q2_m, q3_m])
model_function=t.col_join(parametrization)
delta=target_value-model_function
delta_jacobian=delta.jacobian(all_symbols)

with open("relative_pose_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1(Eigen::Matrix<double, 6, 1> &delta, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double px_m, double py_m, double pz_m, double om_m, double fi_m, double ka_m)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1_jacobian(Eigen::Matrix<double, 6, 12, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void relative_pose_tait_bryan_wc_case1(Eigen::Matrix<double, 6, 1> &relative_pose, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("relative_pose.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(model_function[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
