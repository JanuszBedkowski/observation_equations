from sympy import *
from quaternion_R_utils import *

px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')
px_m, py_m, pz_m = symbols('px_m py_m pz_m')
q0_m, q1_m, q2_m, q3_m = symbols('q0_m q1_m q2_m q3_m')

position_symbols_1 = [px_1, py_1, pz_1]
quaternion_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
position_symbols_2 = [px_2, py_2, pz_2]
quaternion_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
all_symbols = position_symbols_1 + quaternion_symbols_1 + position_symbols_2 + quaternion_symbols_2

RT_wc_1 = matrix44FromQuaternion(px_1, py_1, pz_1, q0_1, q1_1, q2_1, q3_1)
RT_wc_2 = matrix44FromQuaternion(px_2, py_2, pz_2, q0_2, q1_2, q2_2, q3_2)

R_cw_1=RT_wc_1[:-1,:-1].transpose()
T_wc_1=Matrix([px_1, py_1, pz_1]).vec()
T_cw_1=-R_cw_1*T_wc_1
RT_cw_1=Matrix.hstack(R_cw_1, T_cw_1)
RT_cw_1=Matrix.vstack(RT_cw_1, Matrix([[0,0,0,1]]))

relative_pose = RT_cw_1 * RT_wc_2

quatenion = quaternionFromMatrix44(relative_pose)
t = relative_pose[0 : 3, 3]

observation = t.col_join(quatenion)
measurement = Matrix([px_m, py_m, pz_m, q0_m, q1_m, q2_m, q3_m])
delta = (measurement-observation)
delta_jacobian = delta.jacobian(all_symbols)

with open("relative_pose_quaternion_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void relative_pose_obs_eq_quaternion_wc(Eigen::Matrix<double, 7, 1> &delta, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double px_m, double py_m, double pz_m, double q0_m, double q1_m, double q2_m, double q3_m)\n")
    f_cpp.write("{")
    for i in range (7):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_obs_eq_quaternion_wc_jacobian(Eigen::Matrix<double, 7, 14, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2)\n")
    f_cpp.write("{")
    for i in range (7):
        for j in range (14):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void relative_pose_quaternion_wc(Eigen::Matrix<double, 7, 1> &relative_pose, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2)\n")
    f_cpp.write("{")
    for i in range (7):
        f_cpp.write("relative_pose.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(observation[i])))
    f_cpp.write("}")
    f_cpp.write("\n")

