from sympy import *
from quaternion_R_utils import *

px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')
px_3, py_3, pz_3 = symbols('px_3 py_3 pz_3')
q0_3, q1_3, q2_3, q3_3 = symbols('q0_3 q1_3 q2_3 q3_3')

position_symbols_1 = [px_1, py_1, pz_1]
quaternion_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
position_symbols_2 = [px_2, py_2, pz_2]
quaternion_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
position_symbols_3 = [px_3, py_3, pz_3]
quaternion_symbols_3 = [q0_3, q1_3, q2_3, q3_3]
all_symbols = position_symbols_1 + quaternion_symbols_1 + position_symbols_2 + quaternion_symbols_2 + position_symbols_3 + quaternion_symbols_3

RT_wc_1 = matrix44FromQuaternion(px_1,py_1,pz_1,q0_1,q1_1,q2_1,q3_1)
RT_wc_2 = matrix44FromQuaternion(px_2,py_2,pz_2,q0_2,q1_2,q2_2,q3_2)
RT_wc_3 = matrix44FromQuaternion(px_3,py_3,pz_3,q0_3,q1_3,q2_3,q3_3)

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

quaternion12=quaternionFromMatrix44(relative_pose12)
quaternion23=quaternionFromMatrix44(relative_pose23)
t12 = relative_pose12[0 : 3, 3]
t23 = relative_pose23[0 : 3, 3]

rel_pose12 = t12.col_join(quaternion12)
rel_pose23 = t23.col_join(quaternion23)

measurement = Matrix([0, 0, 0, 0, 0, 0, 0])

delta = measurement-(rel_pose12-rel_pose23)
delta_jacobian = delta.jacobian(all_symbols)

with open("smoothness_quaternion_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void smoothness_obs_eq_quaternion_wc(Eigen::Matrix<double, 7, 1> &delta, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double px_3, double py_3, double pz_3, double q0_3, double q1_3, double q2_3, double q3_3)\n")
    f_cpp.write("{")
    for i in range (7):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void smoothness_obs_eq_quaternion_wc_jacobian(Eigen::Matrix<double, 7, 21, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double px_3, double py_3, double pz_3, double q0_3, double q1_3, double q2_3, double q3_3)\n")
    f_cpp.write("{")
    for i in range (7):
        for j in range (21):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    
