from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

pose_from = symbols('pose_from')
pose_to = symbols('pose_to')
pose_m = symbols('pose_m')

position_symbols_1 = [pose_from]
position_symbols_2 = [pose_to]
all_symbols = position_symbols_1 + position_symbols_2

wc_1 = Matrix([[pose_from]])
wc_2 = Matrix([[pose_to]])

cw_1=wc_1.inv()
#T_wc_1=Matrix([px_1, py_1, pz_1]).vec()
#T_cw_1=-R_cw_1*T_wc_1
#RT_cw_1=Matrix.hstack(R_cw_1, T_cw_1)
#RT_cw_1=Matrix.vstack(RT_cw_1, Matrix([[0,0,0,1]]))

relative_pose = cw_1 * wc_2

#parametrization = taitBryanFromMatrix44Case1(relative_pose)
#parametrization = rodriguesFromMatrix44(relative_pose)
#parametrization = quaternionFromMatrix44(relative_pose)
#t = relative_pose[0 : 3, 3]

target_value=Matrix([pose_m])

#target_value=Matrix([px_m, py_m, pz_m, sx_m, sy_m, sz_m])
#target_value=Matrix([px_m, py_m, pz_m, q0_m, q1_m, q2_m, q3_m])
model_function=relative_pose
delta=target_value-model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)


with open("relative_pose_one_dimentional_example_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void relative_pose_obs_eq(Eigen::Matrix<double, 1, 1> &delta, double pose_from, double pose_to, double pose_m)\n")
    f_cpp.write("{")
    for i in range (1):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_obs_eq_jacobian(Eigen::Matrix<double, 1, 2, Eigen::RowMajor> &j, double pose_from, double pose_to)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (2):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void relative_pose(Eigen::Matrix<double, 1, 1> &relative_pose, double pose_from, double pose_to)\n")
    f_cpp.write("{")
    for i in range (1):
        f_cpp.write("relative_pose.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(model_function[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
