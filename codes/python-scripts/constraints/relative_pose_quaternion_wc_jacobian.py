import sympy
from sympy import *
import sys
sys.path.insert(1, '..')
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

### simplified 1
#print("computing delta1_simple")
#substitutions1 = []

#delta1 = sympy.simplify(delta.subs(substitutions1))
#delta1_variables, delta1_simple = sympy.cse(
#        delta1, order='none',symbols=symbols('a0:1000'))
#delta1_simple = delta1_simple[0]

###
#print("computing delta_jacobian1_simple")
#delta_jacobian1 = sympy.simplify(delta_jacobian.subs(substitutions1))
#delta_jacobian1_variables, delta_jacobian1_simple = sympy.cse(
#        delta_jacobian1, order='none')
#delta_jacobian1_simple = delta_jacobian1_simple[0]
###

#print("computing model_function1_simple")
#model_function1 = sympy.simplify(model_function.subs(substitutions1))
#model_function1_variables, model_function1_simple = sympy.cse(
#        model_function1, order='none')
#model_function1_simple = model_function1_simple[0]

### simplified 3
#print("computing AtPA AtPB")
#p_x, p_y, p_z, p_q0, p_q1, p_q2, p_q3 = symbols('p_x p_y p_z p_q0 p_q1 p_q2 p_q3')
#P=Matrix([[p_x, 0, 0, 0, 0, 0, 0],[0, p_y, 0, 0, 0, 0, 0],[0, 0, p_z, 0, 0, 0, 0],[0, 0, 0, p_q0, 0, 0, 0],[0, 0, 0, 0, p_q1, 0, 0],[0, 0, 0, 0, 0, p_q2, 0], [0, 0, 0, 0, 0, 0, p_q3]])

#AtPA = delta_jacobian1_simple.transpose() * P * delta_jacobian1_simple
#AtPB = delta_jacobian1_simple.transpose() * P * delta1_simple



with open("relative_pose_quaternion_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("#ifndef _RELATIVE_POSE_QUATERNION_WC_JACOBIAN_H_\n")
    f_cpp.write("#define _RELATIVE_POSE_QUATERNION_WC_JACOBIAN_H_\n")
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
    f_cpp.write("}\n")
########################################### _simplified_3 #######################################
    #f_cpp.write("inline void relative_pose_obs_eq_quaternion_wc_case1_AtPA_simplified(Eigen::Matrix<double, 14, 14> &AtPA, const double &tx_1, const double &ty_1, const double &tz_1, const double &q0_1, const double &q1_1, const double &q2_1, const double &q3_1, const double &tx_2, const double &ty_2, const double &tz_2, const double &q0_2, const double &q1_2, const double &q2_2, const double &q3_2, const double &p_x, const double &p_y, const double &p_z, const double &p_q0, const double &p_q1, const double &p_q2, const double &p_q3)\n")
    #f_cpp.write("{\n")
    #for name, value_expr in delta_jacobian1_variables:
    #    f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    #for i in range (14):
    #    for j in range (14):
    #        f_cpp.write("AtPA.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(AtPA[i,j])))
    #f_cpp.write("}\n")
    #f_cpp.write("inline void relative_pose_obs_eq_quaternion_wc_case1_AtPB_simplified(Eigen::Matrix<double, 14, 1> &AtPB, const double &tx_1, const double &ty_1, const double &tz_1, const double &q0_1, const double &q1_1, const double &q2_1, const double &q3_1, const double &tx_2, const double &ty_2, const double &tz_2, const double &q0_2, const double &q1_2, const double &q2_2, const double &q3_2, const double &tx_m, const double &ty_m, const double &tz_m, const double &q0_m, const double &q1_m, const double &q2_m, const double &q3_m, const double &p_x, const double &p_y, const double &p_z, const double &p_q0, const double &p_q1, const double &p_q2, const double &p_q3)\n")
    #f_cpp.write("{\n")
    #for name, value_expr in delta1_variables:
    #    f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    #for name, value_expr in delta_jacobian1_variables:
    #    f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    #for i in range (12):
    #    f_cpp.write("AtPB.coeffRef(%d) = %s;\n"%(i, ccode(AtPB[i])))
    #f_cpp.write("}\n")
    f_cpp.write("#endif\n")

