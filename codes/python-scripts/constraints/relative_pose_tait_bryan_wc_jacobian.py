import sympy
from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

tx_1, ty_1, tz_1 = symbols('tx_1 ty_1 tz_1')
om_1, fi_1, ka_1 = symbols('om_1 fi_1 ka_1')
#sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
#q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
tx_2, ty_2, tz_2 = symbols('tx_2 ty_2 tz_2')
om_2, fi_2, ka_2 = symbols('om_2 fi_2 ka_2')
#sx_2, sy_2, sz_2 = symbols('sx_2 sy_2 sz_2')
#q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')
tx_m, ty_m, tz_m = symbols('tx_m ty_m tz_m')
om_m, fi_m, ka_m = symbols('om_m fi_m ka_m')
#sx_m, sy_m, sz_m = symbols('sx_m sy_m sz_m')
#q0_m, q1_m, q2_m, q3_m = symbols('q0_m q1_m q2_m q3_m')

position_symbols_1 = [tx_1, ty_1, tz_1]
orientation_symbols_1 = [om_1, fi_1, ka_1]
#orientation_symbols_1 = [sx_1, sy_1, sz_1]
#orientation_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
position_symbols_2 = [tx_2, ty_2, tz_2]
orientation_symbols_2 = [om_2, fi_2, ka_2]
#orientation_symbols_2 = [sx_2, sy_2, sz_2]
#orientation_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
all_symbols = position_symbols_1 + orientation_symbols_1 + position_symbols_2 + orientation_symbols_2

Rt_wc_1 = matrix44FromTaitBryan(tx_1, ty_1, tz_1, om_1, fi_1, ka_1)
#Rt_wc_1 = matrix44FromRodrigues(tx_1, ty_1, tz_1, sx_1, sy_1, sz_1)
#Rt_wc_1 = matrix44FromQuaternion(tx_1, ty_1, tz_1, q0_1, q1_1, q2_1, q3_1)
Rt_wc_2 = matrix44FromTaitBryan(tx_2, ty_2, tz_2, om_2, fi_2, ka_2)
#Rt_wc_2 = matrix44FromRodrigues(tx_2, ty_2, tz_2, sx_2, sy_2, sz_2)
#Rt_wc_2 = matrix44FromQuaternion(tx_2, ty_2, tz_2, q0_2, q1_2, q2_2, q3_2)

R_cw_1=Rt_wc_1[:-1,:-1].transpose()
t_wc_1=Matrix([tx_1, ty_1, tz_1]).vec()
t_cw_1=-R_cw_1*t_wc_1
Rt_cw_1=Matrix.hstack(R_cw_1, t_cw_1)
Rt_cw_1=Matrix.vstack(Rt_cw_1, Matrix([[0,0,0,1]]))

relative_pose = Rt_cw_1 * Rt_wc_2

parametrization = taitBryanFromMatrix44Case1(relative_pose)
#parametrization = rodriguesFromMatrix44(relative_pose)
#parametrization = quaternionFromMatrix44(relative_pose)
t = relative_pose[0 : 3, 3]

target_value=Matrix([tx_m, ty_m, tz_m, om_m, fi_m, ka_m])
#target_value=Matrix([tx_m, ty_m, tz_m, sx_m, sy_m, sz_m])
#target_value=Matrix([tx_m, ty_m, tz_m, q0_m, q1_m, q2_m, q3_m])
model_function=t.col_join(parametrization)
delta=target_value-model_function
delta_jacobian=delta.jacobian(all_symbols)

beta_symbols = all_symbols
x_symbols = [tx_m, ty_m, tz_m, om_m, fi_m, ka_m]
sum=Matrix([delta[0,0]*delta[0,0]+delta[1,0]*delta[1,0]+delta[2,0]*delta[2,0]+delta[3,0]*delta[3,0]+delta[4,0]*delta[4,0]+delta[5,0]*delta[5,0]]).vec()
d2sum_dbeta2=sum.jacobian(beta_symbols).jacobian(beta_symbols)
d2sum_dbetadx=sum.jacobian(beta_symbols).jacobian(x_symbols)


##################################################
sin_om_1, cos_om_1 = symbols('sin_om_1 cos_om_1')
sin_fi_1, cos_fi_1 = symbols('sin_fi_1 cos_fi_1')
sin_ka_1, cos_ka_1 = symbols('sin_ka_1 cos_ka_1')

sin_om_2, cos_om_2 = symbols('sin_om_2 cos_om_2')
sin_fi_2, cos_fi_2 = symbols('sin_fi_2 cos_fi_2')
sin_ka_2, cos_ka_2 = symbols('sin_ka_2 cos_ka_2')



### simplified 1
substitutions1 = [(sympy.sin(om_1), sin_om_1), (sympy.cos(om_1), cos_om_1), (sympy.sin(
fi_1), sin_fi_1), (sympy.cos(fi_1), cos_fi_1), (sympy.sin(ka_1), sin_ka_1), (sympy.cos(ka_1), cos_ka_1), (sympy.sin(om_2), sin_om_2), (sympy.cos(om_2), cos_om_2), (sympy.sin(
fi_2), sin_fi_2), (sympy.cos(fi_2), cos_fi_2), (sympy.sin(ka_2), sin_ka_2), (sympy.cos(ka_2), cos_ka_2)]

delta1 = sympy.simplify(delta.subs(substitutions1))
delta1_variables, delta1_simple = sympy.cse(
        delta1, order='none',symbols=symbols('a0:1000'))
delta1_simple = delta1_simple[0]

#print(delta1_variables)
###
delta_jacobian1 = sympy.simplify(delta_jacobian.subs(substitutions1))
delta_jacobian1_variables, delta_jacobian1_simple = sympy.cse(
        delta_jacobian1, order='none')
delta_jacobian1_simple = delta_jacobian1_simple[0]
###
model_function1 = sympy.simplify(model_function.subs(substitutions1))
model_function1_variables, model_function1_simple = sympy.cse(
        model_function1, order='none')
model_function1_simple = model_function1_simple[0]

### simplified 2
substitutions2 = [(sympy.sin(om_1), 0), (sympy.cos(om_1), 1), (sympy.sin(
fi_1), 0), (sympy.cos(fi_1), 1), (sympy.sin(ka_1), 0), (sympy.cos(ka_1), 1), (sympy.sin(om_2), 0), (sympy.cos(om_2), 1), (sympy.sin(
fi_2), 0), (sympy.cos(fi_2), 1), (sympy.sin(ka_2), 0), (sympy.cos(ka_2), 1)]

delta_jacobian2 = sympy.simplify(delta_jacobian.subs(substitutions2))
delta_jacobian2_variables, delta_jacobian2_simple = sympy.cse(
        delta_jacobian2, order='none')
delta_jacobian2_simple = delta_jacobian2_simple[0]

### simplified 3
print("computing AtPA AtPB")
p_x, p_y, p_z, p_om, p_fi, p_ka = symbols('p_x p_y p_z p_om p_fi p_ka')
P=Matrix([[p_x, 0, 0, 0, 0, 0],[0, p_y, 0, 0, 0, 0],[0, 0, p_z, 0, 0, 0],[0, 0, 0, p_om, 0, 0],[0, 0, 0, 0, p_fi, 0],[0, 0, 0, 0, 0, p_ka]])

AtPA = delta_jacobian1_simple.transpose() * P * delta_jacobian1_simple
AtPB = delta_jacobian1_simple.transpose() * P * delta1_simple

#print(AtPA)
#print("----------------------------")
#print(AtPB)

#variables, simple = sympy.cse([delta1, delta_jacobian1], order='none')
#print(variables)

#for name, value_expr in variables:
#    print(name)

#AtPA = delta_jacobian.transpose() * P * delta_jacobian
#AtPB = delta_jacobian.transpose() * P * delta

#AtPA = delta_jacobian1_simple.transpose() * P * delta_jacobian1_simple
#AtPB = delta_jacobian1_simple.transpose() * P * delta1_simple


#print(delta1_variables)
#print(delta_jacobian1_variables)

#AtPA = sympy.simplify(AtPA.subs(substitutions1))
#AtPA_variables, AtPA_simple = sympy.cse(
#        AtPA, order='none')
#AtPA_simple = AtPA_simple[0]

#AtPB = sympy.simplify(AtPB.subs(substitutions1))
#AtPB_variables, AtPB_simple = sympy.cse(
#        AtPB, order='none')
#AtPB_simple = AtPB_simple[0]

with open("relative_pose_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("#ifndef _RELATIVE_POSE_TAIT_BRYAN_WC_JACOBIAN_H_\n")
    f_cpp.write("#define _RELATIVE_POSE_TAIT_BRYAN_WC_JACOBIAN_H_\n")
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1(Eigen::Matrix<double, 6, 1> &delta, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double tx_m, double ty_m, double tz_m, double om_m, double fi_m, double ka_m)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1_jacobian(Eigen::Matrix<double, 6, 12, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void relative_pose_tait_bryan_wc_case1(Eigen::Matrix<double, 6, 1> &relative_pose, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("relative_pose.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(model_function[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1_d2sum_dbeta2(Eigen::Matrix<double, 12, 12, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double tx_m, double ty_m, double tz_m, double om_m, double fi_m, double ka_m)\n")
    f_cpp.write("{")
    for i in range (12):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbeta2[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1_d2sum_dbetadx(Eigen::Matrix<double, 12, 6, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double tx_m, double ty_m, double tz_m, double om_m, double fi_m, double ka_m)\n")
    f_cpp.write("{")
    for i in range (12):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbetadx[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")
    ########################################### _simplified_1 #######################################
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1_simplified_1(Eigen::Matrix<double, 6, 1> &delta, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double tx_m, double ty_m, double tz_m, double om_m, double fi_m, double ka_m)\n")
    f_cpp.write("{\n")
    f_cpp.write("double sin_om_1 = sin(om_1);\n")
    f_cpp.write("double cos_om_1 = cos(om_1);\n")
    f_cpp.write("double sin_fi_1 = sin(fi_1);\n")
    f_cpp.write("double cos_fi_1 = cos(fi_1);\n")
    f_cpp.write("double sin_ka_1 = sin(ka_1);\n")
    f_cpp.write("double cos_ka_1 = cos(ka_1);\n")
    f_cpp.write("double sin_om_2 = sin(om_2);\n")
    f_cpp.write("double cos_om_2 = cos(om_2);\n")
    f_cpp.write("double sin_fi_2 = sin(fi_2);\n")
    f_cpp.write("double cos_fi_2 = cos(fi_2);\n")
    f_cpp.write("double sin_ka_2 = sin(ka_2);\n")
    f_cpp.write("double cos_ka_2 = cos(ka_2);\n")
    for name, value_expr in delta1_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    for i in range (6):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta1_simple[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1_jacobian_simplified_1(Eigen::Matrix<double, 6, 12, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{\n")
    f_cpp.write("double sin_om_1 = sin(om_1);\n")
    f_cpp.write("double cos_om_1 = cos(om_1);\n")
    f_cpp.write("double sin_fi_1 = sin(fi_1);\n")
    f_cpp.write("double cos_fi_1 = cos(fi_1);\n")
    f_cpp.write("double sin_ka_1 = sin(ka_1);\n")
    f_cpp.write("double cos_ka_1 = cos(ka_1);\n")
    f_cpp.write("double sin_om_2 = sin(om_2);\n")
    f_cpp.write("double cos_om_2 = cos(om_2);\n")
    f_cpp.write("double sin_fi_2 = sin(fi_2);\n")
    f_cpp.write("double cos_fi_2 = cos(fi_2);\n")
    f_cpp.write("double sin_ka_2 = sin(ka_2);\n")
    f_cpp.write("double cos_ka_2 = cos(ka_2);\n")
    for name, value_expr in delta_jacobian1_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    for i in range (6):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian1_simple[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_tait_bryan_wc_case1_simplified_1(Eigen::Matrix<double, 6, 1> &relative_pose, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{")
    f_cpp.write("double sin_om_1 = sin(om_1);\n")
    f_cpp.write("double cos_om_1 = cos(om_1);\n")
    f_cpp.write("double sin_fi_1 = sin(fi_1);\n")
    f_cpp.write("double cos_fi_1 = cos(fi_1);\n")
    f_cpp.write("double sin_ka_1 = sin(ka_1);\n")
    f_cpp.write("double cos_ka_1 = cos(ka_1);\n")
    f_cpp.write("double sin_om_2 = sin(om_2);\n")
    f_cpp.write("double cos_om_2 = cos(om_2);\n")
    f_cpp.write("double sin_fi_2 = sin(fi_2);\n")
    f_cpp.write("double cos_fi_2 = cos(fi_2);\n")
    f_cpp.write("double sin_ka_2 = sin(ka_2);\n")
    f_cpp.write("double cos_ka_2 = cos(ka_2);\n")
    for name, value_expr in model_function1_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    for i in range (6):
        f_cpp.write("relative_pose.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(model_function1_simple[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    ########################################### _simplified_2 #######################################
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1_jacobian_simplified_2(Eigen::Matrix<double, 6, 12, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{\n")
    for name, value_expr in delta_jacobian2_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    for i in range (6):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian2_simple[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")
    ########################################### _simplified_3 #######################################
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1_AtPA_simplified(Eigen::Matrix<double, 12, 12> &AtPA, const double &tx_1, const double &ty_1, const double &tz_1, const double &om_1, const double &fi_1, const double &ka_1, const double &tx_2, const double &ty_2, const double &tz_2, const double &om_2, const double &fi_2, const double &ka_2, const double &p_x, const double &p_y, const double &p_z, const double &p_om, const double &p_fi, const double &p_ka)\n")
    f_cpp.write("{\n")
    f_cpp.write("double sin_om_1 = sin(om_1);\n")
    f_cpp.write("double cos_om_1 = cos(om_1);\n")
    f_cpp.write("double sin_fi_1 = sin(fi_1);\n")
    f_cpp.write("double cos_fi_1 = cos(fi_1);\n")
    f_cpp.write("double sin_ka_1 = sin(ka_1);\n")
    f_cpp.write("double cos_ka_1 = cos(ka_1);\n")
    f_cpp.write("double sin_om_2 = sin(om_2);\n")
    f_cpp.write("double cos_om_2 = cos(om_2);\n")
    f_cpp.write("double sin_fi_2 = sin(fi_2);\n")
    f_cpp.write("double cos_fi_2 = cos(fi_2);\n")
    f_cpp.write("double sin_ka_2 = sin(ka_2);\n")
    f_cpp.write("double cos_ka_2 = cos(ka_2);\n")
    for name, value_expr in delta_jacobian1_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    for i in range (12):
        for j in range (12):
            f_cpp.write("AtPA.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(AtPA[i,j])))
    f_cpp.write("}\n")
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_wc_case1_AtPB_simplified(Eigen::Matrix<double, 12, 1> &AtPB, const double &tx_1, const double &ty_1, const double &tz_1, const double &om_1, const double &fi_1, const double &ka_1, const double &tx_2, const double &ty_2, const double &tz_2, const double &om_2, const double &fi_2, const double &ka_2, const double &tx_m, const double &ty_m, const double &tz_m, const double &om_m, const double &fi_m, const double &ka_m, const double &p_x, const double &p_y, const double &p_z, const double &p_om, const double &p_fi, const double &p_ka)\n")
    f_cpp.write("{\n")
    f_cpp.write("double sin_om_1 = sin(om_1);\n")
    f_cpp.write("double cos_om_1 = cos(om_1);\n")
    f_cpp.write("double sin_fi_1 = sin(fi_1);\n")
    f_cpp.write("double cos_fi_1 = cos(fi_1);\n")
    f_cpp.write("double sin_ka_1 = sin(ka_1);\n")
    f_cpp.write("double cos_ka_1 = cos(ka_1);\n")
    f_cpp.write("double sin_om_2 = sin(om_2);\n")
    f_cpp.write("double cos_om_2 = cos(om_2);\n")
    f_cpp.write("double sin_fi_2 = sin(fi_2);\n")
    f_cpp.write("double cos_fi_2 = cos(fi_2);\n")
    f_cpp.write("double sin_ka_2 = sin(ka_2);\n")
    f_cpp.write("double cos_ka_2 = cos(ka_2);\n")
    for name, value_expr in delta1_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    for name, value_expr in delta_jacobian1_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    for i in range (12):
        f_cpp.write("AtPB.coeffRef(%d) = %s;\n"%(i, ccode(AtPB[i])))
    f_cpp.write("}\n")
    f_cpp.write("#endif\n")

