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

Rt_cw_1 = matrix44FromTaitBryan(tx_1, ty_1, tz_1, om_1, fi_1, ka_1)
#Rt_wc_1 = matrix44FromRodrigues(tx_1, ty_1, tz_1, sx_1, sy_1, sz_1)
#Rt_wc_1 = matrix44FromQuaternion(tx_1, ty_1, tz_1, q0_1, q1_1, q2_1, q3_1)
Rt_cw_2 = matrix44FromTaitBryan(tx_2, ty_2, tz_2, om_2, fi_2, ka_2)
#Rt_wc_2 = matrix44FromRodrigues(tx_2, ty_2, tz_2, sx_2, sy_2, sz_2)
#Rt_wc_2 = matrix44FromQuaternion(tx_2, ty_2, tz_2, q0_2, q1_2, q2_2, q3_2)

#R_cw_1=Rt_wc_1[:-1,:-1].transpose()
#t_wc_1=Matrix([tx_1, ty_1, tz_1]).vec()
#t_cw_1=-R_cw_1*t_wc_1
#Rt_cw_1=Matrix.hstack(R_cw_1, t_cw_1)
#Rt_cw_1=Matrix.vstack(Rt_cw_1, Matrix([[0,0,0,1]]))
R_wc_2=Rt_cw_2[:-1,:-1].transpose()
t_cw_2=Matrix([tx_2, ty_2, tz_2]).vec()
t_wc_2=-R_wc_2*t_cw_2
Rt_wc_2=Matrix.hstack(R_wc_2, t_wc_2)
Rt_wc_2=Matrix.vstack(Rt_wc_2, Matrix([[0,0,0,1]]))

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

#beta_symbols = all_symbols
#x_symbols = [tx_m, ty_m, tz_m, om_m, fi_m, ka_m]
#sum=Matrix([delta[0,0]*delta[0,0]+delta[1,0]*delta[1,0]+delta[2,0]*delta[2,0]+delta[3,0]*delta[3,0]+delta[4,0]*delta[4,0]+delta[5,0]*delta[5,0]]).vec()
#d2sum_dbeta2=sum.jacobian(beta_symbols).jacobian(beta_symbols)
#d2sum_dbetadx=sum.jacobian(beta_symbols).jacobian(x_symbols)

with open("relative_pose_tait_bryan_cw_jacobian.h",'w') as f_cpp:  
    f_cpp.write("#ifndef _relative_pose_tait_bryan_cw_jacobian_h_")
    f_cpp.write("\n")
    f_cpp.write("#define _relative_pose_tait_bryan_cw_jacobian_h_")
    f_cpp.write("\n")  
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_cw_case1(Eigen::Matrix<double, 6, 1> &delta, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double tx_m, double ty_m, double tz_m, double om_m, double fi_m, double ka_m)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_cw_case1_jacobian(Eigen::Matrix<double, 6, 12, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void relative_pose_tait_bryan_cw_case1(Eigen::Matrix<double, 6, 1> &relative_pose, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("relative_pose.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(model_function[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
#    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_cw_case1_d2sum_dbeta2(Eigen::Matrix<double, 12, 12, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double tx_m, double ty_m, double tz_m, double om_m, double fi_m, double ka_m)\n")
#    f_cpp.write("{")
#    for i in range (12):
#        for j in range (12):
#            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbeta2[i,j])))
#    f_cpp.write("}")
#    f_cpp.write("inline void relative_pose_obs_eq_tait_bryan_cw_case1_d2sum_dbetadx(Eigen::Matrix<double, 12, 6, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double tx_m, double ty_m, double tz_m, double om_m, double fi_m, double ka_m)\n")
#    f_cpp.write("{")
#    for i in range (12):
#        for j in range (6):
#            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbetadx[i,j])))
#    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("#endif")


