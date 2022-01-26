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
tx_3, ty_3, tz_3 = symbols('tx_3 ty_3 tz_3')
om_3, fi_3, ka_3 = symbols('om_3 fi_3 ka_3')
#sx_3, sy_3, sz_3 = symbols('sx_3 sy_3 sz_3')
#q0_3, q1_3, q2_3, q3_3 = symbols('q0_3 q1_3 q2_3 q3_3')

position_symbols_1 = [tx_1, ty_1, tz_1]
orientation_symbols_1 = [om_1, fi_1, ka_1]
#orientation_symbols_1 = [sx_1, sy_1, sz_1]
#orientation_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
position_symbols_2 = [tx_2, ty_2, tz_2]
orientation_symbols_2 = [om_2, fi_2, ka_2]
#orientation_symbols_2 = [sx_2, sy_2, sz_2]
#orientation_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
position_symbols_3 = [tx_3, ty_3, tz_3]
orientation_symbols_3 = [om_3, fi_3, ka_3]
#rodrigues_symbols_3 = [sx_3, sy_3, sz_3]
#quaternion_symbols_3 = [q0_3, q1_3, q2_3, q3_3]
all_symbols = position_symbols_1 + orientation_symbols_1 + position_symbols_2 + orientation_symbols_2 + position_symbols_3 + orientation_symbols_3

Rt_wc_1 = matrix44FromTaitBryan(tx_1, ty_1, tz_1, om_1, fi_1, ka_1)
#Rt_wc_1 = matrix44FromRodrigues(tx_1, ty_1, tz_1, sx_1, sy_1, sz_1)
#Rt_wc_1 = matrix44FromQuaternion(tx_1,ty_1,tz_1,q0_1,q1_1,q2_1,q3_1)
Rt_wc_2 = matrix44FromTaitBryan(tx_2, ty_2, tz_2, om_2, fi_2, ka_2)
#Rt_wc_2 = matrix44FromRodrigues(tx_2, ty_2, tz_2, sx_2, sy_2, sz_2)
#Rt_wc_2 = matrix44FromQuaternion(tx_2,ty_2,tz_2,q0_2,q1_2,q2_2,q3_2)
Rt_wc_3 = matrix44FromTaitBryan(tx_3, ty_3, tz_3, om_3, fi_3, ka_3)
#Rt_wc_3 = matrix44FromRodrigues(tx_3, ty_3, tz_3, sx_3, sy_3, sz_3)
#Rt_wc_3 = matrix44FromQuaternion(tx_3,ty_3,tz_3,q0_3,q1_3,q2_3,q3_3)

R_cw_1=Rt_wc_1[:-1,:-1].transpose()
t_wc_1=Matrix([tx_1, ty_1, tz_1]).vec()
t_cw_1=-R_cw_1*t_wc_1
Rt_cw_1=Matrix.hstack(R_cw_1, t_cw_1)
Rt_cw_1=Matrix.vstack(Rt_cw_1, Matrix([[0,0,0,1]]))

R_cw_2=Rt_wc_2[:-1,:-1].transpose()
t_wc_2=Matrix([tx_2, ty_2, tz_2]).vec()
t_cw_2=-R_cw_2*t_wc_2
Rt_cw_2=Matrix.hstack(R_cw_2, t_cw_2)
Rt_cw_2=Matrix.vstack(Rt_cw_2, Matrix([[0,0,0,1]]))

relative_pose12 = Rt_cw_1 * Rt_wc_2
relative_pose23 = Rt_cw_2 * Rt_wc_3

parametrization12=taitBryanFromMatrix44Case1(relative_pose12)
#parametrization12=rodriguesFromMatrix44(relative_pose12)
#parametrization12=quaternionFromMatrix44(relative_pose12)
parametrization23=taitBryanFromMatrix44Case1(relative_pose23)
#parametrization23=rodriguesFromMatrix44(relative_pose23)
#parametrization23=quaternionFromMatrix44(relative_pose23)
t12 = relative_pose12[0 : 3, 3]
t23 = relative_pose23[0 : 3, 3]

rel_pose12=t12.col_join(parametrization12)
rel_pose23=t23.col_join(parametrization23)

target_value = Matrix([0, 0, 0, 0, 0, 0])
model_function = rel_pose12-rel_pose23

delta = target_value-model_function
delta_jacobian = delta.jacobian(all_symbols)


with open("smoothness_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void smoothness_obs_eq_tait_bryan_wc(Eigen::Matrix<double, 6, 1> &delta, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double tx_3, double ty_3, double tz_3, double om_3, double fi_3, double ka_3)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void smoothness_obs_eq_tait_bryan_wc_jacobian(Eigen::Matrix<double, 6, 18, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double tx_3, double ty_3, double tz_3, double om_3, double fi_3, double ka_3)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (18):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    
