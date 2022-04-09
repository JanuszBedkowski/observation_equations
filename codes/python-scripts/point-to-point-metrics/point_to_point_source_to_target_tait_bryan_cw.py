from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_t, y_t, z_t = symbols('x_t y_t z_t')
x_s, y_s, z_s = symbols('x_s y_s z_s')
tx_cw, ty_cw, tz_cw = symbols('tx_cw ty_cw tz_cw')
om_cw, fi_cw, ka_cw = symbols('om_cw fi_cw ka_cw')
#sx_cw, sy_cw, sz_cw = symbols('sx sy sz')
#q0_cw, q1_cw, q2_cw, q3_cw = symbols('q0_cw q1_cw q2_cw q3_cw')

position_symbols = [tx_cw, ty_cw, tz_cw]
orientation_symbols = [om_cw, fi_cw, ka_cw]
#orientation_symbols = [sx_cw, sy_cw, sz_cw]
#orientation_symbols = [q0_cw, q1_cw, q2_cw, q3_cw]
all_symbols = position_symbols + orientation_symbols

point_t = Matrix([x_t, y_t, z_t]).vec()
point_s = Matrix([x_s, y_s, z_s, 1]).vec()

RT_cw = matrix44FromTaitBryan(tx_cw, ty_cw, tz_cw, om_cw, fi_cw, ka_cw)
#RT_cw = matrix44FromRodrigues(tx_cw, ty_cw, tz_cw, sx_cw, sy_cw, sz_cw)
#RT_cw = matrix44FromQuaternion(tx_cw, ty_cw, tz_cw, q0_cw, q1_cw, q2_cw, q3_cw)
R_wc=RT_cw[:-1,:-1].transpose()
T_cw=Matrix([tx_cw, ty_cw, tz_cw]).vec()
T_wc=-R_wc*T_cw
RT_wc=Matrix.hstack(R_wc, T_wc)
RT_wc=Matrix.vstack(RT_wc, Matrix([[0,0,0,1]]))

transformed_point_s = (RT_wc * point_s)[:-1,:]
#transformed_point_s = (RT_wc * point_s)[:-1,:]
#transformed_point_s = (RT_wc * point_s)[:-1,:]

target_value = Matrix([0,0,0]).vec()
model_function = transformed_point_s - point_t
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

beta_symbols = position_symbols
x_symbols = [x_t, y_t, z_t]
sum=Matrix([delta[0,0]*delta[0,0]+delta[1,0]*delta[1,0]+delta[2,0]*delta[2,0]]).vec()
d2sum_dbeta2=sum.jacobian(beta_symbols).jacobian(beta_symbols)
d2sum_dbetadx=sum.jacobian(beta_symbols).jacobian(x_symbols)

with open("point_to_point_source_to_target_tait_bryan_cw_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_cw(double &delta_x, double &delta_y, double &delta_z, double tx_cw, double ty_cw, double tz_cw, double om_cw, double fi_cw, double ka_cw, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_cw_jacobian(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double tx_cw, double ty_cw, double tz_cw, double om_cw, double fi_cw, double ka_cw, double x_s, double y_s, double z_s)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_cw_d2sum_dbeta2(Eigen::Matrix<double, 3, 3, Eigen::RowMajor> &j, double tx_cw, double ty_cw, double tz_cw, double om_cw, double fi_cw, double ka_cw, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (3):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbeta2[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_cw_d2sum_dbetadx(Eigen::Matrix<double, 3, 3, Eigen::RowMajor> &j, double tx_cw, double ty_cw, double tz_cw, double om_cw, double fi_cw, double ka_cw, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (3):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbetadx[i,j])))
    f_cpp.write("}")


