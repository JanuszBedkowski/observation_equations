from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_1, y_1, z_1 = symbols('x_1 y_1 z_1')
s_1 = symbols('s_1')
tx_1, ty_1, tz_1 = symbols('tx_1 ty_1 tz_1')
om_1, fi_1, ka_1 = symbols('om_1 fi_1 ka_1')
tx_i, ty_i, tz_i = symbols('tx_i ty_i tz_i')
om_i, fi_i, ka_i = symbols('om_i fi_i ka_i')
tx_c, ty_c, tz_c = symbols('tx_c ty_c tz_c')
om_c, fi_c, ka_c = symbols('om_c fi_c ka_c')


#sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
#q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
x_t, y_t, z_t = symbols('x_t y_t z_t')
#s_2 = symbols('s_2')
#tx_2, ty_2, tz_2 = symbols('tx_2 ty_2 tz_2')
#om_2, fi_2, ka_2 = symbols('om_2 fi_2 ka_2')
#sx_2, sy_2, sz_2 = symbols('sx_2 sy_2 sz_2')
#q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')

position_symbols_1 = [tx_1, ty_1, tz_1]
orientation_symbols_1 = [om_1, fi_1, ka_1]
position_symbols_i = [tx_i, ty_i, tz_i]
orientation_symbols_i = [om_i, fi_i, ka_i]
position_symbols_c = [tx_c, ty_c, tz_c]
orientation_symbols_c = [om_c, fi_c, ka_c]

#orientation_symbols_1 = [sx_1, sy_1, sz_1]
#orientation_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
#position_symbols_2 = [tx_2, ty_2, tz_2]
#orientation_symbols_2 = [om_2, fi_2, ka_2]
#orientation_symbols_2 = [sx_2, sy_2, sz_2]
#orientation_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
all_symbols = position_symbols_i + orientation_symbols_i + position_symbols_1 + orientation_symbols_1 + position_symbols_c + orientation_symbols_c + [s_1]

point_1 = Matrix([x_1, y_1, z_1, 1]).vec()
point_2 = Matrix([x_t, y_t, z_t]).vec()
m_scale_1 = Matrix([[s_1, 0, 0, 0], [0, s_1, 0, 0], [0, 0, s_1, 0], [0,0,0,1]])
#m_scale_2 = Matrix([[s_2, 0, 0, 0], [0, s_2, 0, 0], [0, 0, s_2, 0], [0,0,0,1]])

transformed_point_1 = ( matrix44FromTaitBryan(tx_i, ty_i, tz_i, om_i, fi_i, ka_i) * 
matrix44FromTaitBryan(tx_1, ty_1, tz_1, om_1, fi_1, ka_1) * matrix44FromTaitBryan(tx_c, ty_c, tz_c, om_c, fi_c, ka_c) *  m_scale_1 * point_1)[:-1,:]
#transformed_point_1 = ( matrix44FromRodrigues(tx_1, ty_1, tz_1, sx_1, sy_1, sz_1) * m_scale_1 * point_1)[:-1,:]
#transformed_point_1 = ( matrix44FromQuaternion(tx_1, py_1, tz_1, q0_1, q1_1, q2_1, q3_1) * m_scale_1 * point_1)[:-1,:]
transformed_point_2 = point_2
#transformed_point_2 = ( matrix44FromRodrigues(tx_2, ty_2, tz_2, sx_2, sy_2, sz_2) * m_scale_2 * point_2)[:-1,:]
#transformed_point_2 = ( matrix44FromQuaternion(tx_2, ty_2, tz_2, q0_2, q1_2, q2_2, q3_2) * m_scale_2 * point_2)[:-1,:]

target_value = Matrix([0,0,0]).vec()
model_function = transformed_point_1-transformed_point_2
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_point_with_scale_hpp_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_with_scale_hpp_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, double tx_i, double ty_i, double tz_i, double om_i, double fi_i, double ka_i, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_c, double ty_c, double tz_c, double om_c, double fi_c, double ka_c, double x_1, double y_1, double z_1, double x_t, double y_t, double z_t, double s_1)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_with_scale_hpp_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 19, Eigen::RowMajor> &j, double tx_i, double ty_i, double tz_i, double om_i, double fi_i, double ka_i, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_c, double ty_c, double tz_c, double om_c, double fi_c, double ka_c, double x_1, double y_1, double z_1, double x_t, double y_t, double z_t, double s_1)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (19):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")



