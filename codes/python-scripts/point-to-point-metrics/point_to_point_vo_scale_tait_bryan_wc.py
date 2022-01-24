from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_1, y_1, z_1 = symbols('x_1 y_1 z_1')
tx_1, ty_1, tz_1 = symbols('tx_1 ty_1 tz_1')
om_1, fi_1, ka_1 = symbols('om_1 fi_1 ka_1')
#sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
#q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
x_2, y_2, z_2 = symbols('x_2 y_2 z_2')
tx_2, ty_2, tz_2 = symbols('tx_2 ty_2 tz_2')
om_2, fi_2, ka_2 = symbols('om_2 fi_2 ka_2')
#sx_2, sy_2, sz_2 = symbols('sx_2 sy_2 sz_2')
#q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')
scale = symbols('scale')
all_symbols = [scale]

point_1 = Matrix([x_1, y_1, z_1, 1]).vec()
point_2 = Matrix([x_2, y_2, z_2, 1]).vec()

transformed_point_1 = (matrix44FromTaitBryan(scale * tx_1, scale * ty_1, scale * tz_1, om_1, fi_1, ka_1) * point_1)[:-1,:]
transformed_point_2 = (matrix44FromTaitBryan(scale * tx_2, scale * ty_2, scale * tz_2, om_2, fi_2, ka_2) * point_2)[:-1,:]
#transformed_point_1 = (matrix44FromRodrigues(scale * tx_1, scale * ty_1, scale * tz_1, sx_1, sy_1, sz_1) * point_1)[:-1,:]
#transformed_point_2 = (matrix44FromRodrigues(scale * tx_2, scale * ty_2, scale * tz_2, sx_2, sy_2, sz_2) * point_2)[:-1,:]
#transformed_point_1 = (matrix44FromQuaternion(scale * tx_1, scale * ty_1, scale * tz_1, q0_1, q1_1, q2_1, q3_1) * point_1)[:-1,:]
#transformed_point_2 = (matrix44FromQuaternion(scale * tx_2, scale * ty_2, scale * tz_2, q0_1, q1_1, q2_1, q3_1) * point_2)[:-1,:]

target_value = Matrix([0,0,0]).vec()
model_function = transformed_point_1-transformed_point_2
delta = target_value-model_function
delta_jacobian = delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_point_vo_scale_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_vo_scale_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, double scale)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_vo_scale_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 1> &j, double tx_1, double ty_1, double tz_1, double tx_2, double ty_2, double tz_2)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (1):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")



