from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_t, y_t, z_t = symbols('x_t y_t z_t')
x_s, y_s, z_s = symbols('x_s y_s z_s')
tx, ty, tz = symbols('tx ty tz')
om, fi, ka = symbols('om fi ka')
ctx, cty, ctz = symbols('ctx cty ctz')
com, cfi, cka = symbols('com cfi cka')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

cal_position_symbols = [ctx, cty, ctz]
cal_orientation_symbols = [com, cfi, cka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
all_symbols = cal_position_symbols + cal_orientation_symbols

point_t = Matrix([x_t, y_t, z_t]).vec()
point_s = Matrix([x_s, y_s, z_s, 1]).vec()
transformed_point_s = (matrix44FromTaitBryan(tx, ty, tz, om, fi, ka) * matrix44FromTaitBryan(ctx, cty, ctz, com, cfi, cka) * point_s)[:-1,:]
#transformed_point_s = (matrix44FromRodrigues(tx, ty, tz, sx, sy, sz) * point_s)[:-1,:]
#transformed_point_s = (matrix44FromQuaternion(tx, ty, tz, q0, q1, q2, q3) * point_s)[:-1,:]

target_value = Matrix([0,0,0]).vec()
model_function = transformed_point_s-point_t
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)


with open("point_to_point_source_to_target_calibration_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_source_to_target_calibration_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, double tx, double ty, double tz, double om, double fi, double ka, double ctx, double cty, double ctz, double com, double cfi, double cka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_target_calibration_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double ctx, double cty, double ctz, double com, double cfi, double cka, double x_s, double y_s, double z_s)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
