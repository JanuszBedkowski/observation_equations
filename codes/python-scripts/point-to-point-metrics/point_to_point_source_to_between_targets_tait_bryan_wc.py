from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_t1, y_t1, z_t1 = symbols('x_t1 y_t1 z_t1')
x_t2, y_t2, z_t2 = symbols('x_t2 y_t2 z_t2')
x_s, y_s, z_s = symbols('x_s y_s z_s')
tx, ty, tz = symbols('tx ty tz')
om, fi, ka = symbols('om fi ka')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

position_symbols = [tx, ty, tz]
orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + orientation_symbols

point_t1 = Matrix([x_t1, y_t1, z_t1]).vec()
point_t2 = Matrix([x_t2, y_t2, z_t2]).vec()
point_s = Matrix([x_s, y_s, z_s, 1]).vec()
transformed_point_s = (matrix44FromTaitBryan(tx, ty, tz, om, fi, ka) * point_s)[:-1,:]
#transformed_point_s = (matrix44FromRodrigues(tx, ty, tz, sx, sy, sz) * point_s)[:-1,:]
#transformed_point_s = (matrix44FromQuaternion(tx, ty, tz, q0, q1, q2, q3) * point_s)[:-1,:]

target_value = Matrix([0,0,0]).vec()
model_function = transformed_point_s - point_t1 - (point_t2 - transformed_point_s)
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

#beta_symbols = position_symbols + orientation_symbols
#x_symbols = [x_s, y_s, z_s, x_t, y_t, z_t]
#sum=Matrix([delta[0,0]*delta[0,0]+delta[1,0]*delta[1,0]+delta[2,0]*delta[2,0]]).vec()
#d2sum_dbeta2=sum.jacobian(beta_symbols).jacobian(beta_symbols)
#d2sum_dbetadx=sum.jacobian(beta_symbols).jacobian(x_symbols)

with open("point_to_point_source_to_between_targets_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_source_to_between_targets_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t1, double y_t1, double z_t1, double x_t2, double y_t2, double z_t2)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_between_targets_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    #f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_d2sum_dbeta2(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    #f_cpp.write("{")
    #for i in range (6):
    #    for j in range (6):
    #        f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbeta2[i,j])))
    #f_cpp.write("}")
    #f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_d2sum_dbetadx(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    #f_cpp.write("{")
    #for i in range (6):
    #    for j in range (6):
    #        f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbetadx[i,j])))
    #f_cpp.write("}")


