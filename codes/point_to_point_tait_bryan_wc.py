from sympy import *
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_1, y_1, z_1 = symbols('x_1 y_1 z_1')
px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
om_1, fi_1, ka_1 = symbols('om_1 fi_1 ka_1')
#sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
#q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
x_2, y_2, z_2 = symbols('x_2 y_2 z_2')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
om_2, fi_2, ka_2 = symbols('om_2 fi_2 ka_2')
#sx_2, sy_2, sz_2 = symbols('sx_2 sy_2 sz_2')
#q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')

position_symbols_1 = [px_1, py_1, pz_1]
orientation_symbols_1 = [om_1, fi_1, ka_1]
#orientation_symbols_1 = [sx_1, sy_1, sz_1]
#orientation_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
position_symbols_2 = [px_2, py_2, pz_2]
orientation_symbols_2 = [om_2, fi_2, ka_2]
#orientation_symbols_2 = [sx_2, sy_2, sz_2]
#orientation_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
all_symbols = position_symbols_1 + orientation_symbols_1 + position_symbols_2 + orientation_symbols_2

point_1 = Matrix([x_1, y_1, z_1, 1]).vec()
point_2 = Matrix([x_2, y_2, z_2, 1]).vec()

transformed_point_1 = (matrix44FromTaitBryan(px_1, py_1, pz_1, om_1, fi_1, ka_1) * point_1)[:-1,:]
#transformed_point_1 = (matrix44FromRodrigues(px_1, py_1, pz_1, sx_1, sy_1, sz_1) * point_1)[:-1,:]
#transformed_point_1 = (matrix44FromQuaternion(px_1, py_1, pz_1, q0_1, q1_1, q2_1, q3_1) * point_1)[:-1,:]
transformed_point_2 = (matrix44FromTaitBryan(px_2, py_2, pz_2, om_2, fi_2, ka_2) * point_2)[:-1,:]
#transformed_point_2 = (matrix44FromRodrigues(px_2, py_2, pz_2, sx_2, sy_2, sz_2) * point_2)[:-1,:]
#transformed_point_2 = (matrix44FromQuaternion(px_2, py_2, pz_2, q0_2, q1_2, q2_2, q3_2) * point_2)[:-1,:]

target_value = Matrix([0,0,0]).vec()
model_function = transformed_point_1-transformed_point_2
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_point_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 12, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")



