from sympy import *
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
om_1, fi_1, ka_1 = symbols('om_1 fi_1 ka_1')
#sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
#q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
om_2, fi_2, ka_2 = symbols('om_2 fi_2 ka_2')
#sx_2, sy_2, sz_2 = symbols('sx_2 sy_2 sz_2')
#q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')
rows_1, cols_1, rows_2, cols_2 = symbols('rows_1 cols_1 rows_2 cols_2')
u_1, v_1, u_2, v_2 = symbols('u_1 v_1 u_2 v_2')
pi = symbols('pi')

position_symbols_1 = [px_1, py_1, pz_1]
orientation_symbols_1 = [om_1, fi_1, ka_1]
#orientation_symbols_1 = [sx_1, sy_1, sz_1]
#orientation_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
position_symbols_2 = [px_2, py_2, pz_2]
orientation_symbols_2 = [om_2, fi_2, ka_2]
#orientation_symbols_2 = [sx_2, sy_2, sz_2]
#orientation_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
all_symbols = position_symbols_1 + orientation_symbols_1 + position_symbols_2 + orientation_symbols_2

camera_matrix_1 = matrix44FromTaitBryan(px_1, py_1, pz_1, om_1, fi_1, ka_1)
#camera_matrix_1 = matrix44FromRodrigues(px_1, py_1, pz_1, sx_1, sy_1, sz_1)
#camera_matrix_1 = matrix44FromQuaternion(px_1, py_1, pz_1, q0_1, q1_1, q2_1, q3_1)
R_1t=camera_matrix_1[:-1,:-1].transpose()
camera_matrix_2 = matrix44FromTaitBryan(px_2, py_2, pz_2, om_2, fi_2, ka_2)
#camera_matrix_2 = matrix44FromRodrigues(px_2, py_2, pz_2, sx_2, sy_2, sz_2)
#camera_matrix_2 = matrix44FromQuaternion(px_2, py_2, pz_2, q0_2, q1_2, q2_2, q3_2)
R_2=camera_matrix_2[:-1,:-1]

lon1 = (u_1 / cols_1 - 0.5) * (2 * pi)
lat1 = -(v_1 / rows_1 - 0.5) * pi
bearing1_x=cos(lat1) * sin(lon1)
bearing1_y=-sin(lat1)
bearing1_z=cos(lat1) * cos(lon1)
bearing1 = Matrix([bearing1_x, bearing1_y, bearing1_z]).vec()

lon2 = (u_2 / cols_2 - 0.5) * (2 * pi)
lat2 = -(v_2 / rows_2 - 0.5) * pi
bearing2_x=cos(lat2) * sin(lon2)
bearing2_y=-sin(lat2)
bearing2_z=cos(lat2) * cos(lon2)
bearing2 = Matrix([bearing2_x, bearing2_y, bearing2_z]).vec()

bx=px_2-px_1
by=py_2-py_1
bz=pz_2-pz_1
b=Matrix([[0, -bz, by], [bz, 0, -bx], [-by, bx, 0]])

target_value = Matrix([[0]])
model_function = bearing1.transpose() * R_1t * b * R_2 * bearing2
obs_eq = target_value - model_function
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(obs_eq)
print(obs_eq_jacobian)

with open("equirectangular_camera_coplanarity_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void observation_equation_equirectangular_camera_coplanarity_tait_bryan_wc(double &delta, double rows_1, double cols_1, double pi, double u_1, double v_1, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double rows_2, double cols_2, double u_2, double v_2, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{")
    f_cpp.write("delta = %s;\n"%(ccode(obs_eq[0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_equirectangular_camera_coplanarity_tait_bryan_wc_jacobian(Eigen::Matrix<double, 1, 12, Eigen::RowMajor> &j, double rows_1, double cols_1, double pi, double u_1, double v_1, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double rows_2, double cols_2, double u_2, double v_2, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2)\n")
    f_cpp.write("{")
    for i in range (12):
        f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(0,i, ccode(obs_eq_jacobian[0,i])))
    f_cpp.write("}")






