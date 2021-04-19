from sympy import *
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

a_1, b_1, c_1, d_1 = symbols('a_1 b_1 c_1 d_1')
px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
om_1, fi_1, ka_1 = symbols('om_1 fi_1 ka_1')
#sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
#q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
a_2, b_2, c_2, d_2 = symbols('a_2 b_2 c_2 d_2')
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

RT_wc_1 = matrix44FromTaitBryan(px_1, py_1, pz_1, om_1, fi_1, ka_1)
#RT_wc_1 = matrix44FromRodrigues(px_1, py_1, pz_1, sx_1, sy_1, sz_1)
#RT_wc_1 = matrix44FromQuaternion(px_1, py_1, pz_1, q0_1, q1_1, q2_1, q3_1)
RT_wc_2 = matrix44FromTaitBryan(px_2, py_2, pz_2, om_2, fi_2, ka_2)
#RT_wc_2 = matrix44FromRodrigues(px_2, py_2, pz_2, sx_2, sy_2, sz_2)
#RT_wc_2 = matrix44FromQuaternion(px_2, py_2, pz_2, q0_2, q1_2, q2_2, q3_2)

R_cw_1=RT_wc_1[:-1,:-1].transpose()
T_wc_1=Matrix([px_1, py_1, pz_1]).vec()
T_cw_1=-R_cw_1*T_wc_1
RT_cw_1=Matrix.hstack(R_cw_1, T_cw_1)
RT_cw_1=Matrix.vstack(RT_cw_1, Matrix([[0,0,0,1]]))

R_cw_2=RT_wc_2[:-1,:-1].transpose()
T_wc_2=Matrix([px_2, py_2, pz_2]).vec()
T_cw_2=-R_cw_2*T_wc_2
RT_cw_2=Matrix.hstack(R_cw_2, T_cw_2)
RT_cw_2=Matrix.vstack(RT_cw_2, Matrix([[0,0,0,1]]))

plane_1 = Matrix([[a_1, b_1, c_1, d_1]])
plane_2 = Matrix([[a_2, b_2, c_2, d_2]])

target_value = Matrix([0,0,0,0]).vec().transpose()
model_function = plane_1 * RT_cw_1 - plane_2 * RT_cw_2
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)


with open("plane_to_plane_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void plane_to_plane_tait_bryan_wc(Eigen::Matrix<double, 4, 1> &delta, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double a_1, double b_1, double c_1, double d_1, double a_2, double b_2, double c_2, double d_2)\n")
    f_cpp.write("{")
    for i in range (4):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void plane_to_plane_tait_bryan_wc_jacobian(Eigen::Matrix<double, 4, 12, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double a_1, double b_1, double c_1, double d_1, double a_2, double b_2, double c_2, double d_2)\n")
    f_cpp.write("{")
    for i in range (4):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")



