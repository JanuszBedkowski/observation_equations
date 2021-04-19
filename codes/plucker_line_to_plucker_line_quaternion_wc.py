from sympy import *
from quaternion_R_utils import *
from plucker_line_utils import *

mx_1, my_1, mz_1, lx_1, ly_1, lz_1  = symbols('mx_1 my_1 mz_1 lx_1 ly_1 lz_1')
px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
mx_2, my_2, mz_2, lx_2, ly_2, lz_2  = symbols('mx_2 my_2 mz_2 lx_2 ly_2 lz_2')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')

position_symbols_1 = [px_1, py_1, pz_1]
quaternion_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
position_symbols_2 = [px_2, py_2, pz_2]
quaternion_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
all_symbols = position_symbols_1 + quaternion_symbols_1 + position_symbols_2 + quaternion_symbols_2

RT_wc_1 = matrix44FromQuaternion(px_1, py_1, pz_1, q0_1, q1_1, q2_1, q3_1)
RT_wc_2 = matrix44FromQuaternion(px_2, py_2, pz_2, q0_2, q1_2, q2_2, q3_2)
plucker_line_motion_matrix_1=plucker_line_motion_matrix_wc(RT_wc_1)
plucker_line_motion_matrix_2=plucker_line_motion_matrix_wc(RT_wc_2)

plucker_line_local_1 = Matrix([mx_1, my_1, mz_1, lx_1, ly_1, lz_1]).vec()
plucker_line_local_2 = Matrix([mx_2, my_2, mz_2, lx_2, ly_2, lz_2]).vec()

plucker_line_global_1 = plucker_line_motion_matrix_1 * plucker_line_local_1
plucker_line_global_2 = plucker_line_motion_matrix_2 * plucker_line_local_2

delta = Matrix([0,0,0,0,0,0]).vec() - (plucker_line_global_1 - plucker_line_global_2)
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("plucker_line_to_plucker_line_quaternion_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void plucker_line_to_plucker_line_quaternion_wc(Eigen::Matrix<double, 6, 1> &delta, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double mx_1, double my_1, double mz_1, double lx_1, double ly_1, double lz_1, double mx_2, double my_2, double mz_2, double lx_2, double ly_2, double lz_2)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void plucker_line_to_plucker_line_quaternion_wc_jacobian(Eigen::Matrix<double, 6, 14, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double mx_1, double my_1, double mz_1, double lx_1, double ly_1, double lz_1, double mx_2, double my_2, double mz_2, double lx_2, double ly_2, double lz_2)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (14):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")



