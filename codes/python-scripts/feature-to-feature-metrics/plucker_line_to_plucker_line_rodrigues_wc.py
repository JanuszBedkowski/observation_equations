from sympy import *
import sys
sys.path.insert(1, '..')
from rodrigues_R_utils import *
from plucker_line_utils import *

mx_1, my_1, mz_1, lx_1, ly_1, lz_1  = symbols('mx_1 my_1 mz_1 lx_1 ly_1 lz_1')
px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
mx_2, my_2, mz_2, lx_2, ly_2, lz_2  = symbols('mx_2 my_2 mz_2 lx_2 ly_2 lz_2')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
sx_2, sy_2, sz_2 = symbols('sx_2 sy_2 sz_2')

position_symbols_1 = [px_1, py_1, pz_1]
rordigues_symbols_1 = [sx_1, sy_1, sz_1]
position_symbols_2 = [px_2, py_2, pz_2]
rordigues_symbols_2 = [sx_2, sy_2, sz_2]
all_symbols = position_symbols_1 + rordigues_symbols_1 + position_symbols_2 + rordigues_symbols_2

RT_wc_1 = matrix44FromRodrigues(px_1, py_1, pz_1, sx_1, sy_1, sz_1)
RT_wc_2 = matrix44FromRodrigues(px_2, py_2, pz_2, sx_2, sy_2, sz_2)
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

with open("plucker_line_to_plucker_line_rodrigues_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void plucker_line_to_plucker_line_rodrigues_wc(Eigen::Matrix<double, 6, 1> &delta, double px_1, double py_1, double pz_1, double sx_1, double sy_1, double sz_1, double px_2, double py_2, double pz_2, double sx_2, double sy_2, double sz_2, double mx_1, double my_1, double mz_1, double lx_1, double ly_1, double lz_1, double mx_2, double my_2, double mz_2, double lx_2, double ly_2, double lz_2)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void plucker_line_to_plucker_line_rodrigues_wc_jacobian(Eigen::Matrix<double, 6, 12, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double sx_1, double sy_1, double sz_1, double px_2, double py_2, double pz_2, double sx_2, double sy_2, double sz_2, double mx_1, double my_1, double mz_1, double lx_1, double ly_1, double lz_1, double mx_2, double my_2, double mz_2, double lx_2, double ly_2, double lz_2)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")



