import sympy
from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

tx_tb, ty_tb, tz_tb = symbols('tx_tb ty_tb tz_tb')
om_tb, fi_tb, ka_tb = symbols('om_tb fi_tb ka_tb')
position_tb_symbols = [tx_tb, ty_tb, tz_tb]
orientation_tb_symbols = [om_tb, fi_tb, ka_tb]
all_symbols_tb = position_tb_symbols + orientation_tb_symbols
m_tb = matrix44FromTaitBryan(tx_tb, ty_tb, tz_tb, om_tb, fi_tb, ka_tb)
rodrigues = rodriguesFromMatrix44(m_tb)
rodrigues = Matrix([tx_tb, ty_tb, tz_tb, rodrigues[0], rodrigues[1], rodrigues[2]]).vec()

J_tb = rodrigues.jacobian(all_symbols_tb)
print(latex(J_tb))
print(shape(J_tb))
print("-------------------------------------------------")

tx_rodrigues, ty_rodrigues, tz_rodrigues = symbols('tx_rodrigues ty_rodrigues tz_rodrigues')
sx_rodrigues, sy_rodrigues, sz_rodrigues = symbols('sx_rodrigues sy_rodrigues sz_rodrigues')
position_rodrigues_symbols = [tx_rodrigues, ty_rodrigues, tz_rodrigues]
orientation_rodrigues_symbols = [sx_rodrigues, sy_rodrigues, sz_rodrigues]
all_symbols_rodrigues = position_rodrigues_symbols + orientation_rodrigues_symbols

m_rodrigues = matrix44FromRodrigues(tx_rodrigues, ty_rodrigues, tz_rodrigues, sx_rodrigues, sy_rodrigues, sz_rodrigues)
tb = taitBryanFromMatrix44Case1(m_rodrigues)
tb = Matrix([tx_rodrigues, ty_rodrigues, tz_rodrigues, tb[0], tb[1], tb[2]]).vec()
J_rodrigues = tb.jacobian(all_symbols_rodrigues)
print(latex(J_rodrigues))
print(shape(J_rodrigues))

with open("uncertainty_calculator.h",'w') as f_cpp:
    f_cpp.write("#ifndef __UNCERTAINTY_CALCULATOR_H__\n")
    f_cpp.write("#define __UNCERTAINTY_CALCULATOR_H__\n")
    f_cpp.write("inline void uncertainty_pose_tait_bryan_to_rodrigues(Eigen::Matrix<double, 6, 6> &j, double tx_tb, double ty_tb, double tz_tb, double om_tb, double fi_tb, double ka_tb)\n")
    f_cpp.write("{\n")
    for i in range (6):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(J_tb[i,j])))
    f_cpp.write("}\n")
    f_cpp.write("inline void uncertainty_pose_rodrigues_to_tait_bryan(Eigen::Matrix<double, 6, 6> &j, double tx_rodrigues, double ty_rodrigues, double tz_rodrigues, double sx_rodrigues, double sy_rodrigues, double sz_rodrigues)\n")
    f_cpp.write("{\n")
    for i in range (6):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(J_rodrigues[i,j])))
    f_cpp.write("}\n")
    f_cpp.write("#endif")


