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
print("-------------------------------------------------")

tx_quaternion, ty_quaternion, tz_quaternion = symbols('tx_quaternion ty_quaternion tz_quaternion')
q0_quaternion, q1_quaternion, q2_quaternion, q3_quaternion = symbols('q0_quaternion q1_quaternion q2_quaternion q3_quaternion')
position_quaternion_symbols = [tx_quaternion, ty_quaternion, tz_quaternion]
orientation_quaternion_symbols = [q0_quaternion, q1_quaternion, q2_quaternion, q3_quaternion]
all_symbols_quaternion = position_quaternion_symbols + orientation_quaternion_symbols

#m_quaternion = matrix44FromQuaternion(tx_quaternion, ty_quaternion, tz_quaternion, q0_quaternion, q1_quaternion, q2_quaternion, q3_quaternion)
#tb = taitBryanFromMatrix44Case1(m_quaternion)
#tb = Matrix([tx_quaternion, ty_quaternion, tz_quaternion, tb[0], tb[1], tb[2]]).vec()
#J_quaternion = tb.jacobian(all_symbols_quaternion)
#q = Matrix([tx_quaternion, ty_quaternion, tz_quaternion, q0_quaternion, q1_quaternion, q2_quaternion, q3_quaternion]).vec()
#J_quaternion = q.jacobian(all_symbols_quaternion)

m_tb = matrix44FromTaitBryan(tx_tb, ty_tb, tz_tb, om_tb, fi_tb, ka_tb)
q = quaternionFromMatrix44(m_tb)
tb = Matrix([tx_tb, ty_tb, tz_tb, q[0], q[1], q[2], q[3]]).vec()

J_quaternion = tb.jacobian(all_symbols_tb)
#J_quaternion = tb.jacobian(all_symbols_quaternion)
print(latex(J_quaternion))
print(shape(J_quaternion))

m_quaternion = matrix44FromQuaternion(tx_quaternion, ty_quaternion, tz_quaternion, q0_quaternion, q1_quaternion, q2_quaternion, q3_quaternion)
tb_q = taitBryanFromMatrix44Case1(m_quaternion)
tb_q = Matrix([tx_quaternion, ty_quaternion, tz_quaternion, tb_q[0], tb_q[1], tb_q[2]]).vec()
J_quaternion_to_tb = tb_q.jacobian(all_symbols_quaternion)
print(latex(J_quaternion_to_tb))
print(shape(J_quaternion_to_tb))


with open("uncertainty_calculator.h",'w') as f_cpp:
    f_cpp.write("#ifndef __UNCERTAINTY_CALCULATOR_H__\n")
    f_cpp.write("#define __UNCERTAINTY_CALCULATOR_H__\n")
    f_cpp.write("inline void uncertainty_pose_rodrigues_to_tait_bryan(Eigen::Matrix<double, 6, 6> &j, double tx_tb, double ty_tb, double tz_tb, double om_tb, double fi_tb, double ka_tb)\n")
    f_cpp.write("{\n")
    for i in range (6):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(J_tb[i,j])))
    f_cpp.write("}\n")
    f_cpp.write("inline void uncertainty_pose_tait_bryan_to_rodrigues(Eigen::Matrix<double, 6, 6> &j, double tx_rodrigues, double ty_rodrigues, double tz_rodrigues, double sx_rodrigues, double sy_rodrigues, double sz_rodrigues)\n")
    f_cpp.write("{\n")
    for i in range (6):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(J_rodrigues[i,j])))
    f_cpp.write("}\n")
    f_cpp.write("inline void uncertainty_pose_quaternion_to_tait_bryan(Eigen::Matrix<double, 6, 7> &j, const double &tx_tb, const double &ty_tb, const double &tz_tb, const double &q0_quaternion, const double &q1_quaternion, const double &q2_quaternion, const double &q3_quaternion)\n")
    f_cpp.write("{\n")
    for i in range (6):
        for j in range (7):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(J_quaternion_to_tb[i,j])))
    f_cpp.write("}\n")
    f_cpp.write("inline void uncertainty_pose_tait_bryan_to_quaternion(Eigen::Matrix<double, 7, 6> &j, const double &tx_tb, const double &ty_tb, const double &tz_tb, const double &om_tb, const double &fi_tb, const double &ka_tb)\n")
    f_cpp.write("{\n")
    for i in range (7):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(J_quaternion[i,j])))
    f_cpp.write("}\n")
    f_cpp.write("#endif")


