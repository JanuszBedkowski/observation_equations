from sympy import *
import sys
sys.path.insert(1, '..')
from quaternion_R_utils import *

px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
q0_1, q1_1, q2_1, q3_1 = symbols('q0_1 q1_1 q2_1 q3_1')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
q0_2, q1_2, q2_2, q3_2 = symbols('q0_2 q1_2 q2_2 q3_2')
ksi_1, eta_1, ksi_2, eta_2, ksi_01, eta_01, ksi_02, eta_02, c_1, c_2 = symbols('ksi_1 eta_1 ksi_2 eta_2 ksi_01 eta_01 ksi_02 eta_02 c_1 c_2');

position_symbols_1 = [px_1, py_1, pz_1]
quaternion_symbols_1 = [q0_1, q1_1, q2_1, q3_1]
position_symbols_2 = [px_2, py_2, pz_2]
quaternion_symbols_2 = [q0_2, q1_2, q2_2, q3_2]
c_symbols = [c_1, c_2]
all_symbols = position_symbols_1 + quaternion_symbols_1 + position_symbols_2 + quaternion_symbols_2

bx=px_2-px_1
by=py_2-py_1
bz=pz_2-pz_1
b=Matrix([[0, -bz, by], [bz, 0, -bx], [-by, bx, 0]])
C_1t=Matrix([[1, 0, -ksi_01], [0, 1, -eta_01], [0, 0, -c_1]]).transpose()
C_2=Matrix([[1, 0, -ksi_02], [0, 1, -eta_02], [0, 0, -c_2]])

camera_matrix_1 = matrix44FromQuaternion(px_1, py_1, pz_1, q0_1, q1_1, q2_1, q3_1)
R_1t=camera_matrix_1[:-1,:-1].transpose()
camera_matrix_2 = matrix44FromQuaternion(px_2, py_2, pz_2, q0_2, q1_2, q2_2, q3_2)
R_2=camera_matrix_2[:-1,:-1]
ksieta_1=Matrix([[ksi_1, eta_1, 1]])
ksieta_2t=Matrix([[ksi_2, eta_2, 1]]).transpose()

obs_eq = Matrix([[0]]) - ksieta_1 * C_1t * R_1t * b * R_2 * C_2 * ksieta_2t
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(obs_eq)
print(obs_eq_jacobian)

with open("metric_camera_coplanarity_quaternion_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void observation_equation_metric_camera_coplanarity_quaternion_wc(double &delta, double ksi_01, double eta_01, double c_1, double ksi_1, double eta_1, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double ksi_02, double eta_02, double c_2, double ksi_2, double eta_2, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2)\n")
    f_cpp.write("{")
    f_cpp.write("delta = %s;\n"%(ccode(obs_eq[0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_metric_camera_coplanarity_quaternion_wc_jacobian(Eigen::Matrix<double, 1, 14, Eigen::RowMajor> &j, double ksi_01, double eta_01, double c_1, double ksi_1, double eta_1, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double ksi_02, double eta_02, double c_2, double ksi_2, double eta_2, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2)\n")
    f_cpp.write("{")
    for i in range (12):
        f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(0,i, ccode(obs_eq_jacobian[0,i])))
    f_cpp.write("}")






