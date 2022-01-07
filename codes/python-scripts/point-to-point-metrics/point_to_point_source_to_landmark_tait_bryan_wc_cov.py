from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_L, y_L, z_L = symbols('x_L y_L z_L')
x_s, y_s, z_s = symbols('x_s y_s z_s')
px, py, pz = symbols('px py pz')
om, fi, ka = symbols('om fi ka')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

position_symbols = [px, py, pz]
orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
landmark_symbols = [x_L, y_L, z_L]
source_point_symbols = [x_s, y_s, z_s]

beta_symbols = position_symbols + orientation_symbols
x_symbols = source_point_symbols + landmark_symbols

point_Landmark = Matrix([x_L, y_L, z_L]).vec()
point_source = Matrix([x_s, y_s, z_s, 1]).vec()
transformed_point_source = (matrix44FromTaitBryan(px, py, pz, om, fi, ka) * point_source)[:-1,:]
#transformed_point_source = (matrix44FromRodrigues(px, py, pz, sx, sy, sz) * point_source)[:-1,:]
#transformed_point_source = (matrix44FromQuaternion(px, py, pz, q0, q1, q2, q3) * point_source)[:-1,:]

target_value = Matrix([0,0,0]).vec()
model_function = transformed_point_source-point_Landmark
delta = target_value - model_function
sum=Matrix([delta[0,0]*delta[0,0]+delta[1,0]*delta[1,0]+delta[2,0]*delta[2,0]]).vec()
d2sum_dbeta2=sum.jacobian(beta_symbols).jacobian(beta_symbols).transpose()
d2sum_dxdbeta=sum.jacobian(x_symbols).jacobian(beta_symbols).transpose()

with open("point_to_point_source_to_landmark_tait_bryan_wc_cov.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_source_to_landmark_tait_bryan_wc_d2sum_dbeta2(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &d2sum_dbeta2, double px, double py, double pz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_L, double y_L, double z_L)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (6):
            f_cpp.write("d2sum_dbeta2.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbeta2[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_landmark_tait_bryan_wc_d2sum_dxdbeta(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &d2sum_dxdbeta, double px, double py, double pz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_L, double y_L, double z_L)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (6):
            f_cpp.write("d2sum_dxdbeta.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dxdbeta[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")

