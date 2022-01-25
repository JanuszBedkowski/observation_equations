from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

tx, ty, tz = symbols('tx ty tz')
om, fi, ka = symbols('om fi ka')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
px, py, pz = symbols('px py pz'); 
cols, rows = symbols('cols rows');
pi = symbols('pi')
u_kp, v_kp = symbols('u_kp v_kp')

position_symbols = [tx, ty, tz]
orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
point_symbols = [px, py, pz]
all_symbols = position_symbols + orientation_symbols + point_symbols

RT_wc = matrix44FromTaitBryan(tx, ty, tz, om, fi, ka)
#RT_wc = matrix44FromRodrigues(tx, ty, tz, sx, sy, sz)
#RT_wc = matrix44FromQuaternion(tx, ty, tz, q0, q1, q2, q3)
r=RT_wc[:-1,:-1]
t=Matrix([tx, ty, tz]).vec()

pos_w=Matrix([px, py, pz]).vec()
bearing = r * pos_w + t;
norm = sqrt(bearing[0]*bearing[0] + bearing[1]*bearing[1] + bearing[2]*bearing[2])
bearing=bearing/norm
latitude=-asin(bearing[1])
longitude=atan2(bearing[0], bearing[2])

u=cols*(0.5 + longitude / (2.0 * pi))
v=rows*(0.5 - latitude/pi)
uv = Matrix([u, v]).vec()

target_value = Matrix([u_kp, v_kp]).vec()
model_function = Matrix([u, v]).vec()
obs_eq = target_value - model_function
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(obs_eq)
print(obs_eq_jacobian)

with open("equirectangular_camera_colinearity_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void equrectangular_camera_colinearity_tait_bryan_wc(double &u, double &v, double rows, double cols, double pi, double tx, double ty, double tz, double om, double fi, double ka, double px, double py, double pz)\n")
    f_cpp.write("{")
    f_cpp.write("u = %s;\n"%(ccode(uv[0,0])))
    f_cpp.write("v = %s;\n"%(ccode(uv[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_equrectangular_camera_colinearity_tait_bryan_wc(Eigen::Matrix<double, 2, 1> &delta, double rows, double cols, double pi, double tx, double ty, double tz, double om, double fi, double ka, double px, double py, double pz, double u_kp, double v_kp)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(obs_eq[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(obs_eq[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_equrectangular_camera_colinearity_tait_bryan_wc_jacobian(Eigen::Matrix<double, 2, 9, Eigen::RowMajor> &j, double rows, double cols, double pi, double tx, double ty, double tz, double om, double fi, double ka, double px, double py, double pz, double u_kp, double v_kp)\n")
    f_cpp.write("{")
    for i in range (2):
        for j in range (9):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_jacobian[i,j])))
    f_cpp.write("}")






