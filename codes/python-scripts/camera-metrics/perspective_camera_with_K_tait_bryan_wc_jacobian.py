from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

fx,fy,cx,cy = symbols('fx fy cx cy');
tx, ty, tz = symbols('tx ty tz')
om, fi, ka = symbols('om fi ka')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
px, py, pz = symbols('px py pz'); 
u_kp, v_kp = symbols('u_kp v_kp')

position_symbols = [tx, ty, tz]
orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
point_symbols = [px, py, pz]
K_symbols = [fx,fy,cx,cy]
all_symbols = K_symbols + position_symbols + orientation_symbols + point_symbols

Rt_wc = matrix44FromTaitBryan(tx, ty, tz, om, fi, ka)
#Rt_wc = matrix44FromRodrigues(tx, ty, tz, sx, sy, sz)
#Rt_wc = matrix44FromQuaternion(tx, ty, tz, q0, q1, q2, q3)
R_cw=Rt_wc[:-1,:-1].transpose()
t_wc=Matrix([tx, ty, tz]).vec()
t_cw=-R_cw*t_wc

point_global=Matrix([px, py, pz]).vec()
point_local = R_cw * point_global + t_cw;
s = 1.0 / point_local[2];
u = s * fx * point_local[0] + cx
v = s * fy * point_local[1] + cy
u_delta = u_kp - u;
v_delta = v_kp - v;
uv = Matrix([u, v]).vec()

target_value = Matrix([u_kp, v_kp]).vec()
model_function = Matrix([u, v]).vec()
obs_eq = target_value - model_function
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(obs_eq)
print(obs_eq_jacobian)

with open("perspective_camera_with_K_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void observation_equation_perspective_camera_with_K_tait_bryan_wc(Eigen::Matrix<double, 2, 1> &delta, double fx, double fy, double cx, double cy, double tx, double ty, double tz, double om, double fi, double ka, double px, double py, double pz, double u_kp, double v_kp)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(obs_eq[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(obs_eq[1,0])))
    f_cpp.write("}")
    f_cpp.write("inline void projection_perspective_camera_with_K_tait_bryan_wc(double &u, double &v, double fx, double fy, double cx, double cy, double tx, double ty, double tz, double om, double fi, double ka, double px, double py, double pz)\n")
    f_cpp.write("{")
    f_cpp.write("u = %s;\n"%(ccode(uv[0,0])))
    f_cpp.write("v = %s;\n"%(ccode(uv[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_perspective_camera_with_K_tait_bryan_wc_jacobian(Eigen::Matrix<double, 2, 13> &j, double fx, double fy, double cx, double cy, double tx, double ty, double tz, double om, double fi, double ka, double px, double py, double pz)\n")
    f_cpp.write("{")
    for i in range (2):
        for j in range (13):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_jacobian[i,j])))
    f_cpp.write("}")
    
  

