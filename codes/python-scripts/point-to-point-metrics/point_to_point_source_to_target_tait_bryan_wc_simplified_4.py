import sympy
from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_t, y_t, z_t = symbols('x_t y_t z_t')
x_s, y_s, z_s = symbols('x_s y_s z_s')
tx, ty, tz = symbols('tx ty tz')
om, fi, ka = symbols('om fi ka')

p11, p12, p13, p21, p22, p23, p31, p32, p33 = symbols('p11 p12 p13 p21 p22 p23 p31 p32 p33')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

position_symbols = [tx, ty, tz]
orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + orientation_symbols

point_t = Matrix([x_t, y_t, z_t]).vec()
point_s = Matrix([x_s, y_s, z_s, 1]).vec()
transformed_point_s = (matrix44FromTaitBryan(tx, ty, tz, om, fi, ka) * point_s)[:-1,:]
#transformed_point_s = (matrix44FromRodrigues(tx, ty, tz, sx, sy, sz) * point_s)[:-1,:]
#transformed_point_s = (matrix44FromQuaternion(tx, ty, tz, q0, q1, q2, q3) * point_s)[:-1,:]

target_value = Matrix([0,0,0]).vec()
model_function = transformed_point_s-point_t
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

#P=Matrix([[p11, p12, p13],[p21, p22, p23],[p31, p32, p33]])
#AtPA = delta_jacobian.transpose() * P * delta_jacobian
#AtPB = delta_jacobian.transpose() * P * delta

AtPA = delta_jacobian.transpose() * delta_jacobian
AtPB = delta_jacobian.transpose() * delta

#print(delta)
#print(delta_jacobian)
#print(ATPA)

sin_om, cos_om = symbols('sin_om cos_om')
sin_fi, cos_fi = symbols('sin_fi cos_fi')
sin_ka, cos_ka = symbols('sin_ka cos_ka')
#sin_om, cos_om = symbols('0 1')
#sin_fi, cos_fi = symbols('0 1')
#sin_ka, cos_ka = symbols('0 1')

substitutions = [(sympy.sin(om), sin_om), (sympy.cos(om), cos_om), (sympy.sin(
fi), sin_fi), (sympy.cos(fi), cos_fi), (sympy.sin(ka), sin_ka), (sympy.cos(ka), cos_ka)]

delta = sympy.simplify(delta.subs(substitutions))
delta_variables, delta_simple = sympy.cse(
        delta, order='none')
delta_simple = delta_simple[0]

delta_jacobian = sympy.simplify(delta_jacobian.subs(substitutions))
delta_jacobian_variables, delta_jacobian_simple = sympy.cse(
        delta_jacobian, order='none')
delta_jacobian_simple = delta_jacobian_simple[0]

AtPA = sympy.simplify(AtPA.subs(substitutions))
AtPA_variables, AtPA_simple = sympy.cse(
        AtPA, order='none')
AtPA_simple = AtPA_simple[0]


AtPB = sympy.simplify(AtPB.subs(substitutions))
AtPB_variables, AtPB_simple = sympy.cse(
        AtPB, order='none')
AtPB_simple = AtPB_simple[0]


#print("-----------")
#print(delta)
#print(delta_jacobian_simple)
#print(AtPA)
#print("AtPB")
#print(AtPB)


with open("point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_4.h",'w') as f_cpp:  
    f_cpp.write("#ifndef _point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_4_h_\n")
    f_cpp.write("#define _point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_4_h_\n")
    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_simplified(double &delta_x, double &delta_y, double &delta_z, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{\n")
    f_cpp.write("double sin_om = sin(om);\n")
    f_cpp.write("double cos_om = cos(om);\n")
    f_cpp.write("double sin_fi = sin(fi);\n")
    f_cpp.write("double cos_fi = cos(fi);\n")
    f_cpp.write("double sin_ka = sin(ka);\n")
    f_cpp.write("double cos_ka = cos(ka);\n")
    for name, value_expr in delta_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    f_cpp.write("delta_x = %s;\n"%(ccode(delta_simple[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta_simple[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta_simple[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_4(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s)\n")
    f_cpp.write("{\n")
    f_cpp.write("double sin_om = sin(om);\n")
    f_cpp.write("double cos_om = cos(om);\n")
    f_cpp.write("double sin_fi = sin(fi);\n")
    f_cpp.write("double cos_fi = cos(fi);\n")
    f_cpp.write("double sin_ka = sin(ka);\n")
    f_cpp.write("double cos_ka = cos(ka);\n")
    for name, value_expr in delta_jacobian_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    for i in range (3):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian_simple[i,j])))
    f_cpp.write("}\n")
    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_AtPA_simplified_4(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &AtPA, const double &tx, const double &ty, const double &tz, const double &om, const double &fi, const double &ka, const double &x_s, const double &y_s, const double &z_s)\n")
    f_cpp.write("{\n")
    f_cpp.write("double sin_om = sin(om);\n")
    f_cpp.write("double cos_om = cos(om);\n")
    f_cpp.write("double sin_fi = sin(fi);\n")
    f_cpp.write("double cos_fi = cos(fi);\n")
    f_cpp.write("double sin_ka = sin(ka);\n")
    f_cpp.write("double cos_ka = cos(ka);\n")
    for name, value_expr in AtPA_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    for i in range (6):
        for j in range (6):
            f_cpp.write("AtPA.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(AtPA_simple[i,j])))
    f_cpp.write("}\n")
    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_AtPB_simplified_4(Eigen::Matrix<double, 6, 1> &AtPB, const double &tx, const double &ty, const double &tz, const double &om, const double &fi, const double &ka, const double &x_s, const double &y_s, const double &z_s, const double &x_t, const double &y_t, const double &z_t)\n")
    f_cpp.write("{\n")
    f_cpp.write("double sin_om = sin(om);\n")
    f_cpp.write("double cos_om = cos(om);\n")
    f_cpp.write("double sin_fi = sin(fi);\n")
    f_cpp.write("double cos_fi = cos(fi);\n")
    f_cpp.write("double sin_ka = sin(ka);\n")
    f_cpp.write("double cos_ka = cos(ka);\n")
    for name, value_expr in AtPB_variables:
        f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    for i in range (6):
        f_cpp.write("AtPB.coeffRef(%d) = %s;\n"%(i, ccode(AtPB_simple[i])))
    f_cpp.write("}\n")
    f_cpp.write("#endif\n")
#    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_d2sum_dbeta2(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
#    f_cpp.write("{")
#    for i in range (6):
#        for j in range (6):
#            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbeta2[i,j])))
#    f_cpp.write("}")
#    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_d2sum_dbetadx(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
#    f_cpp.write("{")
#    for i in range (6):
#        for j in range (6):
#            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbetadx[i,j])))
#    f_cpp.write("}")



#beta_symbols = position_symbols + orientation_symbols
#x_symbols = [x_s, y_s, z_s, x_t, y_t, z_t]
#sum=Matrix([delta[0,0]*delta[0,0]+delta[1,0]*delta[1,0]+delta[2,0]*delta[2,0]]).vec()
#d2sum_dbeta2=sum.jacobian(beta_symbols).jacobian(beta_symbols)
#d2sum_dbetadx=sum.jacobian(beta_symbols).jacobian(x_symbols)

#with open("point_to_point_source_to_target_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
#    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
#    f_cpp.write("{")
#    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
#    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
#    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
#    f_cpp.write("}")
#    f_cpp.write("\n")
#    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s)\n")
#    f_cpp.write("{")
#    for i in range (3):
#        for j in range (6):
#            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
#    f_cpp.write("}")
#    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_d2sum_dbeta2(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
#    f_cpp.write("{")
#    for i in range (6):
#        for j in range (6):
#            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbeta2[i,j])))
#    f_cpp.write("}")
#    f_cpp.write("inline void point_to_point_source_to_target_tait_bryan_wc_d2sum_dbetadx(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
#    f_cpp.write("{")
#    for i in range (6):
#        for j in range (6):
#            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(d2sum_dbetadx[i,j])))
#    f_cpp.write("}")


