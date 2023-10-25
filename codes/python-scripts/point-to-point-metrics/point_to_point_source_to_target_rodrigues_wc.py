import sympy
from sympy import *
import sys
sys.path.insert(1, '..')
from rodrigues_R_utils import *

x_t, y_t, z_t = symbols('x_t y_t z_t')
x_s, y_s, z_s = symbols('x_s y_s z_s')
px, py, pz = symbols('px py pz')
sx, sy, sz = symbols('sx sy sz')

position_symbols = [px, py, pz]
rodrigues_symbols = [sx, sy, sz]
all_symbols = position_symbols + rodrigues_symbols

point_t = Matrix([x_t, y_t, z_t]).vec()
point_s = Matrix([x_s, y_s, z_s, 1]).vec()
transformed_point_s = (matrix44FromRodrigues(px, py, pz, sx, sy, sz) * point_s)[:-1,:]

delta=Matrix([0,0,0]).vec()-(transformed_point_s-point_t)
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_point_source_to_target_rodrigues_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_point_source_to_target_rodrigues_wc(double &delta_x, double &delta_y, double &delta_z, double px, double py, double pz, double sx, double sy, double sz, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{")
    f_cpp.write("delta_x = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta_y = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("delta_z = %s;\n"%(ccode(delta[2,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_target_rodrigues_wc_jacobian(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double px, double py, double pz, double sx, double sy, double sz, double x_s, double y_s, double z_s)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")

print("---------------------------------------------")
substitutions = [(sx, 0.000000000001), (sy, 0.0), (sz, 0.0)]

delta = sympy.simplify(delta.subs(substitutions))
delta_variables, delta_simple = sympy.cse(
        delta, order='none')
delta_simple = delta_simple[0]

delta_jacobian = sympy.simplify(delta_jacobian.subs(substitutions))
delta_jacobian_variables, delta_jacobian_simple = sympy.cse(
        delta_jacobian, order='none')
delta_jacobian_simple = delta_jacobian_simple[0]

print("---------------------------------------------")
for name, value_expr in delta_jacobian_variables:
        print("%s,%s;\n"%(name,value_expr))

print(delta_jacobian_simple)
