from sympy import *
from sympy.physics.vector.functions import *

a_x,a_y,a_z,b_x,b_y,b_z = symbols('a_x a_y a_z b_x b_y b_z')
a,b,c,d = symbols('a b c d')
line_symbols=[a_x, a_y, a_z, b_x, b_y, b_z]
all_symbols = line_symbols

l = Matrix([b_x-a_x, b_y-a_y, b_z-a_z]).vec()
p = Matrix([a_x, a_y, a_z]).vec()
n = Matrix([a, b, c]).vec()

target_value_norm = 1
model_function_norm = sqrt(l[0]*l[0] + l[1]*l[1] + l[2]*l[2])
residual_norm = Matrix([target_value_norm - model_function_norm]).vec()
residual_norm_jacobian = residual_norm.jacobian(all_symbols)

target_value = Matrix([0,0,0,0]).vec()
model_function = Matrix([a * a_x + b * a_y + c * a_z + d, a * b_x + b * b_y + c * b_z + d, p.dot(l), l.dot(n)]).vec()
residual = Matrix([target_value - model_function]).vec()
residual_jacobian = residual.jacobian(all_symbols)

with open("sheaf_of_planes_observation_equation_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void sheaf_of_planes_observation_equation(Eigen::Matrix<double, 4, 1> &residual, double a_x, double a_y, double a_z, double b_x, double b_y, double b_z, double a, double b, double c, double d)\n")
    f_cpp.write("{")
    f_cpp.write("residual.coeffRef(0,0) = %s;\n"%(ccode(residual[0])))
    f_cpp.write("residual.coeffRef(1,0) = %s;\n"%(ccode(residual[1])))
    f_cpp.write("residual.coeffRef(2,0) = %s;\n"%(ccode(residual[2])))
    f_cpp.write("residual.coeffRef(3,0) = %s;\n"%(ccode(residual[3])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void sheaf_of_planes_observation_equation_jacobian(Eigen::Matrix<double, 4, 6> &j, double a_x, double a_y, double a_z, double b_x, double b_y, double b_z, double a, double b, double c, double d)\n")
    f_cpp.write("{")
    for i in range (4):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(residual_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void line_direction_norm_observation_equation(Eigen::Matrix<double, 1, 1> &residual, double a_x, double a_y, double a_z, double b_x, double b_y, double b_z)\n")
    f_cpp.write("{")
    f_cpp.write("residual.coeffRef(0,0) = %s;\n"%(ccode(residual_norm[0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void line_direction_observation_equation_jacobian(Eigen::Matrix<double, 1, 6> &j, double a_x, double a_y, double a_z, double b_x, double b_y, double b_z)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(residual_norm_jacobian[i,j])))
    f_cpp.write("}")

