from sympy import *

T_x_int, T_y_int, T_z_int = symbols('T_x_int T_y_int T_z_int')
T_x_ln, T_y_ln, T_z_ln = symbols('T_x_ln T_y_ln T_z_ln')
vxx, vxy, vxz = symbols('vxx vxy vxz')
vyx, vyy, vyz = symbols('vyx vyy vyz')

ray_intersection_symbols = [T_x_int, T_y_int, T_z_int]
all_symbols = ray_intersection_symbols

T_int = Matrix([T_x_int, T_y_int, T_z_int])
T_ln = Matrix([T_x_ln, T_y_ln, T_z_ln])
vxt = Matrix([[vxx, vxy, vxz]]).transpose()
vyt = Matrix([[vyx, vyy, vyz]]).transpose()

target_value = Matrix([0,0]).vec()
model_function = Matrix([vxt.dot(T_int - T_ln),vyt.dot(T_int - T_ln)]).vec()
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("ray_intersection_observation_equation_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void ray_intersection_observation_equation(Eigen::Matrix<double, 2, 1> &delta, double T_x_int, double T_y_int, double T_z_int, double T_x_ln, double T_y_ln, double T_z_ln, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void ray_intersection_observation_equation_jacobian(Eigen::Matrix<double, 2, 3> &j, double T_x_int, double T_y_int, double T_z_int, double T_x_ln, double T_y_ln, double T_z_ln, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)\n")
    f_cpp.write("{")
    for i in range (2):
        for j in range (3):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")

