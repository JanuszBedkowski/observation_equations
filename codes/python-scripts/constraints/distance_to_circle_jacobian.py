from sympy import *

x, y = symbols('x y')
cx, cy, cr = symbols('cx cy cr')

circle_symbols = [cx, cy, cr]
all_symbols = circle_symbols

v = Matrix([x, y]).vec()
vc = Matrix([cx, cy]).vec()

diff = v-vc

target_value = 0
model_function = sqrt(diff[0] * diff[0] + diff[1] * diff[1] ) - cr

obs_eq = Matrix([[target_value - model_function]])
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print (obs_eq)
print (obs_eq_jacobian)

with open("distance_to_circle_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void observation_equation_distance_to_circle(double &delta, double x, double y, double cx, double cy, double cr)\n")
    f_cpp.write("{")
    f_cpp.write("delta = %s;\n"%(ccode(obs_eq[0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_distance_to_circle_jacobian(Eigen::Matrix<double, 1, 3, Eigen::RowMajor> &j, double x, double y, double cx, double cy, double cr)\n")
    f_cpp.write("{")
    for i in range (3):
        f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(0,i, ccode(obs_eq_jacobian[0,i])))
    f_cpp.write("}")








