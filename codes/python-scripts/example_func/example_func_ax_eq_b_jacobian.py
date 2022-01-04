from sympy import *

x, y = symbols('x y')
a, b = symbols('a b')
all_symbols = [a,b]

psi = a*x + b
target_value = 0
model_function = y - psi
obs_eq = Matrix([target_value - model_function]).vec()
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(model_function)
print(obs_eq)
print(obs_eq_jacobian)
print(latex(model_function))
print(latex(obs_eq_jacobian))

with open("example_func_ax_eq_b_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void example_func_ax_eq_b(double &psi, double x, double a, double b)\n")
    f_cpp.write("{")
    f_cpp.write("psi = %s;\n"%(ccode(psi)))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_example_func_ax_eq_b(double &delta, double x, double y, double a, double b)\n")
    f_cpp.write("{")
    f_cpp.write("delta = %s;\n"%(ccode(obs_eq[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_example_func_ax_eq_b_jacobian(Eigen::Matrix<double, 1, 2> &j, double x, double y, double a, double b)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (2):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_jacobian[i,j])))
    f_cpp.write("}")
  


  

