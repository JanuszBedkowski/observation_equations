from sympy import *

x,y,z = symbols('x y z')
all_symbols = [x,y]
ro_x=5
ro_y=5

psi=exp(-(((x-1)*(x-1))/(2*ro_x*ro_x) + ((y-2)*(y-2))/(2*ro_y*ro_y)))

target_value = 0
model_function = z - psi
obs_eq = Matrix([target_value - model_function]).vec()
obs_eq_jacobian = obs_eq.jacobian(all_symbols)

print(psi)
print(obs_eq)
print(obs_eq_jacobian)
print(latex(psi))
print(latex(obs_eq_jacobian))

with open("example_func_xy_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void example_func_xy(double &psi, double x, double y)\n")
    f_cpp.write("{")
    f_cpp.write("psi = %s;\n"%(ccode(psi)))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_example_func_xy(double &delta, double x, double y, double z)\n")
    f_cpp.write("{")
    f_cpp.write("delta = %s;\n"%(ccode(obs_eq[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_example_func_xy_jacobian(Eigen::Matrix<double, 1, 2> &j, double x, double y)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (2):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_jacobian[i,j])))
    f_cpp.write("}")
  


  

