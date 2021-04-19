from sympy import *

a,x,x_trg = symbols('a x x_trg')
all_symbols = [x]

obs_eq = Matrix([0 - (a*(x_trg - x))]).vec()
obs_eq_jacobian = obs_eq.jacobian(all_symbols)
obs_eq_sq = Matrix([0 - (a*(x_trg - x))*(a*(x_trg - x))]).vec()
obs_eq_sq_jacobian = obs_eq_sq.jacobian(all_symbols)

print(obs_eq)
print(obs_eq_jacobian)
print(obs_eq_sq)
print(obs_eq_sq_jacobian)

with open("constraints_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void observation_equation_constraint(double &delta, double a, double x, double x_trg)\n")
    f_cpp.write("{")
    f_cpp.write("delta = %s;\n"%(ccode(obs_eq[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_constraint_jacobian(Eigen::Matrix<double, 1, 1> &j, double a, double x, double x_trg)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (1):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void observation_equation_sq_constraint(double &delta, double a, double x, double x_trg)\n")
    f_cpp.write("{")
    f_cpp.write("delta = %s;\n"%(ccode(obs_eq_sq[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void observation_equation_sq_constraint_jacobian(Eigen::Matrix<double, 1, 1> &j, double a, double x, double x_trg)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (1):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(obs_eq_sq_jacobian[i,j])))
    f_cpp.write("}")
  


  

