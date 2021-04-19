from sympy import *

target_value, model_function = symbols('target_value model_function')
all_symbols = [model_function]

residual = Matrix([target_value - model_function]).vec()
residual_jacobian = residual.jacobian(all_symbols)

print(residual)
print(residual_jacobian)

with open("constraint_fixed_parameter_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void residual_constraint_fixed_optimization_parameter(double &residual, double target_value, double model_function)\n")
    f_cpp.write("{")
    f_cpp.write("residual = %s;\n"%(ccode(residual[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void residual_constraint_fixed_optimization_parameter_jacobian(Eigen::Matrix<double, 1, 1> &j, double target_value, double model_function)\n")
    f_cpp.write("{")
    f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(0,0, ccode(residual_jacobian[0,0])))
    f_cpp.write("}")
    

  

