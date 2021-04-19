from sympy import *

q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
quaternion_symbols = [q0, q1, q2, q3]
target_value = 1
model_function = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3)
delta = Matrix([target_value - model_function]).vec()
delta_jacobian = delta.jacobian(quaternion_symbols)

print(delta)
print(delta_jacobian)

with open("quaternion_constraint_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void quaternion_constraint(double &delta, double q0, double q1, double q2, double q3)\n")
    f_cpp.write("{")
    f_cpp.write("delta = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void quaternion_constraint_jacobian(Eigen::Matrix<double, 1, 4> &j, double q0, double q1, double q2, double q3)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (4):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
  


  

