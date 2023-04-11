from sympy import *
import sys

x,y,z = symbols('x y z')
r,polar_angle,azimuthal_angle = symbols('r alpha theta')

x = r * sin(polar_angle) * cos (azimuthal_angle)
y = r * sin(polar_angle) * sin (azimuthal_angle)
z = r * cos(polar_angle)

TLS_symbols = [r, polar_angle, azimuthal_angle]

coordinates = Matrix([x, y, z]).vec()

print(latex(coordinates))

jacobian=coordinates.jacobian(TLS_symbols)

print("----\n")
print(latex(jacobian))

with open("elementary_error_theory_for_terrestrial_laser_scanner_jacobian.h",'w') as f_cpp:
    f_cpp.write("#ifndef __ELEMENTARY_ERROR_THEORY_FOR_TERRESTRIAL_LASER_SCANNER_JACOBIAN_H__\n")
    f_cpp.write("#define __ELEMENTARY_ERROR_THEORY_FOR_TERRESTRIAL_LASER_SCANNER_JACOBIAN_H__\n")
    f_cpp.write("inline void elementary_error_theory_for_terrestrial_laser_scanner_jacobian(Eigen::Matrix<double, 3, 3> &j, double r, double alpha, double theta)\n")
    f_cpp.write("{\n")
    f_cpp.write("j(0,0) = %s;\n"%(ccode(jacobian[0,0])))
    f_cpp.write("j(0,1) = %s;\n"%(ccode(jacobian[0,1])))
    f_cpp.write("j(0,2) = %s;\n"%(ccode(jacobian[0,2])))
    f_cpp.write("j(1,0) = %s;\n"%(ccode(jacobian[1,0])))
    f_cpp.write("j(1,1) = %s;\n"%(ccode(jacobian[1,1])))
    f_cpp.write("j(1,2) = %s;\n"%(ccode(jacobian[1,2])))
    f_cpp.write("j(2,0) = %s;\n"%(ccode(jacobian[2,0])))
    f_cpp.write("j(2,1) = %s;\n"%(ccode(jacobian[2,1])))
    f_cpp.write("j(2,2) = %s;\n"%(ccode(jacobian[2,2])))
    f_cpp.write("}\n")
    f_cpp.write("#endif")