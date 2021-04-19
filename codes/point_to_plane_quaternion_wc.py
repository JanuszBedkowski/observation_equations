from sympy import *
from quaternion_R_utils import *

xsl, ysl, zsl = symbols('xsl ysl zsl')
xtg, ytg, ztg = symbols('xtg ytg ztg')
px, py, pz = symbols('px py pz')
q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
vzx, vzy, vzz = symbols('vzx vzy vzz')

position_symbols = [px, py, pz]
quaternion_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + quaternion_symbols

point_source_local = Matrix([xsl, ysl, zsl,1]).vec()
RT_wc = matrix44FromQuaternion(px, py, pz, q0, q1, q2, q3)[:-1,:]
point_source_global = RT_wc * point_source_local
point_on_line_target_global = Matrix([xtg, ytg, ztg])

vzt = Matrix([[vzx, vzy, vzz]]).transpose()

delta = Matrix([[0 - vzt.dot(point_source_global - point_on_line_target_global)]])
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_plane_quaternion_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_plane_quaternion_wc(Eigen::Matrix<double, 1, 1> &delta, double px, double py, double pz, double q0, double q1, double q2, double q3, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vzx, double vzy, double vzz)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_plane_quaternion_wc_jacobian(Eigen::Matrix<double, 1, 7> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vzx, double vzy, double vzz)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (7):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")

