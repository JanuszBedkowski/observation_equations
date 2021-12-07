from sympy import *
import sys
sys.path.insert(1, '..')
from rodrigues_R_utils import *

xsl, ysl, zsl = symbols('xsl ysl zsl')
xtg, ytg, ztg = symbols('xtg ytg ztg')
px, py, pz = symbols('px py pz')
sx, sy, sz = symbols('sx sy sz')
vxx, vxy, vxz = symbols('vxx vxy vxz')
vyx, vyy, vyz = symbols('vyx vyy vyz')

position_symbols = [px, py, pz]
rodrigues_symbols = [sx, sy, sz]
all_symbols = position_symbols + rodrigues_symbols

point_source_local = Matrix([xsl, ysl, zsl,1]).vec()
RT_wc = matrix44FromRodrigues(px, py, pz, sx, sy, sz)[:-1,:]
point_source_global = RT_wc * point_source_local
point_on_line_target_global = Matrix([xtg, ytg, ztg])

vxt = Matrix([[vxx, vxy, vxz]]).transpose()
vyt = Matrix([[vyx, vyy, vyz]]).transpose()

delta = Matrix([[0 - vxt.dot(point_source_global - point_on_line_target_global)],[0 - vyt.dot(point_source_global - point_on_line_target_global)]])
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_line_rodrigues_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_line_rodrigues_wc(Eigen::Matrix<double, 2, 1> &delta, double px, double py, double pz, double sx, double sy, double sz, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_line_rodrigues_wc_jacobian(Eigen::Matrix<double, 2, 6> &j, double px, double py, double pz, double sx, double sy, double sz, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)\n")
    f_cpp.write("{")
    for i in range (2):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")

