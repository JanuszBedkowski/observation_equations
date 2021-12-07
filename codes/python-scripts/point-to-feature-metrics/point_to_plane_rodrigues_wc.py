from sympy import *
import sys
sys.path.insert(1, '..')
from rodrigues_R_utils import *

xsl, ysl, zsl = symbols('xsl ysl zsl')
xtg, ytg, ztg = symbols('xtg ytg ztg')
px, py, pz = symbols('px py pz')
sx, sy, sz = symbols('sx sy sz')
vzx, vzy, vzz = symbols('vzx vzy vzz')

position_symbols = [px, py, pz]
rodrigues_symbols = [sx, sy, sz]
all_symbols = position_symbols + rodrigues_symbols

point_source_local = Matrix([xsl, ysl, zsl,1]).vec()
RT_wc = matrix44FromRodrigues(px, py, pz, sx, sy, sz)[:-1,:]
point_source_global = RT_wc * point_source_local
point_on_line_target_global = Matrix([xtg, ytg, ztg])

vzt = Matrix([[vzx, vzy, vzz]]).transpose()

delta = Matrix([[0 - vzt.dot(point_source_global - point_on_line_target_global)]])
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_plane_rodrigues_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_plane_rodrigues_wc(Eigen::Matrix<double, 1, 1> &delta, double px, double py, double pz, double sx, double sy, double sz, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vzx, double vzy, double vzz)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_plane_rodrigues_wc_jacobian(Eigen::Matrix<double, 1, 6> &j, double px, double py, double pz, double sx, double sy, double sz, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vzx, double vzy, double vzz)\n")
    f_cpp.write("{")
    for i in range (1):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")

