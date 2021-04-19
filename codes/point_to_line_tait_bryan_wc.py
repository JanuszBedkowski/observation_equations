from sympy import *
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

xsl, ysl, zsl = symbols('xsl ysl zsl')
xtg, ytg, ztg = symbols('xtg ytg ztg')
px, py, pz = symbols('px py pz')
om, fi, ka = symbols('om fi ka')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
vxx, vxy, vxz = symbols('vxx vxy vxz')
vyx, vyy, vyz = symbols('vyx vyy vyz')

position_symbols = [px, py, pz]
orientation_symbols = [om, fi, ka]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + orientation_symbols

point_source_local = Matrix([xsl, ysl, zsl,1]).vec()
RT_wc = matrix44FromTaitBryan(px, py, pz, om, fi, ka)[:-1,:]
#RT_wc = matrix44FromRodrigues(px, py, pz, sx, sy, sz)[:-1,:]
#RT_wc = matrix44FromQuaternion(px, py, pz, q0, q1, q2, q3)[:-1,:]
point_source_global = RT_wc * point_source_local
point_on_line_target_global = Matrix([xtg, ytg, ztg])

vxt = Matrix([[vxx, vxy, vxz]]).transpose()
vyt = Matrix([[vyx, vyy, vyz]]).transpose()

target_value = Matrix([0,0]).vec()
model_function = Matrix([vxt.dot(point_source_global - point_on_line_target_global),vyt.dot(point_source_global - point_on_line_target_global)]).vec()
delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)

with open("point_to_line_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void point_to_line_tait_bryan_wc(Eigen::Matrix<double, 2, 1> &delta, double px, double py, double pz, double om, double fi, double ka, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)\n")
    f_cpp.write("{")
    f_cpp.write("delta.coeffRef(0,0) = %s;\n"%(ccode(delta[0,0])))
    f_cpp.write("delta.coeffRef(1,0) = %s;\n"%(ccode(delta[1,0])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_line_tait_bryan_wc_jacobian(Eigen::Matrix<double, 2, 6> &j, double px, double py, double pz, double om, double fi, double ka, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)\n")
    f_cpp.write("{")
    for i in range (2):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")

