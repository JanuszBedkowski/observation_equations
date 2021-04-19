from sympy import *
from tait_bryan_R_utils import *

x, y, z = symbols('x y z')
px, py, pz = symbols('px py pz')
om, fi, ka = symbols('om fi ka')

point = Matrix([x, y, z, 1]).vec()
position_symbols = [px, py, pz]
orientation_symbols = [om, fi, ka]
all_symbols = position_symbols + orientation_symbols

transformation_matrix=matrix44FromTaitBryan(px, py, pz, om, fi, ka)
transformed_point = (transformation_matrix * point)[:-1,:]
transformed_point_jacobian = (transformed_point).jacobian(all_symbols)
transformed_point_jacobian.simplify()

with open("transform_point_jacobian.h",'w') as f_cpp:
    f_cpp.write("inline void transform_point3D_Tait_Bryan_jacobian(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double x, double y, double z, double px, double py, double pz, double om, double fi, double ka)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (6):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(transformed_point_jacobian[i,j])))
    f_cpp.write("}")





