from sympy import *
import sys
sys.path.insert(1, '..')
from tait_bryan_R_utils import *
from rodrigues_R_utils import *
from quaternion_R_utils import *

x_t, y_t, z_t = symbols('x_t y_t z_t')
#x_s, y_s, z_s = symbols('x_s y_s z_s')
tx, ty, tz = symbols('tx ty tz')
om, fi, ka = symbols('om fi ka')
ray_dir_x, ray_dir_y, ray_dir_z, ray_length = symbols('ray_dir_x ray_dir_y ray_dir_z ray_length')
plane_a, plane_b, plane_c, plane_d = symbols('plane_a plane_b plane_c plane_d')
#sx, sy, sz = symbols('sx sy sz')
#q0, q1, q2, q3 = symbols('q0 q1 q2 q3')

position_symbols = [tx, ty, tz]
orientation_symbols = [om, fi, ka]
plane_symbols = [plane_a, plane_b, plane_c, plane_d]
#orientation_symbols = [sx, sy, sz]
#orientation_symbols = [q0, q1, q2, q3]
all_symbols = position_symbols + orientation_symbols + plane_symbols

a = plane_a * ray_dir_x + plane_b * ray_dir_y + plane_c * ray_dir_z

intersection_x = - ray_dir_x * (plane_d/a)
intersection_y = - ray_dir_y * (plane_d/a)
intersection_z = - ray_dir_z * (plane_d/a)

n=Matrix([plane_a, plane_b, plane_c]).vec()
d=Matrix([ray_dir_x, ray_dir_y, ray_dir_z]).vec()
rd=2*d.dot(n)*n-d 

ll = ray_length - sqrt(intersection_x * intersection_x + intersection_y * intersection_y + intersection_z * intersection_z)
x_s=-(intersection_x + rd[0] * ll)
y_s=-(intersection_y + rd[1] * ll)
z_s=-(intersection_y + rd[2] * ll)

point_source = Matrix([x_s, y_s, z_s, 1]).vec()
point_target = Matrix([x_t, y_t, z_t]).vec()

transformed_point_source = (matrix44FromTaitBryan(tx, ty, tz, om, fi, ka) * point_source)[:-1,:]
target_value = Matrix([0,0,0]).vec()
model_function = transformed_point_source-point_target

delta = target_value - model_function
delta_jacobian=delta.jacobian(all_symbols)
print(delta)
print(delta_jacobian)
print(point_source)

with open("point_to_point_source_to_target_mirror_tait_bryan_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void transform_point_mirror_tait_bryan_wc(double &x, double &y, double &z, double tx, double ty, double tz, double om, double fi, double ka, double ray_dir_x, double ray_dir_y, double ray_dir_z, double ray_length, double plane_a, double plane_b, double plane_c, double plane_d)\n")
    f_cpp.write("{")
    f_cpp.write("x = %s;\n"%(ccode(transformed_point_source[0])))
    f_cpp.write("y = %s;\n"%(ccode(transformed_point_source[1])))
    f_cpp.write("z = %s;\n"%(ccode(transformed_point_source[2])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_target_mirror_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 10, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double ray_dir_x, double ray_dir_y, double ray_dir_z, double ray_length, double plane_a, double plane_b, double plane_c, double plane_d, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (10):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_target_mirror_tait_bryan_wc(Eigen::Matrix<double, 3, 1> delta, double tx, double ty, double tz, double om, double fi, double ka, double ray_dir_x, double ray_dir_y, double ray_dir_z, double ray_length, double plane_a, double plane_b, double plane_c, double plane_d, double x_t, double y_t, double z_t)\n")
    f_cpp.write("{")
    for i in range (3):
        for j in range (1):
            f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta[i,j])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void point_to_point_source_to_target_mirror_tait_bryan_wc_get_intersection(Eigen::Matrix<double, 3, 1> &intersection, double tx, double ty, double tz, double om, double fi, double ka, double ray_dir_x, double ray_dir_y, double ray_dir_z, double ray_length, double plane_a, double plane_b, double plane_c, double plane_d)\n")
    f_cpp.write("{")
    f_cpp.write("intersection.coeffRef(%d,%d) = %s;\n"%(0,0, ccode(intersection_x)))
    f_cpp.write("intersection.coeffRef(%d,%d) = %s;\n"%(1,0, ccode(intersection_y)))
    f_cpp.write("intersection.coeffRef(%d,%d) = %s;\n"%(2,0, ccode(intersection_z)))
    f_cpp.write("}")
    f_cpp.write("\n")



