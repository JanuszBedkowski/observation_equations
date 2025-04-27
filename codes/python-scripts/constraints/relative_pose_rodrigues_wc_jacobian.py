import sympy
from sympy import *
import sys
sys.path.insert(1, '..')
from rodrigues_R_utils import *

px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
sx_1, sy_1, sz_1 = symbols('sx_1 sy_1 sz_1')
px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
sx_2, sy_2, sz_2 = symbols('sx_2 sy_2 sz_2')
px_m, py_m, pz_m = symbols('px_m py_m pz_m')
sx_m, sy_m, sz_m = symbols('sx_m sy_m sz_m')

position_symbols_1 = [px_1, py_1, pz_1]
rodrigues_symbols_1 = [sx_1, sy_1, sz_1]
position_symbols_2 = [px_2, py_2, pz_2]
rodrigues_symbols_2 = [sx_2, sy_2, sz_2]
all_symbols = position_symbols_1 + rodrigues_symbols_1 + position_symbols_2 + rodrigues_symbols_2

RT_wc_1 = matrix44FromRodrigues(px_1, py_1, pz_1, sx_1, sy_1, sz_1)
RT_wc_2 = matrix44FromRodrigues(px_2, py_2, pz_2, sx_2, sy_2, sz_2)

R_cw_1=RT_wc_1[:-1,:-1].transpose()
T_wc_1=Matrix([px_1, py_1, pz_1]).vec()
T_cw_1=-R_cw_1*T_wc_1
RT_cw_1=Matrix.hstack(R_cw_1, T_cw_1)
RT_cw_1=Matrix.vstack(RT_cw_1, Matrix([[0,0,0,1]]))

relative_pose = RT_cw_1 * RT_wc_2

rodrigues = rodriguesFromMatrix44(relative_pose)
t = relative_pose[0 : 3, 3]

observation = t.col_join(rodrigues)
measurement = Matrix([px_m, py_m, pz_m, sx_m, sy_m, sz_m])
delta = (measurement-observation)
delta_jacobian = delta.jacobian(all_symbols)

theta = symbols('theta')
sin_theta, cos_theta = symbols('sin_theta cos_theta')

#substitutions1 = [(sympy.sin(theta), sin_theta), (sympy.cos(theta), cos_theta)]

#print("computing delta1_simple")
#delta1 = sympy.simplify(delta.subs(substitutions1))
#delta1_variables, delta1_simple = sympy.cse(
#        delta1, order='none',symbols=symbols('a0:1000'))
#delta1_simple = delta1_simple[0]

#print("computing delta_jacobian1_simple")
#delta_jacobian1 = sympy.simplify(delta.subs(substitutions1))
#delta_jacobian1_variables, delta_jacobian1_simple = sympy.cse(
#        delta_jacobian1, order='none')
#delta_jacobian1_simple = delta_jacobian1_simple[0]

#print("computing model_function1_simple")
#model_function1 = sympy.simplify(delta.subs(substitutions1))
#model_function1_variables, model_function1_simple = sympy.cse(
#        model_function1, order='none')
#model_function1_simple = model_function1_simple[0]


#print("computing AtPA AtPB")
#p_x, p_y, p_z, p_sx, p_sy, p_sz = symbols('p_x p_y p_z p_sx p_sy p_sz')
#P=Matrix([[p_x, 0, 0, 0, 0, 0],[0, p_y, 0, 0, 0, 0],[0, 0, p_z, 0, 0, 0],[0, 0, 0, p_sx, 0, 0],[0, 0, 0, 0, p_sy, 0],[0, 0, 0, 0, 0, p_sz]])

#AtPA = delta_jacobian1_simple.transpose() * P * delta_jacobian1_simple
#AtPB = delta_jacobian1_simple.transpose() * P * delta1_simple

#AtPA = delta_jacobian.transpose() * P * delta_jacobian
#AtPB = delta_jacobian.transpose() * P * delta

#print("computing AtPA AtPB done")

with open("relative_pose_rodrigues_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("#ifndef _RELATIVE_POSE_RODRIGUES_WC_JACOBIAN_H_\n")
    f_cpp.write("#define _RELATIVE_POSE_RODRIGUES_WC_JACOBIAN_H_\n")
    f_cpp.write("inline void relative_pose_obs_eq_rodrigues_wc(Eigen::Matrix<double, 6, 1> &delta, double px_1, double py_1, double pz_1, double sx_1, double sy_1, double sz_1, double px_2, double py_2, double pz_2, double sx_2, double sy_2, double sz_2, double px_m, double py_m, double pz_m, double sx_m, double sy_m, double sz_m)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_obs_eq_rodrigues_wc_jacobian(Eigen::Matrix<double, 6, 12, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double sx_1, double sy_1, double sz_1, double px_2, double py_2, double pz_2, double sx_2, double sy_2, double sz_2)\n")
    f_cpp.write("{")
    for i in range (6):
        for j in range (12):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void relative_pose_rodrigues_wc(Eigen::Matrix<double, 6, 1> &relative_pose, double px_1, double py_1, double pz_1, double sx_1, double sy_1, double sz_1, double px_2, double py_2, double pz_2, double sx_2, double sy_2, double sz_2)\n")
    f_cpp.write("{")
    for i in range (6):
        f_cpp.write("relative_pose.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(observation[i])))
    f_cpp.write("}")
    ########################################### _simplified_ #######################################
    #f_cpp.write("inline void relative_pose_rodrigues_wc_AtPA_simplified(Eigen::Matrix<double, 12, 12> &AtPA, const double &px_1, const double &py_1, const double &pz_1, const double &sx_1, const double &sy_1, const double &sz_1, const double &px_2, const double &py_2, const double &pz_2, const double &sx_2, const double &sy_2, const double &sz_2, const double &p_x, const double &p_y, const double &p_z, const double &p_sx, const double &p_sy, const double &p_sz)\n")
    #f_cpp.write("{\n")
    #for name, value_expr in delta_jacobian1_variables:
    #    f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    #for i in range (12):
    #    for j in range (12):
    #        f_cpp.write("AtPA.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(AtPA[i,j])))
    #f_cpp.write("}\n")
    #f_cpp.write("inline void relative_pose_rodrigues_wc_AtPB_simplified(Eigen::Matrix<double, 12, 1> &AtPB, const double &px_1, const double &py_1, const double &pz_1, const double &sx_1, const double &sy_1, const double &sz_1, const double &px_2, const double &py_2, const double &pz_2, const double &sx_2, const double &sy_2, const double &sz_2, const double &px_m, const double &py_m, const double &pz_m, const double &sx_m, const double &sy_m, const double &sz_m, const double &p_x, const double &p_y, const double &p_z, const double &p_sx, const double &p_sy, const double &p_sz)\n")
    #f_cpp.write("{\n")
    #for name, value_expr in delta1_variables:
    #    f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    #for name, value_expr in delta_jacobian1_variables:
    #    f_cpp.write("double %s = %s;\n"%(name,ccode(value_expr)))
    #for i in range (12):
    #    f_cpp.write("AtPB.coeffRef(%d) = %s;\n"%(i, ccode(AtPB[i])))
    #f_cpp.write("}\n")
    #f_cpp.write("#endif\n")
    f_cpp.write("\n")

