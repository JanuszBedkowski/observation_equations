from sympy import *

px_1, py_1, pz_1 = symbols('px_1 py_1 pz_1')
r11_1, r12_1, r13_1, r21_1, r22_1, r23_1, r31_1, r32_1, r33_1 = symbols('r11_1 r12_1 r13_1 r21_1 r22_1 r23_1 r31_1 r32_1 r33_1')

px_2, py_2, pz_2 = symbols('px_2 py_2 pz_2')
r11_2, r12_2, r13_2, r21_2, r22_2, r23_2, r31_2, r32_2, r33_2 = symbols('r11_2 r12_2 r13_2 r21_2 r22_2 r23_2 r31_2 r32_2 r33_2')

px_m, py_m, pz_m = symbols('px_m py_m pz_m')
r11_m, r12_m, r13_m, r21_m, r22_m, r23_m, r31_m, r32_m, r33_m = symbols('r11_m r12_m r13_m r21_m r22_m r23_m r31_m r32_m r33_m')

position_symbols_1 = [px_1, py_1, pz_1]
orientation_symbols_1 = [r11_1, r12_1, r13_1, r21_1, r22_1, r23_1, r31_1, r32_1, r33_1]

position_symbols_2 = [px_2, py_2, pz_2]
orientation_symbols_2 = [r11_2, r12_2, r13_2, r21_2, r22_2, r23_2, r31_2, r32_2, r33_2]

all_symbols = position_symbols_1 + orientation_symbols_1 + position_symbols_2 + orientation_symbols_2

RT_wc_1 = Matrix([[r11_1, r12_1, r13_1, px_1],[r21_1, r22_1, r23_1, py_1],[r31_1, r32_1, r33_1, pz_1],[0,0,0,1]])
RT_wc_2 = Matrix([[r11_2, r12_2, r13_2, px_2],[r21_2, r22_2, r23_2, py_2],[r31_2, r32_2, r33_2, pz_2],[0,0,0,1]])



R_cw_1=RT_wc_1[:-1,:-1].transpose()
T_wc_1=Matrix([px_1, py_1, pz_1]).vec()
T_cw_1=-R_cw_1*T_wc_1
RT_cw_1=Matrix.hstack(R_cw_1, T_cw_1)
RT_cw_1=Matrix.vstack(RT_cw_1, Matrix([[0,0,0,1]]))

rp = RT_cw_1 * RT_wc_2

target_value = Matrix([px_m, py_m, pz_m, r11_m, r12_m, r13_m, r21_m, r22_m, r23_m, r31_m, r32_m, r33_m])


model_function = Matrix([rp[0,3], rp[1,3], rp[2,3], rp[0,0], rp[0,1], rp[0,2], rp[1,0], rp[1,1], rp[1,2], rp[2,0], rp[2,1], rp[2,2]])
delta=target_value-model_function
delta_jacobian=delta.jacobian(all_symbols)

print(delta)
print(delta_jacobian)


with open("relative_pose_wc_jacobian.h",'w') as f_cpp:  
    f_cpp.write("inline void relative_pose_obs_eq_wc(Eigen::Matrix<double, 12, 1> &delta, double px_1, double py_1, double pz_1, double r11_1, double r12_1, double r13_1, double r21_1, double r22_1, double r23_1, double r31_1, double r32_1, double r33_1, double px_2, double py_2, double pz_2, double r11_2, double r12_2, double r13_2, double r21_2, double r22_2, double r23_2, double r31_2, double r32_2, double r33_2, double px_m, double py_m, double pz_m, double r11_m, double r12_m, double r13_m, double r21_m, double r22_m, double r23_m, double r31_m, double r32_m, double r33_m)\n")
    f_cpp.write("{")
    for i in range (12):
        f_cpp.write("delta.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(delta[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
    f_cpp.write("inline void relative_pose_obs_eq_wc_jacobian(Eigen::Matrix<double, 12, 24, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double r11_1, double r12_1, double r13_1, double r21_1, double r22_1, double r23_1, double r31_1, double r32_1, double r33_1, double px_2, double py_2, double pz_2, double r11_2, double r12_2, double r13_2, double r21_2, double r22_2, double r23_2, double r31_2, double r32_2, double r33_2)\n")
    f_cpp.write("{")
    for i in range (12):
        for j in range (24):
            f_cpp.write("j.coeffRef(%d,%d) = %s;\n"%(i,j, ccode(delta_jacobian[i,j])))
    f_cpp.write("}")
    f_cpp.write("inline void relative_pose_wc(Eigen::Matrix<double, 12, 1> &relative_pose, double px_1, double py_1, double pz_1, double r11_1, double r12_1, double r13_1, double r21_1, double r22_1, double r23_1, double r31_1, double r32_1, double r33_1, double px_2, double py_2, double pz_2, double r11_2, double r12_2, double r13_2, double r21_2, double r22_2, double r23_2, double r31_2, double r32_2, double r33_2)\n")
    f_cpp.write("{")
    for i in range (12):
        f_cpp.write("relative_pose.coeffRef(%d,%d) = %s;\n"%(i, 0, ccode(model_function[i])))
    f_cpp.write("}")
    f_cpp.write("\n")
