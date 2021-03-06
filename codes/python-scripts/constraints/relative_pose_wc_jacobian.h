inline void relative_pose_obs_eq_wc(Eigen::Matrix<double, 12, 1> &delta, double px_1, double py_1, double pz_1, double r11_1, double r12_1, double r13_1, double r21_1, double r22_1, double r23_1, double r31_1, double r32_1, double r33_1, double px_2, double py_2, double pz_2, double r11_2, double r12_2, double r13_2, double r21_2, double r22_2, double r23_2, double r31_2, double r32_2, double r33_2, double px_m, double py_m, double pz_m, double r11_m, double r12_m, double r13_m, double r21_m, double r22_m, double r23_m, double r31_m, double r32_m, double r33_m)
{delta.coeffRef(0,0) = px_1*r11_1 - px_2*r11_1 + px_m + py_1*r21_1 - py_2*r21_1 + pz_1*r31_1 - pz_2*r31_1;
delta.coeffRef(1,0) = px_1*r12_1 - px_2*r12_1 + py_1*r22_1 - py_2*r22_1 + py_m + pz_1*r32_1 - pz_2*r32_1;
delta.coeffRef(2,0) = px_1*r13_1 - px_2*r13_1 + py_1*r23_1 - py_2*r23_1 + pz_1*r33_1 - pz_2*r33_1 + pz_m;
delta.coeffRef(3,0) = -r11_1*r11_2 + r11_m - r21_1*r21_2 - r31_1*r31_2;
delta.coeffRef(4,0) = -r11_1*r12_2 + r12_m - r21_1*r22_2 - r31_1*r32_2;
delta.coeffRef(5,0) = -r11_1*r13_2 + r13_m - r21_1*r23_2 - r31_1*r33_2;
delta.coeffRef(6,0) = -r11_2*r12_1 - r21_2*r22_1 + r21_m - r31_2*r32_1;
delta.coeffRef(7,0) = -r12_1*r12_2 - r22_1*r22_2 + r22_m - r32_1*r32_2;
delta.coeffRef(8,0) = -r12_1*r13_2 - r22_1*r23_2 + r23_m - r32_1*r33_2;
delta.coeffRef(9,0) = -r11_2*r13_1 - r21_2*r23_1 - r31_2*r33_1 + r31_m;
delta.coeffRef(10,0) = -r12_2*r13_1 - r22_2*r23_1 - r32_2*r33_1 + r32_m;
delta.coeffRef(11,0) = -r13_1*r13_2 - r23_1*r23_2 - r33_1*r33_2 + r33_m;
}
inline void relative_pose_obs_eq_wc_jacobian(Eigen::Matrix<double, 12, 24, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double r11_1, double r12_1, double r13_1, double r21_1, double r22_1, double r23_1, double r31_1, double r32_1, double r33_1, double px_2, double py_2, double pz_2, double r11_2, double r12_2, double r13_2, double r21_2, double r22_2, double r23_2, double r31_2, double r32_2, double r33_2)
{j.coeffRef(0,0) = r11_1;
j.coeffRef(0,1) = r21_1;
j.coeffRef(0,2) = r31_1;
j.coeffRef(0,3) = px_1 - px_2;
j.coeffRef(0,4) = 0;
j.coeffRef(0,5) = 0;
j.coeffRef(0,6) = py_1 - py_2;
j.coeffRef(0,7) = 0;
j.coeffRef(0,8) = 0;
j.coeffRef(0,9) = pz_1 - pz_2;
j.coeffRef(0,10) = 0;
j.coeffRef(0,11) = 0;
j.coeffRef(0,12) = -r11_1;
j.coeffRef(0,13) = -r21_1;
j.coeffRef(0,14) = -r31_1;
j.coeffRef(0,15) = 0;
j.coeffRef(0,16) = 0;
j.coeffRef(0,17) = 0;
j.coeffRef(0,18) = 0;
j.coeffRef(0,19) = 0;
j.coeffRef(0,20) = 0;
j.coeffRef(0,21) = 0;
j.coeffRef(0,22) = 0;
j.coeffRef(0,23) = 0;
j.coeffRef(1,0) = r12_1;
j.coeffRef(1,1) = r22_1;
j.coeffRef(1,2) = r32_1;
j.coeffRef(1,3) = 0;
j.coeffRef(1,4) = px_1 - px_2;
j.coeffRef(1,5) = 0;
j.coeffRef(1,6) = 0;
j.coeffRef(1,7) = py_1 - py_2;
j.coeffRef(1,8) = 0;
j.coeffRef(1,9) = 0;
j.coeffRef(1,10) = pz_1 - pz_2;
j.coeffRef(1,11) = 0;
j.coeffRef(1,12) = -r12_1;
j.coeffRef(1,13) = -r22_1;
j.coeffRef(1,14) = -r32_1;
j.coeffRef(1,15) = 0;
j.coeffRef(1,16) = 0;
j.coeffRef(1,17) = 0;
j.coeffRef(1,18) = 0;
j.coeffRef(1,19) = 0;
j.coeffRef(1,20) = 0;
j.coeffRef(1,21) = 0;
j.coeffRef(1,22) = 0;
j.coeffRef(1,23) = 0;
j.coeffRef(2,0) = r13_1;
j.coeffRef(2,1) = r23_1;
j.coeffRef(2,2) = r33_1;
j.coeffRef(2,3) = 0;
j.coeffRef(2,4) = 0;
j.coeffRef(2,5) = px_1 - px_2;
j.coeffRef(2,6) = 0;
j.coeffRef(2,7) = 0;
j.coeffRef(2,8) = py_1 - py_2;
j.coeffRef(2,9) = 0;
j.coeffRef(2,10) = 0;
j.coeffRef(2,11) = pz_1 - pz_2;
j.coeffRef(2,12) = -r13_1;
j.coeffRef(2,13) = -r23_1;
j.coeffRef(2,14) = -r33_1;
j.coeffRef(2,15) = 0;
j.coeffRef(2,16) = 0;
j.coeffRef(2,17) = 0;
j.coeffRef(2,18) = 0;
j.coeffRef(2,19) = 0;
j.coeffRef(2,20) = 0;
j.coeffRef(2,21) = 0;
j.coeffRef(2,22) = 0;
j.coeffRef(2,23) = 0;
j.coeffRef(3,0) = 0;
j.coeffRef(3,1) = 0;
j.coeffRef(3,2) = 0;
j.coeffRef(3,3) = -r11_2;
j.coeffRef(3,4) = 0;
j.coeffRef(3,5) = 0;
j.coeffRef(3,6) = -r21_2;
j.coeffRef(3,7) = 0;
j.coeffRef(3,8) = 0;
j.coeffRef(3,9) = -r31_2;
j.coeffRef(3,10) = 0;
j.coeffRef(3,11) = 0;
j.coeffRef(3,12) = 0;
j.coeffRef(3,13) = 0;
j.coeffRef(3,14) = 0;
j.coeffRef(3,15) = -r11_1;
j.coeffRef(3,16) = 0;
j.coeffRef(3,17) = 0;
j.coeffRef(3,18) = -r21_1;
j.coeffRef(3,19) = 0;
j.coeffRef(3,20) = 0;
j.coeffRef(3,21) = -r31_1;
j.coeffRef(3,22) = 0;
j.coeffRef(3,23) = 0;
j.coeffRef(4,0) = 0;
j.coeffRef(4,1) = 0;
j.coeffRef(4,2) = 0;
j.coeffRef(4,3) = -r12_2;
j.coeffRef(4,4) = 0;
j.coeffRef(4,5) = 0;
j.coeffRef(4,6) = -r22_2;
j.coeffRef(4,7) = 0;
j.coeffRef(4,8) = 0;
j.coeffRef(4,9) = -r32_2;
j.coeffRef(4,10) = 0;
j.coeffRef(4,11) = 0;
j.coeffRef(4,12) = 0;
j.coeffRef(4,13) = 0;
j.coeffRef(4,14) = 0;
j.coeffRef(4,15) = 0;
j.coeffRef(4,16) = -r11_1;
j.coeffRef(4,17) = 0;
j.coeffRef(4,18) = 0;
j.coeffRef(4,19) = -r21_1;
j.coeffRef(4,20) = 0;
j.coeffRef(4,21) = 0;
j.coeffRef(4,22) = -r31_1;
j.coeffRef(4,23) = 0;
j.coeffRef(5,0) = 0;
j.coeffRef(5,1) = 0;
j.coeffRef(5,2) = 0;
j.coeffRef(5,3) = -r13_2;
j.coeffRef(5,4) = 0;
j.coeffRef(5,5) = 0;
j.coeffRef(5,6) = -r23_2;
j.coeffRef(5,7) = 0;
j.coeffRef(5,8) = 0;
j.coeffRef(5,9) = -r33_2;
j.coeffRef(5,10) = 0;
j.coeffRef(5,11) = 0;
j.coeffRef(5,12) = 0;
j.coeffRef(5,13) = 0;
j.coeffRef(5,14) = 0;
j.coeffRef(5,15) = 0;
j.coeffRef(5,16) = 0;
j.coeffRef(5,17) = -r11_1;
j.coeffRef(5,18) = 0;
j.coeffRef(5,19) = 0;
j.coeffRef(5,20) = -r21_1;
j.coeffRef(5,21) = 0;
j.coeffRef(5,22) = 0;
j.coeffRef(5,23) = -r31_1;
j.coeffRef(6,0) = 0;
j.coeffRef(6,1) = 0;
j.coeffRef(6,2) = 0;
j.coeffRef(6,3) = 0;
j.coeffRef(6,4) = -r11_2;
j.coeffRef(6,5) = 0;
j.coeffRef(6,6) = 0;
j.coeffRef(6,7) = -r21_2;
j.coeffRef(6,8) = 0;
j.coeffRef(6,9) = 0;
j.coeffRef(6,10) = -r31_2;
j.coeffRef(6,11) = 0;
j.coeffRef(6,12) = 0;
j.coeffRef(6,13) = 0;
j.coeffRef(6,14) = 0;
j.coeffRef(6,15) = -r12_1;
j.coeffRef(6,16) = 0;
j.coeffRef(6,17) = 0;
j.coeffRef(6,18) = -r22_1;
j.coeffRef(6,19) = 0;
j.coeffRef(6,20) = 0;
j.coeffRef(6,21) = -r32_1;
j.coeffRef(6,22) = 0;
j.coeffRef(6,23) = 0;
j.coeffRef(7,0) = 0;
j.coeffRef(7,1) = 0;
j.coeffRef(7,2) = 0;
j.coeffRef(7,3) = 0;
j.coeffRef(7,4) = -r12_2;
j.coeffRef(7,5) = 0;
j.coeffRef(7,6) = 0;
j.coeffRef(7,7) = -r22_2;
j.coeffRef(7,8) = 0;
j.coeffRef(7,9) = 0;
j.coeffRef(7,10) = -r32_2;
j.coeffRef(7,11) = 0;
j.coeffRef(7,12) = 0;
j.coeffRef(7,13) = 0;
j.coeffRef(7,14) = 0;
j.coeffRef(7,15) = 0;
j.coeffRef(7,16) = -r12_1;
j.coeffRef(7,17) = 0;
j.coeffRef(7,18) = 0;
j.coeffRef(7,19) = -r22_1;
j.coeffRef(7,20) = 0;
j.coeffRef(7,21) = 0;
j.coeffRef(7,22) = -r32_1;
j.coeffRef(7,23) = 0;
j.coeffRef(8,0) = 0;
j.coeffRef(8,1) = 0;
j.coeffRef(8,2) = 0;
j.coeffRef(8,3) = 0;
j.coeffRef(8,4) = -r13_2;
j.coeffRef(8,5) = 0;
j.coeffRef(8,6) = 0;
j.coeffRef(8,7) = -r23_2;
j.coeffRef(8,8) = 0;
j.coeffRef(8,9) = 0;
j.coeffRef(8,10) = -r33_2;
j.coeffRef(8,11) = 0;
j.coeffRef(8,12) = 0;
j.coeffRef(8,13) = 0;
j.coeffRef(8,14) = 0;
j.coeffRef(8,15) = 0;
j.coeffRef(8,16) = 0;
j.coeffRef(8,17) = -r12_1;
j.coeffRef(8,18) = 0;
j.coeffRef(8,19) = 0;
j.coeffRef(8,20) = -r22_1;
j.coeffRef(8,21) = 0;
j.coeffRef(8,22) = 0;
j.coeffRef(8,23) = -r32_1;
j.coeffRef(9,0) = 0;
j.coeffRef(9,1) = 0;
j.coeffRef(9,2) = 0;
j.coeffRef(9,3) = 0;
j.coeffRef(9,4) = 0;
j.coeffRef(9,5) = -r11_2;
j.coeffRef(9,6) = 0;
j.coeffRef(9,7) = 0;
j.coeffRef(9,8) = -r21_2;
j.coeffRef(9,9) = 0;
j.coeffRef(9,10) = 0;
j.coeffRef(9,11) = -r31_2;
j.coeffRef(9,12) = 0;
j.coeffRef(9,13) = 0;
j.coeffRef(9,14) = 0;
j.coeffRef(9,15) = -r13_1;
j.coeffRef(9,16) = 0;
j.coeffRef(9,17) = 0;
j.coeffRef(9,18) = -r23_1;
j.coeffRef(9,19) = 0;
j.coeffRef(9,20) = 0;
j.coeffRef(9,21) = -r33_1;
j.coeffRef(9,22) = 0;
j.coeffRef(9,23) = 0;
j.coeffRef(10,0) = 0;
j.coeffRef(10,1) = 0;
j.coeffRef(10,2) = 0;
j.coeffRef(10,3) = 0;
j.coeffRef(10,4) = 0;
j.coeffRef(10,5) = -r12_2;
j.coeffRef(10,6) = 0;
j.coeffRef(10,7) = 0;
j.coeffRef(10,8) = -r22_2;
j.coeffRef(10,9) = 0;
j.coeffRef(10,10) = 0;
j.coeffRef(10,11) = -r32_2;
j.coeffRef(10,12) = 0;
j.coeffRef(10,13) = 0;
j.coeffRef(10,14) = 0;
j.coeffRef(10,15) = 0;
j.coeffRef(10,16) = -r13_1;
j.coeffRef(10,17) = 0;
j.coeffRef(10,18) = 0;
j.coeffRef(10,19) = -r23_1;
j.coeffRef(10,20) = 0;
j.coeffRef(10,21) = 0;
j.coeffRef(10,22) = -r33_1;
j.coeffRef(10,23) = 0;
j.coeffRef(11,0) = 0;
j.coeffRef(11,1) = 0;
j.coeffRef(11,2) = 0;
j.coeffRef(11,3) = 0;
j.coeffRef(11,4) = 0;
j.coeffRef(11,5) = -r13_2;
j.coeffRef(11,6) = 0;
j.coeffRef(11,7) = 0;
j.coeffRef(11,8) = -r23_2;
j.coeffRef(11,9) = 0;
j.coeffRef(11,10) = 0;
j.coeffRef(11,11) = -r33_2;
j.coeffRef(11,12) = 0;
j.coeffRef(11,13) = 0;
j.coeffRef(11,14) = 0;
j.coeffRef(11,15) = 0;
j.coeffRef(11,16) = 0;
j.coeffRef(11,17) = -r13_1;
j.coeffRef(11,18) = 0;
j.coeffRef(11,19) = 0;
j.coeffRef(11,20) = -r23_1;
j.coeffRef(11,21) = 0;
j.coeffRef(11,22) = 0;
j.coeffRef(11,23) = -r33_1;
}inline void relative_pose_wc(Eigen::Matrix<double, 12, 1> &relative_pose, double px_1, double py_1, double pz_1, double r11_1, double r12_1, double r13_1, double r21_1, double r22_1, double r23_1, double r31_1, double r32_1, double r33_1, double px_2, double py_2, double pz_2, double r11_2, double r12_2, double r13_2, double r21_2, double r22_2, double r23_2, double r31_2, double r32_2, double r33_2)
{relative_pose.coeffRef(0,0) = -px_1*r11_1 + px_2*r11_1 - py_1*r21_1 + py_2*r21_1 - pz_1*r31_1 + pz_2*r31_1;
relative_pose.coeffRef(1,0) = -px_1*r12_1 + px_2*r12_1 - py_1*r22_1 + py_2*r22_1 - pz_1*r32_1 + pz_2*r32_1;
relative_pose.coeffRef(2,0) = -px_1*r13_1 + px_2*r13_1 - py_1*r23_1 + py_2*r23_1 - pz_1*r33_1 + pz_2*r33_1;
relative_pose.coeffRef(3,0) = r11_1*r11_2 + r21_1*r21_2 + r31_1*r31_2;
relative_pose.coeffRef(4,0) = r11_1*r12_2 + r21_1*r22_2 + r31_1*r32_2;
relative_pose.coeffRef(5,0) = r11_1*r13_2 + r21_1*r23_2 + r31_1*r33_2;
relative_pose.coeffRef(6,0) = r11_2*r12_1 + r21_2*r22_1 + r31_2*r32_1;
relative_pose.coeffRef(7,0) = r12_1*r12_2 + r22_1*r22_2 + r32_1*r32_2;
relative_pose.coeffRef(8,0) = r12_1*r13_2 + r22_1*r23_2 + r32_1*r33_2;
relative_pose.coeffRef(9,0) = r11_2*r13_1 + r21_2*r23_1 + r31_2*r33_1;
relative_pose.coeffRef(10,0) = r12_2*r13_1 + r22_2*r23_1 + r32_2*r33_1;
relative_pose.coeffRef(11,0) = r13_1*r13_2 + r23_1*r23_2 + r33_1*r33_2;
}
