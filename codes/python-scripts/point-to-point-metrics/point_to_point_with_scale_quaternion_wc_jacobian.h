inline void point_to_point_with_scale_quaternion_wc(double &delta_x, double &delta_y, double &delta_z, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, double s_1, double s_2)
{delta_x = px_1 - px_2 + s_1*x_1*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1) + s_1*y_1*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + s_1*z_1*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) - s_2*x_2*(-2*pow(q2_2, 2) - 2*pow(q3_2, 2) + 1) - s_2*y_2*(-2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - s_2*z_2*(2.0*q0_2*q2_2 + 2.0*q1_2*q3_2);
delta_y = py_1 - py_2 + s_1*x_1*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + s_1*y_1*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1) + s_1*z_1*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) - s_2*x_2*(2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - s_2*y_2*(-2*pow(q1_2, 2) - 2*pow(q3_2, 2) + 1) - s_2*z_2*(-2.0*q0_2*q1_2 + 2.0*q2_2*q3_2);
delta_z = pz_1 - pz_2 + s_1*x_1*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + s_1*y_1*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + s_1*z_1*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1) - s_2*x_2*(-2.0*q0_2*q2_2 + 2.0*q1_2*q3_2) - s_2*y_2*(2.0*q0_2*q1_2 + 2.0*q2_2*q3_2) - s_2*z_2*(-2*pow(q1_2, 2) - 2*pow(q2_2, 2) + 1);
}
inline void point_to_point_with_scale_quaternion_wc_jacobian(Eigen::Matrix<double, 3, 16, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, double s_1, double s_2)
{j.coeffRef(0,0) = 1;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = 2.0*q2_1*s_1*z_1 - 2.0*q3_1*s_1*y_1;
j.coeffRef(0,4) = 2.0*q2_1*s_1*y_1 + 2.0*q3_1*s_1*z_1;
j.coeffRef(0,5) = 2.0*q0_1*s_1*z_1 + 2.0*q1_1*s_1*y_1 - 4*q2_1*s_1*x_1;
j.coeffRef(0,6) = -2.0*q0_1*s_1*y_1 + 2.0*q1_1*s_1*z_1 - 4*q3_1*s_1*x_1;
j.coeffRef(0,7) = x_1*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1) + y_1*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + z_1*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1);
j.coeffRef(0,8) = -1;
j.coeffRef(0,9) = 0;
j.coeffRef(0,10) = 0;
j.coeffRef(0,11) = -2.0*q2_2*s_2*z_2 + 2.0*q3_2*s_2*y_2;
j.coeffRef(0,12) = -2.0*q2_2*s_2*y_2 - 2.0*q3_2*s_2*z_2;
j.coeffRef(0,13) = -2.0*q0_2*s_2*z_2 - 2.0*q1_2*s_2*y_2 + 4*q2_2*s_2*x_2;
j.coeffRef(0,14) = 2.0*q0_2*s_2*y_2 - 2.0*q1_2*s_2*z_2 + 4*q3_2*s_2*x_2;
j.coeffRef(0,15) = -x_2*(-2*pow(q2_2, 2) - 2*pow(q3_2, 2) + 1) - y_2*(-2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - z_2*(2.0*q0_2*q2_2 + 2.0*q1_2*q3_2);
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = 1;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = -2.0*q1_1*s_1*z_1 + 2.0*q3_1*s_1*x_1;
j.coeffRef(1,4) = -2.0*q0_1*s_1*z_1 - 4*q1_1*s_1*y_1 + 2.0*q2_1*s_1*x_1;
j.coeffRef(1,5) = 2.0*q1_1*s_1*x_1 + 2.0*q3_1*s_1*z_1;
j.coeffRef(1,6) = 2.0*q0_1*s_1*x_1 + 2.0*q2_1*s_1*z_1 - 4*q3_1*s_1*y_1;
j.coeffRef(1,7) = x_1*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + y_1*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1) + z_1*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1);
j.coeffRef(1,8) = 0;
j.coeffRef(1,9) = -1;
j.coeffRef(1,10) = 0;
j.coeffRef(1,11) = 2.0*q1_2*s_2*z_2 - 2.0*q3_2*s_2*x_2;
j.coeffRef(1,12) = 2.0*q0_2*s_2*z_2 + 4*q1_2*s_2*y_2 - 2.0*q2_2*s_2*x_2;
j.coeffRef(1,13) = -2.0*q1_2*s_2*x_2 - 2.0*q3_2*s_2*z_2;
j.coeffRef(1,14) = -2.0*q0_2*s_2*x_2 - 2.0*q2_2*s_2*z_2 + 4*q3_2*s_2*y_2;
j.coeffRef(1,15) = -x_2*(2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - y_2*(-2*pow(q1_2, 2) - 2*pow(q3_2, 2) + 1) - z_2*(-2.0*q0_2*q1_2 + 2.0*q2_2*q3_2);
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = 1;
j.coeffRef(2,3) = 2.0*q1_1*s_1*y_1 - 2.0*q2_1*s_1*x_1;
j.coeffRef(2,4) = 2.0*q0_1*s_1*y_1 - 4*q1_1*s_1*z_1 + 2.0*q3_1*s_1*x_1;
j.coeffRef(2,5) = -2.0*q0_1*s_1*x_1 - 4*q2_1*s_1*z_1 + 2.0*q3_1*s_1*y_1;
j.coeffRef(2,6) = 2.0*q1_1*s_1*x_1 + 2.0*q2_1*s_1*y_1;
j.coeffRef(2,7) = x_1*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + y_1*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + z_1*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1);
j.coeffRef(2,8) = 0;
j.coeffRef(2,9) = 0;
j.coeffRef(2,10) = -1;
j.coeffRef(2,11) = -2.0*q1_2*s_2*y_2 + 2.0*q2_2*s_2*x_2;
j.coeffRef(2,12) = -2.0*q0_2*s_2*y_2 + 4*q1_2*s_2*z_2 - 2.0*q3_2*s_2*x_2;
j.coeffRef(2,13) = 2.0*q0_2*s_2*x_2 + 4*q2_2*s_2*z_2 - 2.0*q3_2*s_2*y_2;
j.coeffRef(2,14) = -2.0*q1_2*s_2*x_2 - 2.0*q2_2*s_2*y_2;
j.coeffRef(2,15) = -x_2*(-2.0*q0_2*q2_2 + 2.0*q1_2*q3_2) - y_2*(2.0*q0_2*q1_2 + 2.0*q2_2*q3_2) - z_2*(-2*pow(q1_2, 2) - 2*pow(q2_2, 2) + 1);
}