inline void point_to_point_vo_scale_quaternion_wc(double &delta_x, double &delta_y, double &delta_z, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, double scale)
{delta_x = -px_1*scale + px_2*scale - x_1*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1) + x_2*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1) - y_1*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + y_2*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) - z_1*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + z_2*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1);
delta_y = -py_1*scale + py_2*scale - x_1*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + x_2*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) - y_1*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1) + y_2*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1) - z_1*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + z_2*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1);
delta_z = -pz_1*scale + pz_2*scale - x_1*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + x_2*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) - y_1*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + y_2*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) - z_1*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1) + z_2*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1);
}
inline void point_to_point_vo_scale_quaternion_wc_jacobian(Eigen::Matrix<double, 3, 1> &j, double px_1, double py_1, double pz_1, double px_2, double py_2, double pz_2)
{j.coeffRef(0,0) = -px_1 + px_2;
j.coeffRef(1,0) = -py_1 + py_2;
j.coeffRef(2,0) = -pz_1 + pz_2;
}