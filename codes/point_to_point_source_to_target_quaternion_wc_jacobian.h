inline void point_to_point_source_to_target_quaternion_wc(double &delta_x, double &delta_y, double &delta_z, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)
{delta_x = -px - x_s*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) + x_t - y_s*(-2.0*q0*q3 + 2.0*q1*q2) - z_s*(2.0*q0*q2 + 2.0*q1*q3);
delta_y = -py - x_s*(2.0*q0*q3 + 2.0*q1*q2) - y_s*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) + y_t - z_s*(-2.0*q0*q1 + 2.0*q2*q3);
delta_z = -pz - x_s*(-2.0*q0*q2 + 2.0*q1*q3) - y_s*(2.0*q0*q1 + 2.0*q2*q3) - z_s*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1) + z_t;
}
inline void point_to_point_source_to_target_quaternion_wc_jacobian(Eigen::Matrix<double, 3, 7, Eigen::RowMajor> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_s, double y_s, double z_s)
{j.coeffRef(0,0) = -1;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = -2.0*q2*z_s + 2.0*q3*y_s;
j.coeffRef(0,4) = -2.0*q2*y_s - 2.0*q3*z_s;
j.coeffRef(0,5) = -2.0*q0*z_s - 2.0*q1*y_s + 4*q2*x_s;
j.coeffRef(0,6) = 2.0*q0*y_s - 2.0*q1*z_s + 4*q3*x_s;
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = -1;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = 2.0*q1*z_s - 2.0*q3*x_s;
j.coeffRef(1,4) = 2.0*q0*z_s + 4*q1*y_s - 2.0*q2*x_s;
j.coeffRef(1,5) = -2.0*q1*x_s - 2.0*q3*z_s;
j.coeffRef(1,6) = -2.0*q0*x_s - 2.0*q2*z_s + 4*q3*y_s;
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = -1;
j.coeffRef(2,3) = -2.0*q1*y_s + 2.0*q2*x_s;
j.coeffRef(2,4) = -2.0*q0*y_s + 4*q1*z_s - 2.0*q3*x_s;
j.coeffRef(2,5) = 2.0*q0*x_s + 4*q2*z_s - 2.0*q3*y_s;
j.coeffRef(2,6) = -2.0*q1*x_s - 2.0*q2*y_s;
}