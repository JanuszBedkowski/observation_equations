inline void point_to_projection_onto_plane_quaternion_wc(Eigen::Matrix<double, 3, 1> &delta, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double a, double b, double c)
{delta.coeffRef(0,0) = -a*(-a*x_trg_g + a*(px + x_src_l*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) + y_src_l*(-2.0*q0*q3 + 2.0*q1*q2) + z_src_l*(2.0*q0*q2 + 2.0*q1*q3)) - b*y_trg_g + b*(py + x_src_l*(2.0*q0*q3 + 2.0*q1*q2) + y_src_l*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) + z_src_l*(-2.0*q0*q1 + 2.0*q2*q3)) - c*z_trg_g + c*(pz + x_src_l*(-2.0*q0*q2 + 2.0*q1*q3) + y_src_l*(2.0*q0*q1 + 2.0*q2*q3) + z_src_l*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1)));
delta.coeffRef(1,0) = -b*(-a*x_trg_g + a*(px + x_src_l*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) + y_src_l*(-2.0*q0*q3 + 2.0*q1*q2) + z_src_l*(2.0*q0*q2 + 2.0*q1*q3)) - b*y_trg_g + b*(py + x_src_l*(2.0*q0*q3 + 2.0*q1*q2) + y_src_l*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) + z_src_l*(-2.0*q0*q1 + 2.0*q2*q3)) - c*z_trg_g + c*(pz + x_src_l*(-2.0*q0*q2 + 2.0*q1*q3) + y_src_l*(2.0*q0*q1 + 2.0*q2*q3) + z_src_l*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1)));
delta.coeffRef(2,0) = -c*(-a*x_trg_g + a*(px + x_src_l*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) + y_src_l*(-2.0*q0*q3 + 2.0*q1*q2) + z_src_l*(2.0*q0*q2 + 2.0*q1*q3)) - b*y_trg_g + b*(py + x_src_l*(2.0*q0*q3 + 2.0*q1*q2) + y_src_l*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) + z_src_l*(-2.0*q0*q1 + 2.0*q2*q3)) - c*z_trg_g + c*(pz + x_src_l*(-2.0*q0*q2 + 2.0*q1*q3) + y_src_l*(2.0*q0*q1 + 2.0*q2*q3) + z_src_l*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1)));
}
inline void point_to_projection_onto_plane_quaternion_wc_jacobian(Eigen::Matrix<double, 3, 7> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double a, double b, double c)
{j.coeffRef(0,0) = -pow(a, 2);
j.coeffRef(0,1) = -a*b;
j.coeffRef(0,2) = -a*c;
j.coeffRef(0,3) = -a*(a*(2.0*q2*z_src_l - 2.0*q3*y_src_l) + b*(-2.0*q1*z_src_l + 2.0*q3*x_src_l) + c*(2.0*q1*y_src_l - 2.0*q2*x_src_l));
j.coeffRef(0,4) = -a*(a*(2.0*q2*y_src_l + 2.0*q3*z_src_l) + b*(-2.0*q0*z_src_l - 4*q1*y_src_l + 2.0*q2*x_src_l) + c*(2.0*q0*y_src_l - 4*q1*z_src_l + 2.0*q3*x_src_l));
j.coeffRef(0,5) = -a*(a*(2.0*q0*z_src_l + 2.0*q1*y_src_l - 4*q2*x_src_l) + b*(2.0*q1*x_src_l + 2.0*q3*z_src_l) + c*(-2.0*q0*x_src_l - 4*q2*z_src_l + 2.0*q3*y_src_l));
j.coeffRef(0,6) = -a*(a*(-2.0*q0*y_src_l + 2.0*q1*z_src_l - 4*q3*x_src_l) + b*(2.0*q0*x_src_l + 2.0*q2*z_src_l - 4*q3*y_src_l) + c*(2.0*q1*x_src_l + 2.0*q2*y_src_l));
j.coeffRef(1,0) = -a*b;
j.coeffRef(1,1) = -pow(b, 2);
j.coeffRef(1,2) = -b*c;
j.coeffRef(1,3) = -b*(a*(2.0*q2*z_src_l - 2.0*q3*y_src_l) + b*(-2.0*q1*z_src_l + 2.0*q3*x_src_l) + c*(2.0*q1*y_src_l - 2.0*q2*x_src_l));
j.coeffRef(1,4) = -b*(a*(2.0*q2*y_src_l + 2.0*q3*z_src_l) + b*(-2.0*q0*z_src_l - 4*q1*y_src_l + 2.0*q2*x_src_l) + c*(2.0*q0*y_src_l - 4*q1*z_src_l + 2.0*q3*x_src_l));
j.coeffRef(1,5) = -b*(a*(2.0*q0*z_src_l + 2.0*q1*y_src_l - 4*q2*x_src_l) + b*(2.0*q1*x_src_l + 2.0*q3*z_src_l) + c*(-2.0*q0*x_src_l - 4*q2*z_src_l + 2.0*q3*y_src_l));
j.coeffRef(1,6) = -b*(a*(-2.0*q0*y_src_l + 2.0*q1*z_src_l - 4*q3*x_src_l) + b*(2.0*q0*x_src_l + 2.0*q2*z_src_l - 4*q3*y_src_l) + c*(2.0*q1*x_src_l + 2.0*q2*y_src_l));
j.coeffRef(2,0) = -a*c;
j.coeffRef(2,1) = -b*c;
j.coeffRef(2,2) = -pow(c, 2);
j.coeffRef(2,3) = -c*(a*(2.0*q2*z_src_l - 2.0*q3*y_src_l) + b*(-2.0*q1*z_src_l + 2.0*q3*x_src_l) + c*(2.0*q1*y_src_l - 2.0*q2*x_src_l));
j.coeffRef(2,4) = -c*(a*(2.0*q2*y_src_l + 2.0*q3*z_src_l) + b*(-2.0*q0*z_src_l - 4*q1*y_src_l + 2.0*q2*x_src_l) + c*(2.0*q0*y_src_l - 4*q1*z_src_l + 2.0*q3*x_src_l));
j.coeffRef(2,5) = -c*(a*(2.0*q0*z_src_l + 2.0*q1*y_src_l - 4*q2*x_src_l) + b*(2.0*q1*x_src_l + 2.0*q3*z_src_l) + c*(-2.0*q0*x_src_l - 4*q2*z_src_l + 2.0*q3*y_src_l));
j.coeffRef(2,6) = -c*(a*(-2.0*q0*y_src_l + 2.0*q1*z_src_l - 4*q3*x_src_l) + b*(2.0*q0*x_src_l + 2.0*q2*z_src_l - 4*q3*y_src_l) + c*(2.0*q1*x_src_l + 2.0*q2*y_src_l));
}