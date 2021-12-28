inline void point_to_projection_onto_line_quaternion_wc(Eigen::Matrix<double, 3, 1> &delta, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double x_trg_ln, double y_trg_ln, double z_trg_ln)
{delta.coeffRef(0,0) = -px - x_src_l*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) + x_trg_g + x_trg_ln*(x_trg_ln*(px + x_src_l*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) - x_trg_g + y_src_l*(-2.0*q0*q3 + 2.0*q1*q2) + z_src_l*(2.0*q0*q2 + 2.0*q1*q3)) + y_trg_ln*(py + x_src_l*(2.0*q0*q3 + 2.0*q1*q2) + y_src_l*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) - y_trg_g + z_src_l*(-2.0*q0*q1 + 2.0*q2*q3)) + z_trg_ln*(pz + x_src_l*(-2.0*q0*q2 + 2.0*q1*q3) + y_src_l*(2.0*q0*q1 + 2.0*q2*q3) + z_src_l*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1) - z_trg_g))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) - y_src_l*(-2.0*q0*q3 + 2.0*q1*q2) - z_src_l*(2.0*q0*q2 + 2.0*q1*q3);
delta.coeffRef(1,0) = -py - x_src_l*(2.0*q0*q3 + 2.0*q1*q2) - y_src_l*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) + y_trg_g + y_trg_ln*(x_trg_ln*(px + x_src_l*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) - x_trg_g + y_src_l*(-2.0*q0*q3 + 2.0*q1*q2) + z_src_l*(2.0*q0*q2 + 2.0*q1*q3)) + y_trg_ln*(py + x_src_l*(2.0*q0*q3 + 2.0*q1*q2) + y_src_l*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) - y_trg_g + z_src_l*(-2.0*q0*q1 + 2.0*q2*q3)) + z_trg_ln*(pz + x_src_l*(-2.0*q0*q2 + 2.0*q1*q3) + y_src_l*(2.0*q0*q1 + 2.0*q2*q3) + z_src_l*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1) - z_trg_g))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) - z_src_l*(-2.0*q0*q1 + 2.0*q2*q3);
delta.coeffRef(2,0) = -pz - x_src_l*(-2.0*q0*q2 + 2.0*q1*q3) - y_src_l*(2.0*q0*q1 + 2.0*q2*q3) - z_src_l*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1) + z_trg_g + z_trg_ln*(x_trg_ln*(px + x_src_l*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) - x_trg_g + y_src_l*(-2.0*q0*q3 + 2.0*q1*q2) + z_src_l*(2.0*q0*q2 + 2.0*q1*q3)) + y_trg_ln*(py + x_src_l*(2.0*q0*q3 + 2.0*q1*q2) + y_src_l*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) - y_trg_g + z_src_l*(-2.0*q0*q1 + 2.0*q2*q3)) + z_trg_ln*(pz + x_src_l*(-2.0*q0*q2 + 2.0*q1*q3) + y_src_l*(2.0*q0*q1 + 2.0*q2*q3) + z_src_l*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1) - z_trg_g))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
}
inline void point_to_projection_onto_line_quaternion_wc_jacobian(Eigen::Matrix<double, 3, 7> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double x_trg_ln, double y_trg_ln, double z_trg_ln)
{j.coeffRef(0,0) = pow(x_trg_ln, 2)/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) - 1;
j.coeffRef(0,1) = x_trg_ln*y_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(0,2) = x_trg_ln*z_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(0,3) = -2.0*q2*z_src_l + 2.0*q3*y_src_l + x_trg_ln*(x_trg_ln*(2.0*q2*z_src_l - 2.0*q3*y_src_l) + y_trg_ln*(-2.0*q1*z_src_l + 2.0*q3*x_src_l) + z_trg_ln*(2.0*q1*y_src_l - 2.0*q2*x_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(0,4) = -2.0*q2*y_src_l - 2.0*q3*z_src_l + x_trg_ln*(x_trg_ln*(2.0*q2*y_src_l + 2.0*q3*z_src_l) + y_trg_ln*(-2.0*q0*z_src_l - 4*q1*y_src_l + 2.0*q2*x_src_l) + z_trg_ln*(2.0*q0*y_src_l - 4*q1*z_src_l + 2.0*q3*x_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(0,5) = -2.0*q0*z_src_l - 2.0*q1*y_src_l + 4*q2*x_src_l + x_trg_ln*(x_trg_ln*(2.0*q0*z_src_l + 2.0*q1*y_src_l - 4*q2*x_src_l) + y_trg_ln*(2.0*q1*x_src_l + 2.0*q3*z_src_l) + z_trg_ln*(-2.0*q0*x_src_l - 4*q2*z_src_l + 2.0*q3*y_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(0,6) = 2.0*q0*y_src_l - 2.0*q1*z_src_l + 4*q3*x_src_l + x_trg_ln*(x_trg_ln*(-2.0*q0*y_src_l + 2.0*q1*z_src_l - 4*q3*x_src_l) + y_trg_ln*(2.0*q0*x_src_l + 2.0*q2*z_src_l - 4*q3*y_src_l) + z_trg_ln*(2.0*q1*x_src_l + 2.0*q2*y_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(1,0) = x_trg_ln*y_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(1,1) = pow(y_trg_ln, 2)/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) - 1;
j.coeffRef(1,2) = y_trg_ln*z_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(1,3) = 2.0*q1*z_src_l - 2.0*q3*x_src_l + y_trg_ln*(x_trg_ln*(2.0*q2*z_src_l - 2.0*q3*y_src_l) + y_trg_ln*(-2.0*q1*z_src_l + 2.0*q3*x_src_l) + z_trg_ln*(2.0*q1*y_src_l - 2.0*q2*x_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(1,4) = 2.0*q0*z_src_l + 4*q1*y_src_l - 2.0*q2*x_src_l + y_trg_ln*(x_trg_ln*(2.0*q2*y_src_l + 2.0*q3*z_src_l) + y_trg_ln*(-2.0*q0*z_src_l - 4*q1*y_src_l + 2.0*q2*x_src_l) + z_trg_ln*(2.0*q0*y_src_l - 4*q1*z_src_l + 2.0*q3*x_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(1,5) = -2.0*q1*x_src_l - 2.0*q3*z_src_l + y_trg_ln*(x_trg_ln*(2.0*q0*z_src_l + 2.0*q1*y_src_l - 4*q2*x_src_l) + y_trg_ln*(2.0*q1*x_src_l + 2.0*q3*z_src_l) + z_trg_ln*(-2.0*q0*x_src_l - 4*q2*z_src_l + 2.0*q3*y_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(1,6) = -2.0*q0*x_src_l - 2.0*q2*z_src_l + 4*q3*y_src_l + y_trg_ln*(x_trg_ln*(-2.0*q0*y_src_l + 2.0*q1*z_src_l - 4*q3*x_src_l) + y_trg_ln*(2.0*q0*x_src_l + 2.0*q2*z_src_l - 4*q3*y_src_l) + z_trg_ln*(2.0*q1*x_src_l + 2.0*q2*y_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,0) = x_trg_ln*z_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,1) = y_trg_ln*z_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,2) = pow(z_trg_ln, 2)/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) - 1;
j.coeffRef(2,3) = -2.0*q1*y_src_l + 2.0*q2*x_src_l + z_trg_ln*(x_trg_ln*(2.0*q2*z_src_l - 2.0*q3*y_src_l) + y_trg_ln*(-2.0*q1*z_src_l + 2.0*q3*x_src_l) + z_trg_ln*(2.0*q1*y_src_l - 2.0*q2*x_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,4) = -2.0*q0*y_src_l + 4*q1*z_src_l - 2.0*q3*x_src_l + z_trg_ln*(x_trg_ln*(2.0*q2*y_src_l + 2.0*q3*z_src_l) + y_trg_ln*(-2.0*q0*z_src_l - 4*q1*y_src_l + 2.0*q2*x_src_l) + z_trg_ln*(2.0*q0*y_src_l - 4*q1*z_src_l + 2.0*q3*x_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,5) = 2.0*q0*x_src_l + 4*q2*z_src_l - 2.0*q3*y_src_l + z_trg_ln*(x_trg_ln*(2.0*q0*z_src_l + 2.0*q1*y_src_l - 4*q2*x_src_l) + y_trg_ln*(2.0*q1*x_src_l + 2.0*q3*z_src_l) + z_trg_ln*(-2.0*q0*x_src_l - 4*q2*z_src_l + 2.0*q3*y_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,6) = -2.0*q1*x_src_l - 2.0*q2*y_src_l + z_trg_ln*(x_trg_ln*(-2.0*q0*y_src_l + 2.0*q1*z_src_l - 4*q3*x_src_l) + y_trg_ln*(2.0*q0*x_src_l + 2.0*q2*z_src_l - 4*q3*y_src_l) + z_trg_ln*(2.0*q1*x_src_l + 2.0*q2*y_src_l))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
}