inline void delta_distance_point_to_plane_quaternion_wc(Eigen::Matrix<double, 1, 1> &delta, double px, double py, double pz, double q0, double q1, double q2, double q3, double x, double y, double z, double a, double b, double c, double d)
{delta.coeffRef(0,0) = -a*(px + x*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) + y*(-2.0*q0*q3 + 2.0*q1*q2) + z*(2.0*q0*q2 + 2.0*q1*q3)) - b*(py + x*(2.0*q0*q3 + 2.0*q1*q2) + y*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) + z*(-2.0*q0*q1 + 2.0*q2*q3)) - c*(pz + x*(-2.0*q0*q2 + 2.0*q1*q3) + y*(2.0*q0*q1 + 2.0*q2*q3) + z*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1)) - d;
}
inline void delta_distance_point_to_plane_quaternion_wc_jacobian(Eigen::Matrix<double, 1, 7> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double x, double y, double z, double a, double b, double c, double d)
{j.coeffRef(0,0) = -a;
j.coeffRef(0,1) = -b;
j.coeffRef(0,2) = -c;
j.coeffRef(0,3) = -a*(2.0*q2*z - 2.0*q3*y) - b*(-2.0*q1*z + 2.0*q3*x) - c*(2.0*q1*y - 2.0*q2*x);
j.coeffRef(0,4) = -a*(2.0*q2*y + 2.0*q3*z) - b*(-2.0*q0*z - 4*q1*y + 2.0*q2*x) - c*(2.0*q0*y - 4*q1*z + 2.0*q3*x);
j.coeffRef(0,5) = -a*(2.0*q0*z + 2.0*q1*y - 4*q2*x) - b*(2.0*q1*x + 2.0*q3*z) - c*(-2.0*q0*x - 4*q2*z + 2.0*q3*y);
j.coeffRef(0,6) = -a*(-2.0*q0*y + 2.0*q1*z - 4*q3*x) - b*(2.0*q0*x + 2.0*q2*z - 4*q3*y) - c*(2.0*q1*x + 2.0*q2*y);
}