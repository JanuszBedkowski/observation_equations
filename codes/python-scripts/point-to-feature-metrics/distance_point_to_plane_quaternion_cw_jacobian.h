#ifndef _distance_point_to_plane_quaternion_cw_jacobian_h_
#define _distance_point_to_plane_quaternion_cw_jacobian_h_
inline void delta_distance_point_to_plane_quaternion_cw(Eigen::Matrix<double, 1, 1> &delta, double px, double py, double pz, double q0, double q1, double q2, double q3, double x, double y, double z, double a, double b, double c, double d)
{delta.coeffRef(0,0) = -a*(px*(2*pow(q2, 2) + 2*pow(q3, 2) - 1) + py*(-2.0*q0*q3 - 2.0*q1*q2) + pz*(2.0*q0*q2 - 2.0*q1*q3) + x*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) + y*(2.0*q0*q3 + 2.0*q1*q2) + z*(-2.0*q0*q2 + 2.0*q1*q3)) - b*(px*(2.0*q0*q3 - 2.0*q1*q2) + py*(2*pow(q1, 2) + 2*pow(q3, 2) - 1) + pz*(-2.0*q0*q1 - 2.0*q2*q3) + x*(-2.0*q0*q3 + 2.0*q1*q2) + y*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) + z*(2.0*q0*q1 + 2.0*q2*q3)) - c*(px*(-2.0*q0*q2 - 2.0*q1*q3) + py*(2.0*q0*q1 - 2.0*q2*q3) + pz*(2*pow(q1, 2) + 2*pow(q2, 2) - 1) + x*(2.0*q0*q2 + 2.0*q1*q3) + y*(-2.0*q0*q1 + 2.0*q2*q3) + z*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1)) - d;
}
inline void delta_distance_point_to_plane_quaternion_cw_jacobian(Eigen::Matrix<double, 1, 7> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double x, double y, double z, double a, double b, double c, double d)
{j.coeffRef(0,0) = -a*(2*pow(q2, 2) + 2*pow(q3, 2) - 1) - b*(2.0*q0*q3 - 2.0*q1*q2) - c*(-2.0*q0*q2 - 2.0*q1*q3);
j.coeffRef(0,1) = -a*(-2.0*q0*q3 - 2.0*q1*q2) - b*(2*pow(q1, 2) + 2*pow(q3, 2) - 1) - c*(2.0*q0*q1 - 2.0*q2*q3);
j.coeffRef(0,2) = -a*(2.0*q0*q2 - 2.0*q1*q3) - b*(-2.0*q0*q1 - 2.0*q2*q3) - c*(2*pow(q1, 2) + 2*pow(q2, 2) - 1);
j.coeffRef(0,3) = -a*(-2.0*py*q3 + 2.0*pz*q2 - 2.0*q2*z + 2.0*q3*y) - b*(2.0*px*q3 - 2.0*pz*q1 + 2.0*q1*z - 2.0*q3*x) - c*(-2.0*px*q2 + 2.0*py*q1 - 2.0*q1*y + 2.0*q2*x);
j.coeffRef(0,4) = -a*(-2.0*py*q2 - 2.0*pz*q3 + 2.0*q2*y + 2.0*q3*z) - b*(-2.0*px*q2 + 4*py*q1 - 2.0*pz*q0 + 2.0*q0*z - 4*q1*y + 2.0*q2*x) - c*(-2.0*px*q3 + 2.0*py*q0 + 4*pz*q1 - 2.0*q0*y - 4*q1*z + 2.0*q3*x);
j.coeffRef(0,5) = -a*(4*px*q2 - 2.0*py*q1 + 2.0*pz*q0 - 2.0*q0*z + 2.0*q1*y - 4*q2*x) - b*(-2.0*px*q1 - 2.0*pz*q3 + 2.0*q1*x + 2.0*q3*z) - c*(-2.0*px*q0 - 2.0*py*q3 + 4*pz*q2 + 2.0*q0*x - 4*q2*z + 2.0*q3*y);
j.coeffRef(0,6) = -a*(4*px*q3 - 2.0*py*q0 - 2.0*pz*q1 + 2.0*q0*y + 2.0*q1*z - 4*q3*x) - b*(2.0*px*q0 + 4*py*q3 - 2.0*pz*q2 - 2.0*q0*x + 2.0*q2*z - 4*q3*y) - c*(-2.0*px*q1 - 2.0*py*q2 + 2.0*q1*x + 2.0*q2*y);
}
#endif