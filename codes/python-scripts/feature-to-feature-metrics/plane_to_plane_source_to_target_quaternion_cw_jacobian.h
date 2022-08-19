#ifndef _plane_to_plane_source_to_target_quaternion_cw_jacobian_h_
#define _plane_to_plane_source_to_target_quaternion_cw_jacobian_h_
inline void plane_to_plane_source_to_target_quaternion_cw(Eigen::Matrix<double, 4, 1> &delta, double tx_1, double ty_1, double tz_1, double q0_1, double q1_1, double q2_1, double q3_1, double a_1, double b_1, double c_1, double d_1, double a_2, double b_2, double c_2, double d_2)
{delta.coeffRef(0,0) = -a_1*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1) + a_2 - b_1*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) - c_1*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1);
delta.coeffRef(1,0) = -a_1*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) - b_1*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1) + b_2 - c_1*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1);
delta.coeffRef(2,0) = -a_1*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) - b_1*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) - c_1*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1) + c_2;
delta.coeffRef(3,0) = -a_1*tx_1 - b_1*ty_1 - c_1*tz_1 - d_1 + d_2;
}
inline void plane_to_plane_source_to_target_quaternion_cw_jacobian(Eigen::Matrix<double, 4, 7, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double q0_1, double q1_1, double q2_1, double q3_1, double a_1, double b_1, double c_1, double d_1, double a_2, double b_2, double c_2, double d_2)
{j.coeffRef(0,0) = 0;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = -2.0*b_1*q3_1 + 2.0*c_1*q2_1;
j.coeffRef(0,4) = -2.0*b_1*q2_1 - 2.0*c_1*q3_1;
j.coeffRef(0,5) = 4*a_1*q2_1 - 2.0*b_1*q1_1 + 2.0*c_1*q0_1;
j.coeffRef(0,6) = 4*a_1*q3_1 - 2.0*b_1*q0_1 - 2.0*c_1*q1_1;
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = 0;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = 2.0*a_1*q3_1 - 2.0*c_1*q1_1;
j.coeffRef(1,4) = -2.0*a_1*q2_1 + 4*b_1*q1_1 - 2.0*c_1*q0_1;
j.coeffRef(1,5) = -2.0*a_1*q1_1 - 2.0*c_1*q3_1;
j.coeffRef(1,6) = 2.0*a_1*q0_1 + 4*b_1*q3_1 - 2.0*c_1*q2_1;
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = 0;
j.coeffRef(2,3) = -2.0*a_1*q2_1 + 2.0*b_1*q1_1;
j.coeffRef(2,4) = -2.0*a_1*q3_1 + 2.0*b_1*q0_1 + 4*c_1*q1_1;
j.coeffRef(2,5) = -2.0*a_1*q0_1 - 2.0*b_1*q3_1 + 4*c_1*q2_1;
j.coeffRef(2,6) = -2.0*a_1*q1_1 - 2.0*b_1*q2_1;
j.coeffRef(3,0) = -a_1;
j.coeffRef(3,1) = -b_1;
j.coeffRef(3,2) = -c_1;
j.coeffRef(3,3) = 0;
j.coeffRef(3,4) = 0;
j.coeffRef(3,5) = 0;
j.coeffRef(3,6) = 0;
}
#endif