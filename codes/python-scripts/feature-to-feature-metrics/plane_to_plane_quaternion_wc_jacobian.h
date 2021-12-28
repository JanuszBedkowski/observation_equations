inline void plane_to_plane_quaternion_wc(Eigen::Matrix<double, 4, 1> &delta, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double a_1, double b_1, double c_1, double d_1, double a_2, double b_2, double c_2, double d_2)
{delta.coeffRef(0,0) = -a_1*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1) + a_2*(-2*pow(q2_2, 2) - 2*pow(q3_2, 2) + 1) - b_1*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + b_2*(-2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - c_1*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + c_2*(2.0*q0_2*q2_2 + 2.0*q1_2*q3_2);
delta.coeffRef(1,0) = -a_1*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + a_2*(2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - b_1*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1) + b_2*(-2*pow(q1_2, 2) - 2*pow(q3_2, 2) + 1) - c_1*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + c_2*(-2.0*q0_2*q1_2 + 2.0*q2_2*q3_2);
delta.coeffRef(2,0) = -a_1*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + a_2*(-2.0*q0_2*q2_2 + 2.0*q1_2*q3_2) - b_1*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + b_2*(2.0*q0_2*q1_2 + 2.0*q2_2*q3_2) - c_1*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1) + c_2*(-2*pow(q1_2, 2) - 2*pow(q2_2, 2) + 1);
delta.coeffRef(3,0) = -a_1*(px_1*(2*pow(q2_1, 2) + 2*pow(q3_1, 2) - 1) + py_1*(-2.0*q0_1*q3_1 - 2.0*q1_1*q2_1) + pz_1*(2.0*q0_1*q2_1 - 2.0*q1_1*q3_1)) + a_2*(px_2*(2*pow(q2_2, 2) + 2*pow(q3_2, 2) - 1) + py_2*(-2.0*q0_2*q3_2 - 2.0*q1_2*q2_2) + pz_2*(2.0*q0_2*q2_2 - 2.0*q1_2*q3_2)) - b_1*(px_1*(2.0*q0_1*q3_1 - 2.0*q1_1*q2_1) + py_1*(2*pow(q1_1, 2) + 2*pow(q3_1, 2) - 1) + pz_1*(-2.0*q0_1*q1_1 - 2.0*q2_1*q3_1)) + b_2*(px_2*(2.0*q0_2*q3_2 - 2.0*q1_2*q2_2) + py_2*(2*pow(q1_2, 2) + 2*pow(q3_2, 2) - 1) + pz_2*(-2.0*q0_2*q1_2 - 2.0*q2_2*q3_2)) - c_1*(px_1*(-2.0*q0_1*q2_1 - 2.0*q1_1*q3_1) + py_1*(2.0*q0_1*q1_1 - 2.0*q2_1*q3_1) + pz_1*(2*pow(q1_1, 2) + 2*pow(q2_1, 2) - 1)) + c_2*(px_2*(-2.0*q0_2*q2_2 - 2.0*q1_2*q3_2) + py_2*(2.0*q0_2*q1_2 - 2.0*q2_2*q3_2) + pz_2*(2*pow(q1_2, 2) + 2*pow(q2_2, 2) - 1)) - d_1 + d_2;
}
inline void plane_to_plane_quaternion_wc_jacobian(Eigen::Matrix<double, 4, 14, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double a_1, double b_1, double c_1, double d_1, double a_2, double b_2, double c_2, double d_2)
{j.coeffRef(0,0) = 0;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = 2.0*b_1*q3_1 - 2.0*c_1*q2_1;
j.coeffRef(0,4) = -2.0*b_1*q2_1 - 2.0*c_1*q3_1;
j.coeffRef(0,5) = 4*a_1*q2_1 - 2.0*b_1*q1_1 - 2.0*c_1*q0_1;
j.coeffRef(0,6) = 4*a_1*q3_1 + 2.0*b_1*q0_1 - 2.0*c_1*q1_1;
j.coeffRef(0,7) = 0;
j.coeffRef(0,8) = 0;
j.coeffRef(0,9) = 0;
j.coeffRef(0,10) = -2.0*b_2*q3_2 + 2.0*c_2*q2_2;
j.coeffRef(0,11) = 2.0*b_2*q2_2 + 2.0*c_2*q3_2;
j.coeffRef(0,12) = -4*a_2*q2_2 + 2.0*b_2*q1_2 + 2.0*c_2*q0_2;
j.coeffRef(0,13) = -4*a_2*q3_2 - 2.0*b_2*q0_2 + 2.0*c_2*q1_2;
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = 0;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = -2.0*a_1*q3_1 + 2.0*c_1*q1_1;
j.coeffRef(1,4) = -2.0*a_1*q2_1 + 4*b_1*q1_1 + 2.0*c_1*q0_1;
j.coeffRef(1,5) = -2.0*a_1*q1_1 - 2.0*c_1*q3_1;
j.coeffRef(1,6) = -2.0*a_1*q0_1 + 4*b_1*q3_1 - 2.0*c_1*q2_1;
j.coeffRef(1,7) = 0;
j.coeffRef(1,8) = 0;
j.coeffRef(1,9) = 0;
j.coeffRef(1,10) = 2.0*a_2*q3_2 - 2.0*c_2*q1_2;
j.coeffRef(1,11) = 2.0*a_2*q2_2 - 4*b_2*q1_2 - 2.0*c_2*q0_2;
j.coeffRef(1,12) = 2.0*a_2*q1_2 + 2.0*c_2*q3_2;
j.coeffRef(1,13) = 2.0*a_2*q0_2 - 4*b_2*q3_2 + 2.0*c_2*q2_2;
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = 0;
j.coeffRef(2,3) = 2.0*a_1*q2_1 - 2.0*b_1*q1_1;
j.coeffRef(2,4) = -2.0*a_1*q3_1 - 2.0*b_1*q0_1 + 4*c_1*q1_1;
j.coeffRef(2,5) = 2.0*a_1*q0_1 - 2.0*b_1*q3_1 + 4*c_1*q2_1;
j.coeffRef(2,6) = -2.0*a_1*q1_1 - 2.0*b_1*q2_1;
j.coeffRef(2,7) = 0;
j.coeffRef(2,8) = 0;
j.coeffRef(2,9) = 0;
j.coeffRef(2,10) = -2.0*a_2*q2_2 + 2.0*b_2*q1_2;
j.coeffRef(2,11) = 2.0*a_2*q3_2 + 2.0*b_2*q0_2 - 4*c_2*q1_2;
j.coeffRef(2,12) = -2.0*a_2*q0_2 + 2.0*b_2*q3_2 - 4*c_2*q2_2;
j.coeffRef(2,13) = 2.0*a_2*q1_2 + 2.0*b_2*q2_2;
j.coeffRef(3,0) = -a_1*(2*pow(q2_1, 2) + 2*pow(q3_1, 2) - 1) - b_1*(2.0*q0_1*q3_1 - 2.0*q1_1*q2_1) - c_1*(-2.0*q0_1*q2_1 - 2.0*q1_1*q3_1);
j.coeffRef(3,1) = -a_1*(-2.0*q0_1*q3_1 - 2.0*q1_1*q2_1) - b_1*(2*pow(q1_1, 2) + 2*pow(q3_1, 2) - 1) - c_1*(2.0*q0_1*q1_1 - 2.0*q2_1*q3_1);
j.coeffRef(3,2) = -a_1*(2.0*q0_1*q2_1 - 2.0*q1_1*q3_1) - b_1*(-2.0*q0_1*q1_1 - 2.0*q2_1*q3_1) - c_1*(2*pow(q1_1, 2) + 2*pow(q2_1, 2) - 1);
j.coeffRef(3,3) = -a_1*(-2.0*py_1*q3_1 + 2.0*pz_1*q2_1) - b_1*(2.0*px_1*q3_1 - 2.0*pz_1*q1_1) - c_1*(-2.0*px_1*q2_1 + 2.0*py_1*q1_1);
j.coeffRef(3,4) = -a_1*(-2.0*py_1*q2_1 - 2.0*pz_1*q3_1) - b_1*(-2.0*px_1*q2_1 + 4*py_1*q1_1 - 2.0*pz_1*q0_1) - c_1*(-2.0*px_1*q3_1 + 2.0*py_1*q0_1 + 4*pz_1*q1_1);
j.coeffRef(3,5) = -a_1*(4*px_1*q2_1 - 2.0*py_1*q1_1 + 2.0*pz_1*q0_1) - b_1*(-2.0*px_1*q1_1 - 2.0*pz_1*q3_1) - c_1*(-2.0*px_1*q0_1 - 2.0*py_1*q3_1 + 4*pz_1*q2_1);
j.coeffRef(3,6) = -a_1*(4*px_1*q3_1 - 2.0*py_1*q0_1 - 2.0*pz_1*q1_1) - b_1*(2.0*px_1*q0_1 + 4*py_1*q3_1 - 2.0*pz_1*q2_1) - c_1*(-2.0*px_1*q1_1 - 2.0*py_1*q2_1);
j.coeffRef(3,7) = a_2*(2*pow(q2_2, 2) + 2*pow(q3_2, 2) - 1) + b_2*(2.0*q0_2*q3_2 - 2.0*q1_2*q2_2) + c_2*(-2.0*q0_2*q2_2 - 2.0*q1_2*q3_2);
j.coeffRef(3,8) = a_2*(-2.0*q0_2*q3_2 - 2.0*q1_2*q2_2) + b_2*(2*pow(q1_2, 2) + 2*pow(q3_2, 2) - 1) + c_2*(2.0*q0_2*q1_2 - 2.0*q2_2*q3_2);
j.coeffRef(3,9) = a_2*(2.0*q0_2*q2_2 - 2.0*q1_2*q3_2) + b_2*(-2.0*q0_2*q1_2 - 2.0*q2_2*q3_2) + c_2*(2*pow(q1_2, 2) + 2*pow(q2_2, 2) - 1);
j.coeffRef(3,10) = a_2*(-2.0*py_2*q3_2 + 2.0*pz_2*q2_2) + b_2*(2.0*px_2*q3_2 - 2.0*pz_2*q1_2) + c_2*(-2.0*px_2*q2_2 + 2.0*py_2*q1_2);
j.coeffRef(3,11) = a_2*(-2.0*py_2*q2_2 - 2.0*pz_2*q3_2) + b_2*(-2.0*px_2*q2_2 + 4*py_2*q1_2 - 2.0*pz_2*q0_2) + c_2*(-2.0*px_2*q3_2 + 2.0*py_2*q0_2 + 4*pz_2*q1_2);
j.coeffRef(3,12) = a_2*(4*px_2*q2_2 - 2.0*py_2*q1_2 + 2.0*pz_2*q0_2) + b_2*(-2.0*px_2*q1_2 - 2.0*pz_2*q3_2) + c_2*(-2.0*px_2*q0_2 - 2.0*py_2*q3_2 + 4*pz_2*q2_2);
j.coeffRef(3,13) = a_2*(4*px_2*q3_2 - 2.0*py_2*q0_2 - 2.0*pz_2*q1_2) + b_2*(2.0*px_2*q0_2 + 4*py_2*q3_2 - 2.0*pz_2*q2_2) + c_2*(-2.0*px_2*q1_2 - 2.0*py_2*q2_2);
}