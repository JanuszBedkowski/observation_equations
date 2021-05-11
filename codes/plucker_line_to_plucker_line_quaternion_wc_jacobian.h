inline void plucker_line_to_plucker_line_quaternion_wc(Eigen::Matrix<double, 6, 1> &delta, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double mx_1, double my_1, double mz_1, double lx_1, double ly_1, double lz_1, double mx_2, double my_2, double mz_2, double lx_2, double ly_2, double lz_2)
{delta.coeffRef(0,0) = -lx_1*(py_1*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) - pz_1*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1)) + lx_2*(py_2*(-2.0*q0_2*q2_2 + 2.0*q1_2*q3_2) - pz_2*(2.0*q0_2*q3_2 + 2.0*q1_2*q2_2)) - ly_1*(py_1*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) - pz_1*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1)) + ly_2*(py_2*(2.0*q0_2*q1_2 + 2.0*q2_2*q3_2) - pz_2*(-2*pow(q1_2, 2) - 2*pow(q3_2, 2) + 1)) - lz_1*(py_1*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1) - pz_1*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1)) + lz_2*(py_2*(-2*pow(q1_2, 2) - 2*pow(q2_2, 2) + 1) - pz_2*(-2.0*q0_2*q1_2 + 2.0*q2_2*q3_2)) - mx_1*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1) + mx_2*(-2*pow(q2_2, 2) - 2*pow(q3_2, 2) + 1) - my_1*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + my_2*(-2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - mz_1*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + mz_2*(2.0*q0_2*q2_2 + 2.0*q1_2*q3_2);
delta.coeffRef(1,0) = -lx_1*(-px_1*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + pz_1*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1)) + lx_2*(-px_2*(-2.0*q0_2*q2_2 + 2.0*q1_2*q3_2) + pz_2*(-2*pow(q2_2, 2) - 2*pow(q3_2, 2) + 1)) - ly_1*(-px_1*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + pz_1*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1)) + ly_2*(-px_2*(2.0*q0_2*q1_2 + 2.0*q2_2*q3_2) + pz_2*(-2.0*q0_2*q3_2 + 2.0*q1_2*q2_2)) - lz_1*(-px_1*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1) + pz_1*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1)) + lz_2*(-px_2*(-2*pow(q1_2, 2) - 2*pow(q2_2, 2) + 1) + pz_2*(2.0*q0_2*q2_2 + 2.0*q1_2*q3_2)) - mx_1*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + mx_2*(2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - my_1*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1) + my_2*(-2*pow(q1_2, 2) - 2*pow(q3_2, 2) + 1) - mz_1*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + mz_2*(-2.0*q0_2*q1_2 + 2.0*q2_2*q3_2);
delta.coeffRef(2,0) = -lx_1*(px_1*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) - py_1*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1)) + lx_2*(px_2*(2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - py_2*(-2*pow(q2_2, 2) - 2*pow(q3_2, 2) + 1)) - ly_1*(px_1*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1) - py_1*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1)) + ly_2*(px_2*(-2*pow(q1_2, 2) - 2*pow(q3_2, 2) + 1) - py_2*(-2.0*q0_2*q3_2 + 2.0*q1_2*q2_2)) - lz_1*(px_1*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) - py_1*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1)) + lz_2*(px_2*(-2.0*q0_2*q1_2 + 2.0*q2_2*q3_2) - py_2*(2.0*q0_2*q2_2 + 2.0*q1_2*q3_2)) - mx_1*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + mx_2*(-2.0*q0_2*q2_2 + 2.0*q1_2*q3_2) - my_1*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + my_2*(2.0*q0_2*q1_2 + 2.0*q2_2*q3_2) - mz_1*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1) + mz_2*(-2*pow(q1_2, 2) - 2*pow(q2_2, 2) + 1);
delta.coeffRef(3,0) = -lx_1*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1) + lx_2*(-2*pow(q2_2, 2) - 2*pow(q3_2, 2) + 1) - ly_1*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + ly_2*(-2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - lz_1*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + lz_2*(2.0*q0_2*q2_2 + 2.0*q1_2*q3_2);
delta.coeffRef(4,0) = -lx_1*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) + lx_2*(2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) - ly_1*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1) + ly_2*(-2*pow(q1_2, 2) - 2*pow(q3_2, 2) + 1) - lz_1*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + lz_2*(-2.0*q0_2*q1_2 + 2.0*q2_2*q3_2);
delta.coeffRef(5,0) = -lx_1*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) + lx_2*(-2.0*q0_2*q2_2 + 2.0*q1_2*q3_2) - ly_1*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) + ly_2*(2.0*q0_2*q1_2 + 2.0*q2_2*q3_2) - lz_1*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1) + lz_2*(-2*pow(q1_2, 2) - 2*pow(q2_2, 2) + 1);
}
inline void plucker_line_to_plucker_line_quaternion_wc_jacobian(Eigen::Matrix<double, 6, 14, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double q0_1, double q1_1, double q2_1, double q3_1, double px_2, double py_2, double pz_2, double q0_2, double q1_2, double q2_2, double q3_2, double mx_1, double my_1, double mz_1, double lx_1, double ly_1, double lz_1, double mx_2, double my_2, double mz_2, double lx_2, double ly_2, double lz_2)
{j.coeffRef(0,0) = 0;
j.coeffRef(0,1) = -lx_1*(-2.0*q0_1*q2_1 + 2.0*q1_1*q3_1) - ly_1*(2.0*q0_1*q1_1 + 2.0*q2_1*q3_1) - lz_1*(-2*pow(q1_1, 2) - 2*pow(q2_1, 2) + 1);
j.coeffRef(0,2) = -lx_1*(-2.0*q0_1*q3_1 - 2.0*q1_1*q2_1) - ly_1*(2*pow(q1_1, 2) + 2*pow(q3_1, 2) - 1) - lz_1*(2.0*q0_1*q1_1 - 2.0*q2_1*q3_1);
j.coeffRef(0,3) = -lx_1*(-2.0*py_1*q2_1 - 2.0*pz_1*q3_1) - 2.0*ly_1*py_1*q1_1 - 2.0*lz_1*pz_1*q1_1 + 2.0*my_1*q3_1 - 2.0*mz_1*q2_1;
j.coeffRef(0,4) = -lx_1*(2.0*py_1*q3_1 - 2.0*pz_1*q2_1) - ly_1*(2.0*py_1*q0_1 + 4*pz_1*q1_1) - lz_1*(-4*py_1*q1_1 + 2.0*pz_1*q0_1) - 2.0*my_1*q2_1 - 2.0*mz_1*q3_1;
j.coeffRef(0,5) = -lx_1*(-2.0*py_1*q0_1 - 2.0*pz_1*q1_1) - 2.0*ly_1*py_1*q3_1 - lz_1*(-4*py_1*q2_1 - 2.0*pz_1*q3_1) + 4*mx_1*q2_1 - 2.0*my_1*q1_1 - 2.0*mz_1*q0_1;
j.coeffRef(0,6) = -lx_1*(2.0*py_1*q1_1 - 2.0*pz_1*q0_1) - ly_1*(2.0*py_1*q2_1 + 4*pz_1*q3_1) + 2.0*lz_1*pz_1*q2_1 + 4*mx_1*q3_1 + 2.0*my_1*q0_1 - 2.0*mz_1*q1_1;
j.coeffRef(0,7) = 0;
j.coeffRef(0,8) = lx_2*(-2.0*q0_2*q2_2 + 2.0*q1_2*q3_2) + ly_2*(2.0*q0_2*q1_2 + 2.0*q2_2*q3_2) + lz_2*(-2*pow(q1_2, 2) - 2*pow(q2_2, 2) + 1);
j.coeffRef(0,9) = lx_2*(-2.0*q0_2*q3_2 - 2.0*q1_2*q2_2) + ly_2*(2*pow(q1_2, 2) + 2*pow(q3_2, 2) - 1) + lz_2*(2.0*q0_2*q1_2 - 2.0*q2_2*q3_2);
j.coeffRef(0,10) = lx_2*(-2.0*py_2*q2_2 - 2.0*pz_2*q3_2) + 2.0*ly_2*py_2*q1_2 + 2.0*lz_2*pz_2*q1_2 - 2.0*my_2*q3_2 + 2.0*mz_2*q2_2;
j.coeffRef(0,11) = lx_2*(2.0*py_2*q3_2 - 2.0*pz_2*q2_2) + ly_2*(2.0*py_2*q0_2 + 4*pz_2*q1_2) + lz_2*(-4*py_2*q1_2 + 2.0*pz_2*q0_2) + 2.0*my_2*q2_2 + 2.0*mz_2*q3_2;
j.coeffRef(0,12) = lx_2*(-2.0*py_2*q0_2 - 2.0*pz_2*q1_2) + 2.0*ly_2*py_2*q3_2 + lz_2*(-4*py_2*q2_2 - 2.0*pz_2*q3_2) - 4*mx_2*q2_2 + 2.0*my_2*q1_2 + 2.0*mz_2*q0_2;
j.coeffRef(0,13) = lx_2*(2.0*py_2*q1_2 - 2.0*pz_2*q0_2) + ly_2*(2.0*py_2*q2_2 + 4*pz_2*q3_2) - 2.0*lz_2*pz_2*q2_2 - 4*mx_2*q3_2 - 2.0*my_2*q0_2 + 2.0*mz_2*q1_2;
j.coeffRef(1,0) = -lx_1*(2.0*q0_1*q2_1 - 2.0*q1_1*q3_1) - ly_1*(-2.0*q0_1*q1_1 - 2.0*q2_1*q3_1) - lz_1*(2*pow(q1_1, 2) + 2*pow(q2_1, 2) - 1);
j.coeffRef(1,1) = 0;
j.coeffRef(1,2) = -lx_1*(-2*pow(q2_1, 2) - 2*pow(q3_1, 2) + 1) - ly_1*(-2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) - lz_1*(2.0*q0_1*q2_1 + 2.0*q1_1*q3_1);
j.coeffRef(1,3) = -2.0*lx_1*px_1*q2_1 - ly_1*(-2.0*px_1*q1_1 - 2.0*pz_1*q3_1) - 2.0*lz_1*pz_1*q2_1 - 2.0*mx_1*q3_1 + 2.0*mz_1*q1_1;
j.coeffRef(1,4) = 2.0*lx_1*px_1*q3_1 - ly_1*(-2.0*px_1*q0_1 + 2.0*pz_1*q2_1) - lz_1*(4*px_1*q1_1 + 2.0*pz_1*q3_1) - 2.0*mx_1*q2_1 + 4*my_1*q1_1 + 2.0*mz_1*q0_1;
j.coeffRef(1,5) = -lx_1*(2.0*px_1*q0_1 - 4*pz_1*q2_1) - ly_1*(-2.0*px_1*q3_1 + 2.0*pz_1*q1_1) - lz_1*(4*px_1*q2_1 + 2.0*pz_1*q0_1) - 2.0*mx_1*q1_1 - 2.0*mz_1*q3_1;
j.coeffRef(1,6) = -lx_1*(-2.0*px_1*q1_1 - 4*pz_1*q3_1) - ly_1*(-2.0*px_1*q2_1 - 2.0*pz_1*q0_1) - 2.0*lz_1*pz_1*q1_1 - 2.0*mx_1*q0_1 + 4*my_1*q3_1 - 2.0*mz_1*q2_1;
j.coeffRef(1,7) = lx_2*(2.0*q0_2*q2_2 - 2.0*q1_2*q3_2) + ly_2*(-2.0*q0_2*q1_2 - 2.0*q2_2*q3_2) + lz_2*(2*pow(q1_2, 2) + 2*pow(q2_2, 2) - 1);
j.coeffRef(1,8) = 0;
j.coeffRef(1,9) = lx_2*(-2*pow(q2_2, 2) - 2*pow(q3_2, 2) + 1) + ly_2*(-2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) + lz_2*(2.0*q0_2*q2_2 + 2.0*q1_2*q3_2);
j.coeffRef(1,10) = 2.0*lx_2*px_2*q2_2 + ly_2*(-2.0*px_2*q1_2 - 2.0*pz_2*q3_2) + 2.0*lz_2*pz_2*q2_2 + 2.0*mx_2*q3_2 - 2.0*mz_2*q1_2;
j.coeffRef(1,11) = -2.0*lx_2*px_2*q3_2 + ly_2*(-2.0*px_2*q0_2 + 2.0*pz_2*q2_2) + lz_2*(4*px_2*q1_2 + 2.0*pz_2*q3_2) + 2.0*mx_2*q2_2 - 4*my_2*q1_2 - 2.0*mz_2*q0_2;
j.coeffRef(1,12) = lx_2*(2.0*px_2*q0_2 - 4*pz_2*q2_2) + ly_2*(-2.0*px_2*q3_2 + 2.0*pz_2*q1_2) + lz_2*(4*px_2*q2_2 + 2.0*pz_2*q0_2) + 2.0*mx_2*q1_2 + 2.0*mz_2*q3_2;
j.coeffRef(1,13) = lx_2*(-2.0*px_2*q1_2 - 4*pz_2*q3_2) + ly_2*(-2.0*px_2*q2_2 - 2.0*pz_2*q0_2) + 2.0*lz_2*pz_2*q1_2 + 2.0*mx_2*q0_2 - 4*my_2*q3_2 + 2.0*mz_2*q2_2;
j.coeffRef(2,0) = -lx_1*(2.0*q0_1*q3_1 + 2.0*q1_1*q2_1) - ly_1*(-2*pow(q1_1, 2) - 2*pow(q3_1, 2) + 1) - lz_1*(-2.0*q0_1*q1_1 + 2.0*q2_1*q3_1);
j.coeffRef(2,1) = -lx_1*(2*pow(q2_1, 2) + 2*pow(q3_1, 2) - 1) - ly_1*(2.0*q0_1*q3_1 - 2.0*q1_1*q2_1) - lz_1*(-2.0*q0_1*q2_1 - 2.0*q1_1*q3_1);
j.coeffRef(2,2) = 0;
j.coeffRef(2,3) = -2.0*lx_1*px_1*q3_1 - 2.0*ly_1*py_1*q3_1 - lz_1*(-2.0*px_1*q1_1 - 2.0*py_1*q2_1) + 2.0*mx_1*q2_1 - 2.0*my_1*q1_1;
j.coeffRef(2,4) = -2.0*lx_1*px_1*q2_1 - ly_1*(-4*px_1*q1_1 - 2.0*py_1*q2_1) - lz_1*(-2.0*px_1*q0_1 - 2.0*py_1*q3_1) - 2.0*mx_1*q3_1 - 2.0*my_1*q0_1 + 4*mz_1*q1_1;
j.coeffRef(2,5) = -lx_1*(2.0*px_1*q1_1 + 4*py_1*q2_1) + 2.0*ly_1*py_1*q1_1 - lz_1*(2.0*px_1*q3_1 - 2.0*py_1*q0_1) + 2.0*mx_1*q0_1 - 2.0*my_1*q3_1 + 4*mz_1*q2_1;
j.coeffRef(2,6) = -lx_1*(2.0*px_1*q0_1 + 4*py_1*q3_1) - ly_1*(-4*px_1*q3_1 + 2.0*py_1*q0_1) - lz_1*(2.0*px_1*q2_1 - 2.0*py_1*q1_1) - 2.0*mx_1*q1_1 - 2.0*my_1*q2_1;
j.coeffRef(2,7) = lx_2*(2.0*q0_2*q3_2 + 2.0*q1_2*q2_2) + ly_2*(-2*pow(q1_2, 2) - 2*pow(q3_2, 2) + 1) + lz_2*(-2.0*q0_2*q1_2 + 2.0*q2_2*q3_2);
j.coeffRef(2,8) = lx_2*(2*pow(q2_2, 2) + 2*pow(q3_2, 2) - 1) + ly_2*(2.0*q0_2*q3_2 - 2.0*q1_2*q2_2) + lz_2*(-2.0*q0_2*q2_2 - 2.0*q1_2*q3_2);
j.coeffRef(2,9) = 0;
j.coeffRef(2,10) = 2.0*lx_2*px_2*q3_2 + 2.0*ly_2*py_2*q3_2 + lz_2*(-2.0*px_2*q1_2 - 2.0*py_2*q2_2) - 2.0*mx_2*q2_2 + 2.0*my_2*q1_2;
j.coeffRef(2,11) = 2.0*lx_2*px_2*q2_2 + ly_2*(-4*px_2*q1_2 - 2.0*py_2*q2_2) + lz_2*(-2.0*px_2*q0_2 - 2.0*py_2*q3_2) + 2.0*mx_2*q3_2 + 2.0*my_2*q0_2 - 4*mz_2*q1_2;
j.coeffRef(2,12) = lx_2*(2.0*px_2*q1_2 + 4*py_2*q2_2) - 2.0*ly_2*py_2*q1_2 + lz_2*(2.0*px_2*q3_2 - 2.0*py_2*q0_2) - 2.0*mx_2*q0_2 + 2.0*my_2*q3_2 - 4*mz_2*q2_2;
j.coeffRef(2,13) = lx_2*(2.0*px_2*q0_2 + 4*py_2*q3_2) + ly_2*(-4*px_2*q3_2 + 2.0*py_2*q0_2) + lz_2*(2.0*px_2*q2_2 - 2.0*py_2*q1_2) + 2.0*mx_2*q1_2 + 2.0*my_2*q2_2;
j.coeffRef(3,0) = 0;
j.coeffRef(3,1) = 0;
j.coeffRef(3,2) = 0;
j.coeffRef(3,3) = 2.0*ly_1*q3_1 - 2.0*lz_1*q2_1;
j.coeffRef(3,4) = -2.0*ly_1*q2_1 - 2.0*lz_1*q3_1;
j.coeffRef(3,5) = 4*lx_1*q2_1 - 2.0*ly_1*q1_1 - 2.0*lz_1*q0_1;
j.coeffRef(3,6) = 4*lx_1*q3_1 + 2.0*ly_1*q0_1 - 2.0*lz_1*q1_1;
j.coeffRef(3,7) = 0;
j.coeffRef(3,8) = 0;
j.coeffRef(3,9) = 0;
j.coeffRef(3,10) = -2.0*ly_2*q3_2 + 2.0*lz_2*q2_2;
j.coeffRef(3,11) = 2.0*ly_2*q2_2 + 2.0*lz_2*q3_2;
j.coeffRef(3,12) = -4*lx_2*q2_2 + 2.0*ly_2*q1_2 + 2.0*lz_2*q0_2;
j.coeffRef(3,13) = -4*lx_2*q3_2 - 2.0*ly_2*q0_2 + 2.0*lz_2*q1_2;
j.coeffRef(4,0) = 0;
j.coeffRef(4,1) = 0;
j.coeffRef(4,2) = 0;
j.coeffRef(4,3) = -2.0*lx_1*q3_1 + 2.0*lz_1*q1_1;
j.coeffRef(4,4) = -2.0*lx_1*q2_1 + 4*ly_1*q1_1 + 2.0*lz_1*q0_1;
j.coeffRef(4,5) = -2.0*lx_1*q1_1 - 2.0*lz_1*q3_1;
j.coeffRef(4,6) = -2.0*lx_1*q0_1 + 4*ly_1*q3_1 - 2.0*lz_1*q2_1;
j.coeffRef(4,7) = 0;
j.coeffRef(4,8) = 0;
j.coeffRef(4,9) = 0;
j.coeffRef(4,10) = 2.0*lx_2*q3_2 - 2.0*lz_2*q1_2;
j.coeffRef(4,11) = 2.0*lx_2*q2_2 - 4*ly_2*q1_2 - 2.0*lz_2*q0_2;
j.coeffRef(4,12) = 2.0*lx_2*q1_2 + 2.0*lz_2*q3_2;
j.coeffRef(4,13) = 2.0*lx_2*q0_2 - 4*ly_2*q3_2 + 2.0*lz_2*q2_2;
j.coeffRef(5,0) = 0;
j.coeffRef(5,1) = 0;
j.coeffRef(5,2) = 0;
j.coeffRef(5,3) = 2.0*lx_1*q2_1 - 2.0*ly_1*q1_1;
j.coeffRef(5,4) = -2.0*lx_1*q3_1 - 2.0*ly_1*q0_1 + 4*lz_1*q1_1;
j.coeffRef(5,5) = 2.0*lx_1*q0_1 - 2.0*ly_1*q3_1 + 4*lz_1*q2_1;
j.coeffRef(5,6) = -2.0*lx_1*q1_1 - 2.0*ly_1*q2_1;
j.coeffRef(5,7) = 0;
j.coeffRef(5,8) = 0;
j.coeffRef(5,9) = 0;
j.coeffRef(5,10) = -2.0*lx_2*q2_2 + 2.0*ly_2*q1_2;
j.coeffRef(5,11) = 2.0*lx_2*q3_2 + 2.0*ly_2*q0_2 - 4*lz_2*q1_2;
j.coeffRef(5,12) = -2.0*lx_2*q0_2 + 2.0*ly_2*q3_2 - 4*lz_2*q2_2;
j.coeffRef(5,13) = 2.0*lx_2*q1_2 + 2.0*ly_2*q2_2;
}