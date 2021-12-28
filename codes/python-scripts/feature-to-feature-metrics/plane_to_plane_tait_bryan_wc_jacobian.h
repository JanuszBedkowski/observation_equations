inline void plane_to_plane_tait_bryan_wc(Eigen::Matrix<double, 4, 1> &delta, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double a_1, double b_1, double c_1, double d_1, double a_2, double b_2, double c_2, double d_2)
{delta.coeffRef(0,0) = -a_1*cos(fi_1)*cos(ka_1) + a_2*cos(fi_2)*cos(ka_2) + b_1*sin(ka_1)*cos(fi_1) - b_2*sin(ka_2)*cos(fi_2) - c_1*sin(fi_1) + c_2*sin(fi_2);
delta.coeffRef(1,0) = -a_1*(sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1)) + a_2*(sin(fi_2)*sin(om_2)*cos(ka_2) + sin(ka_2)*cos(om_2)) - b_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1)) + b_2*(-sin(fi_2)*sin(ka_2)*sin(om_2) + cos(ka_2)*cos(om_2)) + c_1*sin(om_1)*cos(fi_1) - c_2*sin(om_2)*cos(fi_2);
delta.coeffRef(2,0) = -a_1*(-sin(fi_1)*cos(ka_1)*cos(om_1) + sin(ka_1)*sin(om_1)) + a_2*(-sin(fi_2)*cos(ka_2)*cos(om_2) + sin(ka_2)*sin(om_2)) - b_1*(sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1)) + b_2*(sin(fi_2)*sin(ka_2)*cos(om_2) + sin(om_2)*cos(ka_2)) - c_1*cos(fi_1)*cos(om_1) + c_2*cos(fi_2)*cos(om_2);
delta.coeffRef(3,0) = -a_1*(-px_1*cos(fi_1)*cos(ka_1) + py_1*(-sin(fi_1)*sin(om_1)*cos(ka_1) - sin(ka_1)*cos(om_1)) + pz_1*(sin(fi_1)*cos(ka_1)*cos(om_1) - sin(ka_1)*sin(om_1))) + a_2*(-px_2*cos(fi_2)*cos(ka_2) + py_2*(-sin(fi_2)*sin(om_2)*cos(ka_2) - sin(ka_2)*cos(om_2)) + pz_2*(sin(fi_2)*cos(ka_2)*cos(om_2) - sin(ka_2)*sin(om_2))) - b_1*(px_1*sin(ka_1)*cos(fi_1) + py_1*(sin(fi_1)*sin(ka_1)*sin(om_1) - cos(ka_1)*cos(om_1)) + pz_1*(-sin(fi_1)*sin(ka_1)*cos(om_1) - sin(om_1)*cos(ka_1))) + b_2*(px_2*sin(ka_2)*cos(fi_2) + py_2*(sin(fi_2)*sin(ka_2)*sin(om_2) - cos(ka_2)*cos(om_2)) + pz_2*(-sin(fi_2)*sin(ka_2)*cos(om_2) - sin(om_2)*cos(ka_2))) - c_1*(-px_1*sin(fi_1) + py_1*sin(om_1)*cos(fi_1) - pz_1*cos(fi_1)*cos(om_1)) + c_2*(-px_2*sin(fi_2) + py_2*sin(om_2)*cos(fi_2) - pz_2*cos(fi_2)*cos(om_2)) - d_1 + d_2;
}
inline void plane_to_plane_tait_bryan_wc_jacobian(Eigen::Matrix<double, 4, 12, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double a_1, double b_1, double c_1, double d_1, double a_2, double b_2, double c_2, double d_2)
{j.coeffRef(0,0) = 0;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = 0;
j.coeffRef(0,4) = a_1*sin(fi_1)*cos(ka_1) - b_1*sin(fi_1)*sin(ka_1) - c_1*cos(fi_1);
j.coeffRef(0,5) = a_1*sin(ka_1)*cos(fi_1) + b_1*cos(fi_1)*cos(ka_1);
j.coeffRef(0,6) = 0;
j.coeffRef(0,7) = 0;
j.coeffRef(0,8) = 0;
j.coeffRef(0,9) = 0;
j.coeffRef(0,10) = -a_2*sin(fi_2)*cos(ka_2) + b_2*sin(fi_2)*sin(ka_2) + c_2*cos(fi_2);
j.coeffRef(0,11) = -a_2*sin(ka_2)*cos(fi_2) - b_2*cos(fi_2)*cos(ka_2);
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = 0;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = -a_1*(sin(fi_1)*cos(ka_1)*cos(om_1) - sin(ka_1)*sin(om_1)) - b_1*(-sin(fi_1)*sin(ka_1)*cos(om_1) - sin(om_1)*cos(ka_1)) + c_1*cos(fi_1)*cos(om_1);
j.coeffRef(1,4) = -a_1*sin(om_1)*cos(fi_1)*cos(ka_1) + b_1*sin(ka_1)*sin(om_1)*cos(fi_1) - c_1*sin(fi_1)*sin(om_1);
j.coeffRef(1,5) = -a_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1)) - b_1*(-sin(fi_1)*sin(om_1)*cos(ka_1) - sin(ka_1)*cos(om_1));
j.coeffRef(1,6) = 0;
j.coeffRef(1,7) = 0;
j.coeffRef(1,8) = 0;
j.coeffRef(1,9) = a_2*(sin(fi_2)*cos(ka_2)*cos(om_2) - sin(ka_2)*sin(om_2)) + b_2*(-sin(fi_2)*sin(ka_2)*cos(om_2) - sin(om_2)*cos(ka_2)) - c_2*cos(fi_2)*cos(om_2);
j.coeffRef(1,10) = a_2*sin(om_2)*cos(fi_2)*cos(ka_2) - b_2*sin(ka_2)*sin(om_2)*cos(fi_2) + c_2*sin(fi_2)*sin(om_2);
j.coeffRef(1,11) = a_2*(-sin(fi_2)*sin(ka_2)*sin(om_2) + cos(ka_2)*cos(om_2)) + b_2*(-sin(fi_2)*sin(om_2)*cos(ka_2) - sin(ka_2)*cos(om_2));
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = 0;
j.coeffRef(2,3) = -a_1*(sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1)) - b_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1)) + c_1*sin(om_1)*cos(fi_1);
j.coeffRef(2,4) = a_1*cos(fi_1)*cos(ka_1)*cos(om_1) - b_1*sin(ka_1)*cos(fi_1)*cos(om_1) + c_1*sin(fi_1)*cos(om_1);
j.coeffRef(2,5) = -a_1*(sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1)) - b_1*(sin(fi_1)*cos(ka_1)*cos(om_1) - sin(ka_1)*sin(om_1));
j.coeffRef(2,6) = 0;
j.coeffRef(2,7) = 0;
j.coeffRef(2,8) = 0;
j.coeffRef(2,9) = a_2*(sin(fi_2)*sin(om_2)*cos(ka_2) + sin(ka_2)*cos(om_2)) + b_2*(-sin(fi_2)*sin(ka_2)*sin(om_2) + cos(ka_2)*cos(om_2)) - c_2*sin(om_2)*cos(fi_2);
j.coeffRef(2,10) = -a_2*cos(fi_2)*cos(ka_2)*cos(om_2) + b_2*sin(ka_2)*cos(fi_2)*cos(om_2) - c_2*sin(fi_2)*cos(om_2);
j.coeffRef(2,11) = a_2*(sin(fi_2)*sin(ka_2)*cos(om_2) + sin(om_2)*cos(ka_2)) + b_2*(sin(fi_2)*cos(ka_2)*cos(om_2) - sin(ka_2)*sin(om_2));
j.coeffRef(3,0) = a_1*cos(fi_1)*cos(ka_1) - b_1*sin(ka_1)*cos(fi_1) + c_1*sin(fi_1);
j.coeffRef(3,1) = -a_1*(-sin(fi_1)*sin(om_1)*cos(ka_1) - sin(ka_1)*cos(om_1)) - b_1*(sin(fi_1)*sin(ka_1)*sin(om_1) - cos(ka_1)*cos(om_1)) - c_1*sin(om_1)*cos(fi_1);
j.coeffRef(3,2) = -a_1*(sin(fi_1)*cos(ka_1)*cos(om_1) - sin(ka_1)*sin(om_1)) - b_1*(-sin(fi_1)*sin(ka_1)*cos(om_1) - sin(om_1)*cos(ka_1)) + c_1*cos(fi_1)*cos(om_1);
j.coeffRef(3,3) = -a_1*(py_1*(-sin(fi_1)*cos(ka_1)*cos(om_1) + sin(ka_1)*sin(om_1)) + pz_1*(-sin(fi_1)*sin(om_1)*cos(ka_1) - sin(ka_1)*cos(om_1))) - b_1*(py_1*(sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1)) + pz_1*(sin(fi_1)*sin(ka_1)*sin(om_1) - cos(ka_1)*cos(om_1))) - c_1*(py_1*cos(fi_1)*cos(om_1) + pz_1*sin(om_1)*cos(fi_1));
j.coeffRef(3,4) = -a_1*(px_1*sin(fi_1)*cos(ka_1) - py_1*sin(om_1)*cos(fi_1)*cos(ka_1) + pz_1*cos(fi_1)*cos(ka_1)*cos(om_1)) - b_1*(-px_1*sin(fi_1)*sin(ka_1) + py_1*sin(ka_1)*sin(om_1)*cos(fi_1) - pz_1*sin(ka_1)*cos(fi_1)*cos(om_1)) - c_1*(-px_1*cos(fi_1) - py_1*sin(fi_1)*sin(om_1) + pz_1*sin(fi_1)*cos(om_1));
j.coeffRef(3,5) = -a_1*(px_1*sin(ka_1)*cos(fi_1) + py_1*(sin(fi_1)*sin(ka_1)*sin(om_1) - cos(ka_1)*cos(om_1)) + pz_1*(-sin(fi_1)*sin(ka_1)*cos(om_1) - sin(om_1)*cos(ka_1))) - b_1*(px_1*cos(fi_1)*cos(ka_1) + py_1*(sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1)) + pz_1*(-sin(fi_1)*cos(ka_1)*cos(om_1) + sin(ka_1)*sin(om_1)));
j.coeffRef(3,6) = -a_2*cos(fi_2)*cos(ka_2) + b_2*sin(ka_2)*cos(fi_2) - c_2*sin(fi_2);
j.coeffRef(3,7) = a_2*(-sin(fi_2)*sin(om_2)*cos(ka_2) - sin(ka_2)*cos(om_2)) + b_2*(sin(fi_2)*sin(ka_2)*sin(om_2) - cos(ka_2)*cos(om_2)) + c_2*sin(om_2)*cos(fi_2);
j.coeffRef(3,8) = a_2*(sin(fi_2)*cos(ka_2)*cos(om_2) - sin(ka_2)*sin(om_2)) + b_2*(-sin(fi_2)*sin(ka_2)*cos(om_2) - sin(om_2)*cos(ka_2)) - c_2*cos(fi_2)*cos(om_2);
j.coeffRef(3,9) = a_2*(py_2*(-sin(fi_2)*cos(ka_2)*cos(om_2) + sin(ka_2)*sin(om_2)) + pz_2*(-sin(fi_2)*sin(om_2)*cos(ka_2) - sin(ka_2)*cos(om_2))) + b_2*(py_2*(sin(fi_2)*sin(ka_2)*cos(om_2) + sin(om_2)*cos(ka_2)) + pz_2*(sin(fi_2)*sin(ka_2)*sin(om_2) - cos(ka_2)*cos(om_2))) + c_2*(py_2*cos(fi_2)*cos(om_2) + pz_2*sin(om_2)*cos(fi_2));
j.coeffRef(3,10) = a_2*(px_2*sin(fi_2)*cos(ka_2) - py_2*sin(om_2)*cos(fi_2)*cos(ka_2) + pz_2*cos(fi_2)*cos(ka_2)*cos(om_2)) + b_2*(-px_2*sin(fi_2)*sin(ka_2) + py_2*sin(ka_2)*sin(om_2)*cos(fi_2) - pz_2*sin(ka_2)*cos(fi_2)*cos(om_2)) + c_2*(-px_2*cos(fi_2) - py_2*sin(fi_2)*sin(om_2) + pz_2*sin(fi_2)*cos(om_2));
j.coeffRef(3,11) = a_2*(px_2*sin(ka_2)*cos(fi_2) + py_2*(sin(fi_2)*sin(ka_2)*sin(om_2) - cos(ka_2)*cos(om_2)) + pz_2*(-sin(fi_2)*sin(ka_2)*cos(om_2) - sin(om_2)*cos(ka_2))) + b_2*(px_2*cos(fi_2)*cos(ka_2) + py_2*(sin(fi_2)*sin(om_2)*cos(ka_2) + sin(ka_2)*cos(om_2)) + pz_2*(-sin(fi_2)*cos(ka_2)*cos(om_2) + sin(ka_2)*sin(om_2)));
}