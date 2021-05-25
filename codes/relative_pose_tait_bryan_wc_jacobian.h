inline void relative_pose_obs_eq_tait_bryan_wc_case1(Eigen::Matrix<double, 6, 1> &delta, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double px_m, double py_m, double pz_m, double om_m, double fi_m, double ka_m)
{delta.coeffRef(0,0) = pose_m - pose_to/pose_from;
