inline void metric_camera_colinearity_euler_ka_fi_om_wc(double &ksi, double &eta, double ksi_0, double eta_0, double c, double px, double py, double pz, double om, double fi, double ka, double tie_px, double tie_py, double tie_pz)
{ksi = -c*((-px + tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + (-pz + tie_pz)*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + ksi_0;
eta = -c*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + eta_0;
}
inline void observation_equation_metric_camera_colinearity_euler_ka_fi_om_wc(Eigen::Matrix<double, 2, 1> &delta, double ksi_0, double eta_0, double c, double px, double py, double pz, double om, double fi, double ka, double tie_px, double tie_py, double tie_pz, double ksi_kp, double eta_kp)
{delta.coeffRef(0,0) = c*((-px + tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + (-pz + tie_pz)*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) - ksi_0 + ksi_kp;
delta.coeffRef(1,0) = c*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) - eta_0 + eta_kp;
}
inline void observation_equation_metric_camera_colinearity_euler_ka_fi_om_wc_jacobian(Eigen::Matrix<double, 2, 9, Eigen::RowMajor> &j, double ksi_0, double eta_0, double c, double px, double py, double pz, double om, double fi, double ka, double tie_px, double tie_py, double tie_pz, double ksi_kp, double eta_kp)
{j.coeffRef(0,0) = -c*cos(fi)*cos(om)/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + c*((-px + tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + (-pz + tie_pz)*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)))*sin(fi)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(0,1) = c*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) - c*((-px + tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + (-pz + tie_pz)*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)))*sin(ka)*cos(fi)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(0,2) = c*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + c*((-px + tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + (-pz + tie_pz)*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)))*cos(fi)*cos(ka)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(0,3) = c*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka));
j.coeffRef(0,4) = c*(-(-px + tie_px)*sin(fi)*cos(om) + (-py + tie_py)*sin(ka)*cos(fi)*cos(om) - (-pz + tie_pz)*cos(fi)*cos(ka)*cos(om))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + c*(-(-px + tie_px)*cos(fi) + (py - tie_py)*sin(fi)*sin(ka) + (-pz + tie_pz)*sin(fi)*cos(ka))*((-px + tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + (-pz + tie_pz)*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)))/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(0,5) = c*((-py + tie_py)*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + (-pz + tie_pz)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + c*(-(py - tie_py)*cos(fi)*cos(ka) + (-pz + tie_pz)*sin(ka)*cos(fi))*((-px + tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + (-pz + tie_pz)*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)))/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(0,6) = c*cos(fi)*cos(om)/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) - c*((-px + tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + (-pz + tie_pz)*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)))*sin(fi)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(0,7) = c*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + c*((-px + tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + (-pz + tie_pz)*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)))*sin(ka)*cos(fi)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(0,8) = c*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) - c*((-px + tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + (-pz + tie_pz)*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)))*cos(fi)*cos(ka)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(1,0) = c*sin(om)*cos(fi)/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + c*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))*sin(fi)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(1,1) = c*(sin(fi)*sin(ka)*sin(om) - cos(ka)*cos(om))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) - c*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))*sin(ka)*cos(fi)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(1,2) = c*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + c*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))*cos(fi)*cos(ka)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(1,3) = c*((px - tie_px)*cos(fi)*cos(om) + (-py + tie_py)*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) + (-pz + tie_pz)*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka));
j.coeffRef(1,4) = c*(-(px - tie_px)*sin(fi)*sin(om) - (-py + tie_py)*sin(ka)*sin(om)*cos(fi) + (-pz + tie_pz)*sin(om)*cos(fi)*cos(ka))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + c*(-(-px + tie_px)*cos(fi) + (py - tie_py)*sin(fi)*sin(ka) + (-pz + tie_pz)*sin(fi)*cos(ka))*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(1,5) = c*((-py + tie_py)*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om)) + (-pz + tie_pz)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + c*(-(py - tie_py)*cos(fi)*cos(ka) + (-pz + tie_pz)*sin(ka)*cos(fi))*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(1,6) = -c*sin(om)*cos(fi)/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) - c*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))*sin(fi)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(1,7) = c*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) + c*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))*sin(ka)*cos(fi)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
j.coeffRef(1,8) = c*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om))/((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka)) - c*(-(-px + tie_px)*sin(om)*cos(fi) + (-py + tie_py)*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + (-pz + tie_pz)*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)))*cos(fi)*cos(ka)/pow((-px + tie_px)*sin(fi) - (-py + tie_py)*sin(ka)*cos(fi) + (-pz + tie_pz)*cos(fi)*cos(ka), 2);
}