inline void point_to_point_vo_scale_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, double scale)
{delta_x = -px_1*scale + px_2*scale - x_1*cos(fi_1)*cos(ka_1) + x_2*cos(fi_2)*cos(ka_2) + y_1*sin(ka_1)*cos(fi_1) - y_2*sin(ka_2)*cos(fi_2) - z_1*sin(fi_1) + z_2*sin(fi_2);
delta_y = -py_1*scale + py_2*scale - x_1*(sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1)) + x_2*(sin(fi_2)*sin(om_2)*cos(ka_2) + sin(ka_2)*cos(om_2)) - y_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1)) + y_2*(-sin(fi_2)*sin(ka_2)*sin(om_2) + cos(ka_2)*cos(om_2)) + z_1*sin(om_1)*cos(fi_1) - z_2*sin(om_2)*cos(fi_2);
delta_z = -pz_1*scale + pz_2*scale - x_1*(-sin(fi_1)*cos(ka_1)*cos(om_1) + sin(ka_1)*sin(om_1)) + x_2*(-sin(fi_2)*cos(ka_2)*cos(om_2) + sin(ka_2)*sin(om_2)) - y_1*(sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1)) + y_2*(sin(fi_2)*sin(ka_2)*cos(om_2) + sin(om_2)*cos(ka_2)) - z_1*cos(fi_1)*cos(om_1) + z_2*cos(fi_2)*cos(om_2);
}
inline void point_to_point_vo_scale_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 1> &j, double px_1, double py_1, double pz_1, double px_2, double py_2, double pz_2)
{j.coeffRef(0,0) = -px_1 + px_2;
j.coeffRef(1,0) = -py_1 + py_2;
j.coeffRef(2,0) = -pz_1 + pz_2;
}