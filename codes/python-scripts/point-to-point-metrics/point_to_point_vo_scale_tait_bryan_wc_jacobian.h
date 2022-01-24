inline void point_to_point_vo_scale_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double tx_2, double ty_2, double tz_2, double om_2, double fi_2, double ka_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, double scale)
{delta_x = -scale*tx_1 + scale*tx_2 - x_1*cos(fi_1)*cos(ka_1) + x_2*cos(fi_2)*cos(ka_2) + y_1*sin(ka_1)*cos(fi_1) - y_2*sin(ka_2)*cos(fi_2) - z_1*sin(fi_1) + z_2*sin(fi_2);
delta_y = -scale*ty_1 + scale*ty_2 - x_1*(sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1)) + x_2*(sin(fi_2)*sin(om_2)*cos(ka_2) + sin(ka_2)*cos(om_2)) - y_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1)) + y_2*(-sin(fi_2)*sin(ka_2)*sin(om_2) + cos(ka_2)*cos(om_2)) + z_1*sin(om_1)*cos(fi_1) - z_2*sin(om_2)*cos(fi_2);
delta_z = -scale*tz_1 + scale*tz_2 - x_1*(-sin(fi_1)*cos(ka_1)*cos(om_1) + sin(ka_1)*sin(om_1)) + x_2*(-sin(fi_2)*cos(ka_2)*cos(om_2) + sin(ka_2)*sin(om_2)) - y_1*(sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1)) + y_2*(sin(fi_2)*sin(ka_2)*cos(om_2) + sin(om_2)*cos(ka_2)) - z_1*cos(fi_1)*cos(om_1) + z_2*cos(fi_2)*cos(om_2);
}
inline void point_to_point_vo_scale_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 1> &j, double tx_1, double ty_1, double tz_1, double tx_2, double ty_2, double tz_2)
{j.coeffRef(0,0) = -tx_1 + tx_2;
j.coeffRef(1,0) = -ty_1 + ty_2;
j.coeffRef(2,0) = -tz_1 + tz_2;
}