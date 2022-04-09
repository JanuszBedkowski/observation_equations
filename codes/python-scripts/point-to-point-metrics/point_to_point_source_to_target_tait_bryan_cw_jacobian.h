inline void point_to_point_source_to_target_tait_bryan_cw(double &delta_x, double &delta_y, double &delta_z, double tx_cw, double ty_cw, double tz_cw, double om_cw, double fi_cw, double ka_cw, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)
{delta_x = tx_cw*cos(fi_cw)*cos(ka_cw) - ty_cw*(-sin(fi_cw)*sin(om_cw)*cos(ka_cw) - sin(ka_cw)*cos(om_cw)) - tz_cw*(sin(fi_cw)*cos(ka_cw)*cos(om_cw) - sin(ka_cw)*sin(om_cw)) - x_s*cos(fi_cw)*cos(ka_cw) + x_t - y_s*(sin(fi_cw)*sin(om_cw)*cos(ka_cw) + sin(ka_cw)*cos(om_cw)) - z_s*(-sin(fi_cw)*cos(ka_cw)*cos(om_cw) + sin(ka_cw)*sin(om_cw));
delta_y = -tx_cw*sin(ka_cw)*cos(fi_cw) - ty_cw*(sin(fi_cw)*sin(ka_cw)*sin(om_cw) - cos(ka_cw)*cos(om_cw)) - tz_cw*(-sin(fi_cw)*sin(ka_cw)*cos(om_cw) - sin(om_cw)*cos(ka_cw)) + x_s*sin(ka_cw)*cos(fi_cw) - y_s*(-sin(fi_cw)*sin(ka_cw)*sin(om_cw) + cos(ka_cw)*cos(om_cw)) + y_t - z_s*(sin(fi_cw)*sin(ka_cw)*cos(om_cw) + sin(om_cw)*cos(ka_cw));
delta_z = tx_cw*sin(fi_cw) - ty_cw*sin(om_cw)*cos(fi_cw) + tz_cw*cos(fi_cw)*cos(om_cw) - x_s*sin(fi_cw) + y_s*sin(om_cw)*cos(fi_cw) - z_s*cos(fi_cw)*cos(om_cw) + z_t;
}
inline void point_to_point_source_to_target_tait_bryan_cw_jacobian(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double tx_cw, double ty_cw, double tz_cw, double om_cw, double fi_cw, double ka_cw, double x_s, double y_s, double z_s)
{j.coeffRef(0,0) = cos(fi_cw)*cos(ka_cw);
j.coeffRef(0,1) = sin(fi_cw)*sin(om_cw)*cos(ka_cw) + sin(ka_cw)*cos(om_cw);
j.coeffRef(0,2) = -sin(fi_cw)*cos(ka_cw)*cos(om_cw) + sin(ka_cw)*sin(om_cw);
j.coeffRef(0,3) = -ty_cw*(-sin(fi_cw)*cos(ka_cw)*cos(om_cw) + sin(ka_cw)*sin(om_cw)) - tz_cw*(-sin(fi_cw)*sin(om_cw)*cos(ka_cw) - sin(ka_cw)*cos(om_cw)) - y_s*(sin(fi_cw)*cos(ka_cw)*cos(om_cw) - sin(ka_cw)*sin(om_cw)) - z_s*(sin(fi_cw)*sin(om_cw)*cos(ka_cw) + sin(ka_cw)*cos(om_cw));
j.coeffRef(0,4) = -tx_cw*sin(fi_cw)*cos(ka_cw) + ty_cw*sin(om_cw)*cos(fi_cw)*cos(ka_cw) - tz_cw*cos(fi_cw)*cos(ka_cw)*cos(om_cw) + x_s*sin(fi_cw)*cos(ka_cw) - y_s*sin(om_cw)*cos(fi_cw)*cos(ka_cw) + z_s*cos(fi_cw)*cos(ka_cw)*cos(om_cw);
j.coeffRef(0,5) = -tx_cw*sin(ka_cw)*cos(fi_cw) - ty_cw*(sin(fi_cw)*sin(ka_cw)*sin(om_cw) - cos(ka_cw)*cos(om_cw)) - tz_cw*(-sin(fi_cw)*sin(ka_cw)*cos(om_cw) - sin(om_cw)*cos(ka_cw)) + x_s*sin(ka_cw)*cos(fi_cw) - y_s*(-sin(fi_cw)*sin(ka_cw)*sin(om_cw) + cos(ka_cw)*cos(om_cw)) - z_s*(sin(fi_cw)*sin(ka_cw)*cos(om_cw) + sin(om_cw)*cos(ka_cw));
j.coeffRef(1,0) = -sin(ka_cw)*cos(fi_cw);
j.coeffRef(1,1) = -sin(fi_cw)*sin(ka_cw)*sin(om_cw) + cos(ka_cw)*cos(om_cw);
j.coeffRef(1,2) = sin(fi_cw)*sin(ka_cw)*cos(om_cw) + sin(om_cw)*cos(ka_cw);
j.coeffRef(1,3) = -ty_cw*(sin(fi_cw)*sin(ka_cw)*cos(om_cw) + sin(om_cw)*cos(ka_cw)) - tz_cw*(sin(fi_cw)*sin(ka_cw)*sin(om_cw) - cos(ka_cw)*cos(om_cw)) - y_s*(-sin(fi_cw)*sin(ka_cw)*cos(om_cw) - sin(om_cw)*cos(ka_cw)) - z_s*(-sin(fi_cw)*sin(ka_cw)*sin(om_cw) + cos(ka_cw)*cos(om_cw));
j.coeffRef(1,4) = tx_cw*sin(fi_cw)*sin(ka_cw) - ty_cw*sin(ka_cw)*sin(om_cw)*cos(fi_cw) + tz_cw*sin(ka_cw)*cos(fi_cw)*cos(om_cw) - x_s*sin(fi_cw)*sin(ka_cw) + y_s*sin(ka_cw)*sin(om_cw)*cos(fi_cw) - z_s*sin(ka_cw)*cos(fi_cw)*cos(om_cw);
j.coeffRef(1,5) = -tx_cw*cos(fi_cw)*cos(ka_cw) - ty_cw*(sin(fi_cw)*sin(om_cw)*cos(ka_cw) + sin(ka_cw)*cos(om_cw)) - tz_cw*(-sin(fi_cw)*cos(ka_cw)*cos(om_cw) + sin(ka_cw)*sin(om_cw)) + x_s*cos(fi_cw)*cos(ka_cw) - y_s*(-sin(fi_cw)*sin(om_cw)*cos(ka_cw) - sin(ka_cw)*cos(om_cw)) - z_s*(sin(fi_cw)*cos(ka_cw)*cos(om_cw) - sin(ka_cw)*sin(om_cw));
j.coeffRef(2,0) = sin(fi_cw);
j.coeffRef(2,1) = -sin(om_cw)*cos(fi_cw);
j.coeffRef(2,2) = cos(fi_cw)*cos(om_cw);
j.coeffRef(2,3) = -ty_cw*cos(fi_cw)*cos(om_cw) - tz_cw*sin(om_cw)*cos(fi_cw) + y_s*cos(fi_cw)*cos(om_cw) + z_s*sin(om_cw)*cos(fi_cw);
j.coeffRef(2,4) = tx_cw*cos(fi_cw) + ty_cw*sin(fi_cw)*sin(om_cw) - tz_cw*sin(fi_cw)*cos(om_cw) - x_s*cos(fi_cw) - y_s*sin(fi_cw)*sin(om_cw) + z_s*sin(fi_cw)*cos(om_cw);
j.coeffRef(2,5) = 0;
}inline void point_to_point_source_to_target_tait_bryan_cw_d2sum_dbeta2(Eigen::Matrix<double, 3, 3, Eigen::RowMajor> &j, double tx_cw, double ty_cw, double tz_cw, double om_cw, double fi_cw, double ka_cw, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)
{j.coeffRef(0,0) = 2*pow(sin(fi_cw), 2) + 2*pow(sin(ka_cw), 2)*pow(cos(fi_cw), 2) + 2*pow(cos(fi_cw), 2)*pow(cos(ka_cw), 2);
j.coeffRef(0,1) = (2*sin(fi_cw)*sin(ka_cw)*sin(om_cw) - 2*cos(ka_cw)*cos(om_cw))*sin(ka_cw)*cos(fi_cw) + (2*sin(fi_cw)*sin(om_cw)*cos(ka_cw) + 2*sin(ka_cw)*cos(om_cw))*cos(fi_cw)*cos(ka_cw) - 2*sin(fi_cw)*sin(om_cw)*cos(fi_cw);
j.coeffRef(0,2) = (-2*sin(fi_cw)*sin(ka_cw)*cos(om_cw) - 2*sin(om_cw)*cos(ka_cw))*sin(ka_cw)*cos(fi_cw) + (-2*sin(fi_cw)*cos(ka_cw)*cos(om_cw) + 2*sin(ka_cw)*sin(om_cw))*cos(fi_cw)*cos(ka_cw) + 2*sin(fi_cw)*cos(fi_cw)*cos(om_cw);
j.coeffRef(1,0) = -(-2*sin(fi_cw)*sin(ka_cw)*sin(om_cw) + 2*cos(ka_cw)*cos(om_cw))*sin(ka_cw)*cos(fi_cw) + (2*sin(fi_cw)*sin(om_cw)*cos(ka_cw) + 2*sin(ka_cw)*cos(om_cw))*cos(fi_cw)*cos(ka_cw) - 2*sin(fi_cw)*sin(om_cw)*cos(fi_cw);
j.coeffRef(1,1) = (-2*sin(fi_cw)*sin(ka_cw)*sin(om_cw) + 2*cos(ka_cw)*cos(om_cw))*(-sin(fi_cw)*sin(ka_cw)*sin(om_cw) + cos(ka_cw)*cos(om_cw)) + (sin(fi_cw)*sin(om_cw)*cos(ka_cw) + sin(ka_cw)*cos(om_cw))*(2*sin(fi_cw)*sin(om_cw)*cos(ka_cw) + 2*sin(ka_cw)*cos(om_cw)) + 2*pow(sin(om_cw), 2)*pow(cos(fi_cw), 2);
j.coeffRef(1,2) = (-2*sin(fi_cw)*sin(ka_cw)*sin(om_cw) + 2*cos(ka_cw)*cos(om_cw))*(sin(fi_cw)*sin(ka_cw)*cos(om_cw) + sin(om_cw)*cos(ka_cw)) + (2*sin(fi_cw)*sin(om_cw)*cos(ka_cw) + 2*sin(ka_cw)*cos(om_cw))*(-sin(fi_cw)*cos(ka_cw)*cos(om_cw) + sin(ka_cw)*sin(om_cw)) - 2*sin(om_cw)*pow(cos(fi_cw), 2)*cos(om_cw);
j.coeffRef(2,0) = -(2*sin(fi_cw)*sin(ka_cw)*cos(om_cw) + 2*sin(om_cw)*cos(ka_cw))*sin(ka_cw)*cos(fi_cw) + (-2*sin(fi_cw)*cos(ka_cw)*cos(om_cw) + 2*sin(ka_cw)*sin(om_cw))*cos(fi_cw)*cos(ka_cw) + 2*sin(fi_cw)*cos(fi_cw)*cos(om_cw);
j.coeffRef(2,1) = (-sin(fi_cw)*sin(ka_cw)*sin(om_cw) + cos(ka_cw)*cos(om_cw))*(2*sin(fi_cw)*sin(ka_cw)*cos(om_cw) + 2*sin(om_cw)*cos(ka_cw)) + (sin(fi_cw)*sin(om_cw)*cos(ka_cw) + sin(ka_cw)*cos(om_cw))*(-2*sin(fi_cw)*cos(ka_cw)*cos(om_cw) + 2*sin(ka_cw)*sin(om_cw)) - 2*sin(om_cw)*pow(cos(fi_cw), 2)*cos(om_cw);
j.coeffRef(2,2) = (sin(fi_cw)*sin(ka_cw)*cos(om_cw) + sin(om_cw)*cos(ka_cw))*(2*sin(fi_cw)*sin(ka_cw)*cos(om_cw) + 2*sin(om_cw)*cos(ka_cw)) + (-2*sin(fi_cw)*cos(ka_cw)*cos(om_cw) + 2*sin(ka_cw)*sin(om_cw))*(-sin(fi_cw)*cos(ka_cw)*cos(om_cw) + sin(ka_cw)*sin(om_cw)) + 2*pow(cos(fi_cw), 2)*pow(cos(om_cw), 2);
}inline void point_to_point_source_to_target_tait_bryan_cw_d2sum_dbetadx(Eigen::Matrix<double, 3, 3, Eigen::RowMajor> &j, double tx_cw, double ty_cw, double tz_cw, double om_cw, double fi_cw, double ka_cw, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)
{j.coeffRef(0,0) = 2*cos(fi_cw)*cos(ka_cw);
j.coeffRef(0,1) = -2*sin(ka_cw)*cos(fi_cw);
j.coeffRef(0,2) = 2*sin(fi_cw);
j.coeffRef(1,0) = 2*sin(fi_cw)*sin(om_cw)*cos(ka_cw) + 2*sin(ka_cw)*cos(om_cw);
j.coeffRef(1,1) = -2*sin(fi_cw)*sin(ka_cw)*sin(om_cw) + 2*cos(ka_cw)*cos(om_cw);
j.coeffRef(1,2) = -2*sin(om_cw)*cos(fi_cw);
j.coeffRef(2,0) = -2*sin(fi_cw)*cos(ka_cw)*cos(om_cw) + 2*sin(ka_cw)*sin(om_cw);
j.coeffRef(2,1) = 2*sin(fi_cw)*sin(ka_cw)*cos(om_cw) + 2*sin(om_cw)*cos(ka_cw);
j.coeffRef(2,2) = 2*cos(fi_cw)*cos(om_cw);
}