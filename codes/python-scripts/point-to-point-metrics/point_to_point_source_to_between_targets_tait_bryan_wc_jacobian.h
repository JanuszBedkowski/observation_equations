inline void point_to_point_source_to_between_targets_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t1, double y_t1, double z_t1, double x_t2, double y_t2, double z_t2)
{delta_x = -2*tx - 2*x_s*cos(fi)*cos(ka) + x_t1 + x_t2 + 2*y_s*sin(ka)*cos(fi) - 2*z_s*sin(fi);
delta_y = -2*ty - 2*x_s*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) - 2*y_s*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + y_t1 + y_t2 + 2*z_s*sin(om)*cos(fi);
delta_z = -2*tz - 2*x_s*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) - 2*y_s*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) - 2*z_s*cos(fi)*cos(om) + z_t1 + z_t2;
}
inline void point_to_point_source_to_between_targets_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s)
{j.coeffRef(0,0) = -2;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = 0;
j.coeffRef(0,4) = 2*x_s*sin(fi)*cos(ka) - 2*y_s*sin(fi)*sin(ka) - 2*z_s*cos(fi);
j.coeffRef(0,5) = 2*x_s*sin(ka)*cos(fi) + 2*y_s*cos(fi)*cos(ka);
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = -2;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = -2*x_s*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) - 2*y_s*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) + 2*z_s*cos(fi)*cos(om);
j.coeffRef(1,4) = -2*x_s*sin(om)*cos(fi)*cos(ka) + 2*y_s*sin(ka)*sin(om)*cos(fi) - 2*z_s*sin(fi)*sin(om);
j.coeffRef(1,5) = -2*x_s*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - 2*y_s*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om));
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = -2;
j.coeffRef(2,3) = -2*x_s*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) - 2*y_s*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + 2*z_s*sin(om)*cos(fi);
j.coeffRef(2,4) = 2*x_s*cos(fi)*cos(ka)*cos(om) - 2*y_s*sin(ka)*cos(fi)*cos(om) + 2*z_s*sin(fi)*cos(om);
j.coeffRef(2,5) = -2*x_s*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) - 2*y_s*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om));
}