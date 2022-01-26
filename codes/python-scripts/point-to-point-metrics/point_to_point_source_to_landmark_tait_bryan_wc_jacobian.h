inline void point_to_point_source_to_landmark_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_L, double y_L, double z_L)
{delta_x = -tx + x_L - x_s*cos(fi)*cos(ka) + y_s*sin(ka)*cos(fi) - z_s*sin(fi);
delta_y = -ty - x_s*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + y_L - y_s*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + z_s*sin(om)*cos(fi);
delta_z = -tz - x_s*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) - y_s*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + z_L - z_s*cos(fi)*cos(om);
}
inline void point_to_point_source_to_landmark_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 9, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s)
{j.coeffRef(0,0) = -1;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = 0;
j.coeffRef(0,4) = x_s*sin(fi)*cos(ka) - y_s*sin(fi)*sin(ka) - z_s*cos(fi);
j.coeffRef(0,5) = x_s*sin(ka)*cos(fi) + y_s*cos(fi)*cos(ka);
j.coeffRef(0,6) = 1;
j.coeffRef(0,7) = 0;
j.coeffRef(0,8) = 0;
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = -1;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = -x_s*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) - y_s*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) + z_s*cos(fi)*cos(om);
j.coeffRef(1,4) = -x_s*sin(om)*cos(fi)*cos(ka) + y_s*sin(ka)*sin(om)*cos(fi) - z_s*sin(fi)*sin(om);
j.coeffRef(1,5) = -x_s*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - y_s*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om));
j.coeffRef(1,6) = 0;
j.coeffRef(1,7) = 1;
j.coeffRef(1,8) = 0;
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = -1;
j.coeffRef(2,3) = -x_s*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) - y_s*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + z_s*sin(om)*cos(fi);
j.coeffRef(2,4) = x_s*cos(fi)*cos(ka)*cos(om) - y_s*sin(ka)*cos(fi)*cos(om) + z_s*sin(fi)*cos(om);
j.coeffRef(2,5) = -x_s*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) - y_s*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om));
j.coeffRef(2,6) = 0;
j.coeffRef(2,7) = 0;
j.coeffRef(2,8) = 1;
}