inline void point_to_point_with_scale_source_to_target_tait_bryan_wc(double &delta_x, double &delta_y, double &delta_z, [[maybe_unused]] double tx, [[maybe_unused]] double ty, [[maybe_unused]] double tz, [[maybe_unused]] double om, [[maybe_unused]] double fi, [[maybe_unused]] double ka, [[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] double x_target, [[maybe_unused]] double y_target, [[maybe_unused]] double z_target, [[maybe_unused]] double s)
{delta_x = -s*x*cos(fi)*cos(ka) + s*y*sin(ka)*cos(fi) - s*z*sin(fi) - tx + x_target;
delta_y = -s*x*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) - s*y*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + s*z*sin(om)*cos(fi) - ty + y_target;
delta_z = -s*x*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) - s*y*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) - s*z*cos(fi)*cos(om) - tz + z_target;
}
inline void point_to_point_with_scale_source_to_target_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 7, Eigen::RowMajor> &j, [[maybe_unused]] double tx, [[maybe_unused]] double ty, [[maybe_unused]] double tz, [[maybe_unused]] double om, [[maybe_unused]] double fi, [[maybe_unused]] double ka, [[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] double x_target, [[maybe_unused]] double y_target, [[maybe_unused]] double z_target, [[maybe_unused]] double s)
{j.coeffRef(0,0) = -1;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = 0;
j.coeffRef(0,4) = s*x*sin(fi)*cos(ka) - s*y*sin(fi)*sin(ka) - s*z*cos(fi);
j.coeffRef(0,5) = s*x*sin(ka)*cos(fi) + s*y*cos(fi)*cos(ka);
j.coeffRef(0,6) = -x*cos(fi)*cos(ka) + y*sin(ka)*cos(fi) - z*sin(fi);
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = -1;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = -s*x*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) - s*y*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) + s*z*cos(fi)*cos(om);
j.coeffRef(1,4) = -s*x*sin(om)*cos(fi)*cos(ka) + s*y*sin(ka)*sin(om)*cos(fi) - s*z*sin(fi)*sin(om);
j.coeffRef(1,5) = -s*x*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - s*y*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om));
j.coeffRef(1,6) = -x*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) - y*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + z*sin(om)*cos(fi);
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = -1;
j.coeffRef(2,3) = -s*x*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) - s*y*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + s*z*sin(om)*cos(fi);
j.coeffRef(2,4) = s*x*cos(fi)*cos(ka)*cos(om) - s*y*sin(ka)*cos(fi)*cos(om) + s*z*sin(fi)*cos(om);
j.coeffRef(2,5) = -s*x*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) - s*y*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om));
j.coeffRef(2,6) = -x*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) - y*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) - z*cos(fi)*cos(om);
}