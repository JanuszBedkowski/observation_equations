#ifndef _distance_point_to_plane_tait_bryan_cw_jacobian_h_
#define _distance_point_to_plane_tait_bryan_cw_jacobian_h_
inline void delta_distance_point_to_plane_tait_bryan_cw(Eigen::Matrix<double, 1, 1> &delta, double tx, double ty, double tz, double om, double fi, double ka, double x, double y, double z, double a, double b, double c, double d)
{delta.coeffRef(0,0) = -a*(-tx*cos(fi)*cos(ka) + ty*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om)) + tz*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + x*cos(fi)*cos(ka) + y*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + z*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om))) - b*(tx*sin(ka)*cos(fi) + ty*(sin(fi)*sin(ka)*sin(om) - cos(ka)*cos(om)) + tz*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) - x*sin(ka)*cos(fi) + y*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + z*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka))) - c*(-tx*sin(fi) + ty*sin(om)*cos(fi) - tz*cos(fi)*cos(om) + x*sin(fi) - y*sin(om)*cos(fi) + z*cos(fi)*cos(om)) - d;
}
inline void delta_distance_point_to_plane_tait_bryan_cw_jacobian(Eigen::Matrix<double, 1, 6> &j, double tx, double ty, double tz, double om, double fi, double ka, double x, double y, double z, double a, double b, double c, double d)
{j.coeffRef(0,0) = a*cos(fi)*cos(ka) - b*sin(ka)*cos(fi) + c*sin(fi);
j.coeffRef(0,1) = -a*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om)) - b*(sin(fi)*sin(ka)*sin(om) - cos(ka)*cos(om)) - c*sin(om)*cos(fi);
j.coeffRef(0,2) = -a*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) - b*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) + c*cos(fi)*cos(om);
j.coeffRef(0,3) = -a*(ty*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) + tz*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om)) + y*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + z*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om))) - b*(ty*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + tz*(sin(fi)*sin(ka)*sin(om) - cos(ka)*cos(om)) + y*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) + z*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om))) - c*(ty*cos(fi)*cos(om) + tz*sin(om)*cos(fi) - y*cos(fi)*cos(om) - z*sin(om)*cos(fi));
j.coeffRef(0,4) = -a*(tx*sin(fi)*cos(ka) - ty*sin(om)*cos(fi)*cos(ka) + tz*cos(fi)*cos(ka)*cos(om) - x*sin(fi)*cos(ka) + y*sin(om)*cos(fi)*cos(ka) - z*cos(fi)*cos(ka)*cos(om)) - b*(-tx*sin(fi)*sin(ka) + ty*sin(ka)*sin(om)*cos(fi) - tz*sin(ka)*cos(fi)*cos(om) + x*sin(fi)*sin(ka) - y*sin(ka)*sin(om)*cos(fi) + z*sin(ka)*cos(fi)*cos(om)) - c*(-tx*cos(fi) - ty*sin(fi)*sin(om) + tz*sin(fi)*cos(om) + x*cos(fi) + y*sin(fi)*sin(om) - z*sin(fi)*cos(om));
j.coeffRef(0,5) = -a*(tx*sin(ka)*cos(fi) + ty*(sin(fi)*sin(ka)*sin(om) - cos(ka)*cos(om)) + tz*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) - x*sin(ka)*cos(fi) + y*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + z*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka))) - b*(tx*cos(fi)*cos(ka) + ty*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + tz*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) - x*cos(fi)*cos(ka) + y*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om)) + z*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)));
}
#endif