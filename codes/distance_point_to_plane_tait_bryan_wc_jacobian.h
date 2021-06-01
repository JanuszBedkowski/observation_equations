inline void delta_distance_point_to_plane_tait_bryan_wc(Eigen::Matrix<double, 1, 1> &delta, double px, double py, double pz, double om, double fi, double ka, double x, double y, double z, double a, double b, double c, double d)
{delta.coeffRef(0,0) = -a*(px + x*cos(fi)*cos(ka) - y*sin(ka)*cos(fi) + z*sin(fi)) - b*(py + x*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + y*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - z*sin(om)*cos(fi)) - c*(pz + x*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) + y*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + z*cos(fi)*cos(om)) - d;
}
inline void delta_distance_point_to_plane_tait_bryan_wc_jacobian(Eigen::Matrix<double, 1, 6> &j, double px, double py, double pz, double om, double fi, double ka, double x, double y, double z, double a, double b, double c, double d)
{j.coeffRef(0,0) = -a;
j.coeffRef(0,1) = -b;
j.coeffRef(0,2) = -c;
j.coeffRef(0,3) = -b*(x*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + y*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) - z*cos(fi)*cos(om)) - c*(x*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + y*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - z*sin(om)*cos(fi));
j.coeffRef(0,4) = -a*(-x*sin(fi)*cos(ka) + y*sin(fi)*sin(ka) + z*cos(fi)) - b*(x*sin(om)*cos(fi)*cos(ka) - y*sin(ka)*sin(om)*cos(fi) + z*sin(fi)*sin(om)) - c*(-x*cos(fi)*cos(ka)*cos(om) + y*sin(ka)*cos(fi)*cos(om) - z*sin(fi)*cos(om));
j.coeffRef(0,5) = -a*(-x*sin(ka)*cos(fi) - y*cos(fi)*cos(ka)) - b*(x*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + y*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om))) - c*(x*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + y*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)));
}