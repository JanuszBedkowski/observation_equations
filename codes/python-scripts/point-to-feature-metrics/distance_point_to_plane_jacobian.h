inline void delta_distance_point_to_plane(Eigen::Matrix<double, 1, 1> &delta, double x, double y, double z, double a, double b, double c, double d)
{delta.coeffRef(0,0) = -a*x - b*y - c*z - d;
}
inline void delta_distance_point_to_plane_jacobian(Eigen::Matrix<double, 1, 3> &j, double x, double y, double z, double a, double b, double c, double d)
{j.coeffRef(0,0) = -a;
j.coeffRef(0,1) = -b;
j.coeffRef(0,2) = -c;
}