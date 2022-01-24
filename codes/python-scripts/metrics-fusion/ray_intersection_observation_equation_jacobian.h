inline void ray_intersection_observation_equation(Eigen::Matrix<double, 2, 1> &delta, double t_x_int, double t_y_int, double t_z_int, double t_x_ln, double t_y_ln, double t_z_ln, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)
{delta.coeffRef(0,0) = -vxx*(t_x_int - t_x_ln) - vxy*(t_y_int - t_y_ln) - vxz*(t_z_int - t_z_ln);
delta.coeffRef(1,0) = -vyx*(t_x_int - t_x_ln) - vyy*(t_y_int - t_y_ln) - vyz*(t_z_int - t_z_ln);
}
inline void ray_intersection_observation_equation_jacobian(Eigen::Matrix<double, 2, 3> &j, double t_x_int, double t_y_int, double t_z_int, double t_x_ln, double t_y_ln, double t_z_ln, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)
{j.coeffRef(0,0) = -vxx;
j.coeffRef(0,1) = -vxy;
j.coeffRef(0,2) = -vxz;
j.coeffRef(1,0) = -vyx;
j.coeffRef(1,1) = -vyy;
j.coeffRef(1,2) = -vyz;
}