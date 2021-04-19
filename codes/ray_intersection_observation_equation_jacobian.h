inline void ray_intersection_observation_equation(Eigen::Matrix<double, 2, 1> &delta, double T_x_int, double T_y_int, double T_z_int, double T_x_ln, double T_y_ln, double T_z_ln, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)
{delta.coeffRef(0,0) = -vxx*(T_x_int - T_x_ln) - vxy*(T_y_int - T_y_ln) - vxz*(T_z_int - T_z_ln);
delta.coeffRef(1,0) = -vyx*(T_x_int - T_x_ln) - vyy*(T_y_int - T_y_ln) - vyz*(T_z_int - T_z_ln);
}
inline void ray_intersection_observation_equation_jacobian(Eigen::Matrix<double, 2, 3> &j, double T_x_int, double T_y_int, double T_z_int, double T_x_ln, double T_y_ln, double T_z_ln, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)
{j.coeffRef(0,0) = -vxx;
j.coeffRef(0,1) = -vxy;
j.coeffRef(0,2) = -vxz;
j.coeffRef(1,0) = -vyx;
j.coeffRef(1,1) = -vyy;
j.coeffRef(1,2) = -vyz;
}