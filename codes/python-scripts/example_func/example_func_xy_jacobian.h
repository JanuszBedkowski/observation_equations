inline void example_func_xy(double &psi, double x, double y)
{psi = exp(-1.0/50.0*pow(x - 1, 2) - 1.0/50.0*pow(y - 2, 2));
}
inline void observation_equation_example_func_xy(double &delta, double x, double y, double z)
{delta = -z + exp(-1.0/50.0*pow(x - 1, 2) - 1.0/50.0*pow(y - 2, 2));
}
inline void observation_equation_example_func_xy_jacobian(Eigen::Matrix<double, 1, 2> &j, double x, double y)
{j.coeffRef(0,0) = (1.0/25.0 - 1.0/25.0*x)*exp(-1.0/50.0*pow(x - 1, 2) - 1.0/50.0*pow(y - 2, 2));
j.coeffRef(0,1) = (2.0/25.0 - 1.0/25.0*y)*exp(-1.0/50.0*pow(x - 1, 2) - 1.0/50.0*pow(y - 2, 2));
}