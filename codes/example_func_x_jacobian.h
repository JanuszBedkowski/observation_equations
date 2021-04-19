inline void example_func_x(double &psi, double x)
{psi = 0.10000000000000001*pow(x, 3) + 0.5*pow(x, 2);
}
inline void observation_equation_example_func_x(double &delta, double x, double y)
{delta = 0.10000000000000001*pow(x, 3) + 0.5*pow(x, 2) - y;
}
inline void observation_equation_example_func_x_jacobian(Eigen::Matrix<double, 1, 1> &j, double x)
{j.coeffRef(0,0) = 0.30000000000000004*pow(x, 2) + 1.0*x;
}