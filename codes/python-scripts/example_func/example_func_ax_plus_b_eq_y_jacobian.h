inline void example_func_ax_plus_b(double &y, double x, double a, double b)
{y = a*x + b;
}
inline void observation_equation_example_func_ax_plus_b_eq_y(double &delta, double x, double y, double a, double b)
{delta = -a*x - b + y;
}
inline void observation_equation_example_func_ax_plus_b_eq_y_jacobian(Eigen::Matrix<double, 1, 2> &j, double x, double y, double a, double b)
{j.coeffRef(0,0) = -x;
j.coeffRef(0,1) = -1;
}