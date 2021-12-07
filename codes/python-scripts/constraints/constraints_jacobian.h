inline void observation_equation_constraint(double &delta, double a, double x, double x_trg)
{delta = -a*(-x + x_trg);
}
inline void observation_equation_constraint_jacobian(Eigen::Matrix<double, 1, 1> &j, double a, double x, double x_trg)
{j.coeffRef(0,0) = a;
}inline void observation_equation_sq_constraint(double &delta, double a, double x, double x_trg)
{delta = -pow(a, 2)*pow(-x + x_trg, 2);
}
inline void observation_equation_sq_constraint_jacobian(Eigen::Matrix<double, 1, 1> &j, double a, double x, double x_trg)
{j.coeffRef(0,0) = -pow(a, 2)*(2*x - 2*x_trg);
}