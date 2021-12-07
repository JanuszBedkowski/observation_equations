inline void residual_constraint_fixed_optimization_parameter(double &residual, double target_value, double model_function)
{residual = -model_function + target_value;
}
inline void residual_constraint_fixed_optimization_parameter_jacobian(Eigen::Matrix<double, 1, 1> &j, double target_value, double model_function)
{j.coeffRef(0,0) = -1;
}