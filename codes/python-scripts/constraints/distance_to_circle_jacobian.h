inline void observation_equation_distance_to_circle(double &delta, double x, double y, double cx, double cy, double cr)
{delta = cr - sqrt(pow(-cx + x, 2) + pow(-cy + y, 2));
}
inline void observation_equation_distance_to_circle_jacobian(Eigen::Matrix<double, 1, 3, Eigen::RowMajor> &j, double x, double y, double cx, double cy, double cr)
{j.coeffRef(0,0) = -(cx - x)/sqrt(pow(-cx + x, 2) + pow(-cy + y, 2));
j.coeffRef(0,1) = -(cy - y)/sqrt(pow(-cx + x, 2) + pow(-cy + y, 2));
j.coeffRef(0,2) = 1;
}