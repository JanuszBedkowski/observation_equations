inline void sheaf_of_planes_observation_equation(Eigen::Matrix<double, 4, 1> &residual, double a_x, double a_y, double a_z, double b_x, double b_y, double b_z, double a, double b, double c, double d)
{residual.coeffRef(0,0) = -a*a_x - a_y*b - a_z*c - d;
residual.coeffRef(1,0) = -a*b_x - b*b_y - b_z*c - d;
residual.coeffRef(2,0) = -a_x*(-a_x + b_x) - a_y*(-a_y + b_y) - a_z*(-a_z + b_z);
residual.coeffRef(3,0) = -a*(-a_x + b_x) - b*(-a_y + b_y) - c*(-a_z + b_z);
}
inline void sheaf_of_planes_observation_equation_jacobian(Eigen::Matrix<double, 4, 6> &j, double a_x, double a_y, double a_z, double b_x, double b_y, double b_z, double a, double b, double c, double d)
{j.coeffRef(0,0) = -a;
j.coeffRef(0,1) = -b;
j.coeffRef(0,2) = -c;
j.coeffRef(0,3) = 0;
j.coeffRef(0,4) = 0;
j.coeffRef(0,5) = 0;
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = 0;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = -a;
j.coeffRef(1,4) = -b;
j.coeffRef(1,5) = -c;
j.coeffRef(2,0) = 2*a_x - b_x;
j.coeffRef(2,1) = 2*a_y - b_y;
j.coeffRef(2,2) = 2*a_z - b_z;
j.coeffRef(2,3) = -a_x;
j.coeffRef(2,4) = -a_y;
j.coeffRef(2,5) = -a_z;
j.coeffRef(3,0) = a;
j.coeffRef(3,1) = b;
j.coeffRef(3,2) = c;
j.coeffRef(3,3) = -a;
j.coeffRef(3,4) = -b;
j.coeffRef(3,5) = -c;
}inline void line_direction_norm_observation_equation(Eigen::Matrix<double, 1, 1> &residual, double a_x, double a_y, double a_z, double b_x, double b_y, double b_z)
{residual.coeffRef(0,0) = 1 - sqrt(pow(-a_x + b_x, 2) + pow(-a_y + b_y, 2) + pow(-a_z + b_z, 2));
}
inline void line_direction_observation_equation_jacobian(Eigen::Matrix<double, 1, 6> &j, double a_x, double a_y, double a_z, double b_x, double b_y, double b_z)
{j.coeffRef(0,0) = -(a_x - b_x)/sqrt(pow(-a_x + b_x, 2) + pow(-a_y + b_y, 2) + pow(-a_z + b_z, 2));
j.coeffRef(0,1) = -(a_y - b_y)/sqrt(pow(-a_x + b_x, 2) + pow(-a_y + b_y, 2) + pow(-a_z + b_z, 2));
j.coeffRef(0,2) = -(a_z - b_z)/sqrt(pow(-a_x + b_x, 2) + pow(-a_y + b_y, 2) + pow(-a_z + b_z, 2));
j.coeffRef(0,3) = -(-a_x + b_x)/sqrt(pow(-a_x + b_x, 2) + pow(-a_y + b_y, 2) + pow(-a_z + b_z, 2));
j.coeffRef(0,4) = -(-a_y + b_y)/sqrt(pow(-a_x + b_x, 2) + pow(-a_y + b_y, 2) + pow(-a_z + b_z, 2));
j.coeffRef(0,5) = -(-a_z + b_z)/sqrt(pow(-a_x + b_x, 2) + pow(-a_y + b_y, 2) + pow(-a_z + b_z, 2));
}