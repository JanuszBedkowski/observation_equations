inline void quaternion_constraint(double &delta, double q0, double q1, double q2, double q3)
{delta = 1 - sqrt(pow(q0, 2) + pow(q1, 2) + pow(q2, 2) + pow(q3, 2));
}
inline void quaternion_constraint_jacobian(Eigen::Matrix<double, 1, 4> &j, double q0, double q1, double q2, double q3)
{j.coeffRef(0,0) = -q0/sqrt(pow(q0, 2) + pow(q1, 2) + pow(q2, 2) + pow(q3, 2));
j.coeffRef(0,1) = -q1/sqrt(pow(q0, 2) + pow(q1, 2) + pow(q2, 2) + pow(q3, 2));
j.coeffRef(0,2) = -q2/sqrt(pow(q0, 2) + pow(q1, 2) + pow(q2, 2) + pow(q3, 2));
j.coeffRef(0,3) = -q3/sqrt(pow(q0, 2) + pow(q1, 2) + pow(q2, 2) + pow(q3, 2));
}