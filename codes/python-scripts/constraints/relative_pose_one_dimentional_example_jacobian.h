inline void relative_pose_obs_eq(Eigen::Matrix<double, 1, 1> &delta, double pose_from, double pose_to, double pose_m)
{delta.coeffRef(0,0) = pose_m - pose_to/pose_from;
}
inline void relative_pose_obs_eq_jacobian(Eigen::Matrix<double, 1, 2, Eigen::RowMajor> &j, double pose_from, double pose_to)
{j.coeffRef(0,0) = pose_to/pow(pose_from, 2);
j.coeffRef(0,1) = -1/pose_from;
}inline void relative_pose(Eigen::Matrix<double, 1, 1> &relative_pose, double pose_from, double pose_to)
{relative_pose.coeffRef(0,0) = pose_to/pose_from;
}
