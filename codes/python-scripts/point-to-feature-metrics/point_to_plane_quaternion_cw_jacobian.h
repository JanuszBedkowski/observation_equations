#ifndef _point_to_plane_quaternion_cw_jacobian_h_
#define _point_to_plane_quaternion_cw_jacobian_h_
inline void point_to_plane_quaternion_cw(Eigen::Matrix<double, 1, 1> &delta, double px, double py, double pz, double q0, double q1, double q2, double q3, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vzx, double vzy, double vzz)
{delta.coeffRef(0,0) = -vzx*(px*(2*pow(q2, 2) + 2*pow(q3, 2) - 1) + py*(-2.0*q0*q3 - 2.0*q1*q2) + pz*(2.0*q0*q2 - 2.0*q1*q3) + xsl*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) - xtg + ysl*(2.0*q0*q3 + 2.0*q1*q2) + zsl*(-2.0*q0*q2 + 2.0*q1*q3)) - vzy*(px*(2.0*q0*q3 - 2.0*q1*q2) + py*(2*pow(q1, 2) + 2*pow(q3, 2) - 1) + pz*(-2.0*q0*q1 - 2.0*q2*q3) + xsl*(-2.0*q0*q3 + 2.0*q1*q2) + ysl*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) - ytg + zsl*(2.0*q0*q1 + 2.0*q2*q3)) - vzz*(px*(-2.0*q0*q2 - 2.0*q1*q3) + py*(2.0*q0*q1 - 2.0*q2*q3) + pz*(2*pow(q1, 2) + 2*pow(q2, 2) - 1) + xsl*(2.0*q0*q2 + 2.0*q1*q3) + ysl*(-2.0*q0*q1 + 2.0*q2*q3) + zsl*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1) - ztg);
}
inline void point_to_plane_quaternion_cw_jacobian(Eigen::Matrix<double, 1, 7> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vzx, double vzy, double vzz)
{j.coeffRef(0,0) = -vzx*(2*pow(q2, 2) + 2*pow(q3, 2) - 1) - vzy*(2.0*q0*q3 - 2.0*q1*q2) - vzz*(-2.0*q0*q2 - 2.0*q1*q3);
j.coeffRef(0,1) = -vzx*(-2.0*q0*q3 - 2.0*q1*q2) - vzy*(2*pow(q1, 2) + 2*pow(q3, 2) - 1) - vzz*(2.0*q0*q1 - 2.0*q2*q3);
j.coeffRef(0,2) = -vzx*(2.0*q0*q2 - 2.0*q1*q3) - vzy*(-2.0*q0*q1 - 2.0*q2*q3) - vzz*(2*pow(q1, 2) + 2*pow(q2, 2) - 1);
j.coeffRef(0,3) = -vzx*(-2.0*py*q3 + 2.0*pz*q2 - 2.0*q2*zsl + 2.0*q3*ysl) - vzy*(2.0*px*q3 - 2.0*pz*q1 + 2.0*q1*zsl - 2.0*q3*xsl) - vzz*(-2.0*px*q2 + 2.0*py*q1 - 2.0*q1*ysl + 2.0*q2*xsl);
j.coeffRef(0,4) = -vzx*(-2.0*py*q2 - 2.0*pz*q3 + 2.0*q2*ysl + 2.0*q3*zsl) - vzy*(-2.0*px*q2 + 4*py*q1 - 2.0*pz*q0 + 2.0*q0*zsl - 4*q1*ysl + 2.0*q2*xsl) - vzz*(-2.0*px*q3 + 2.0*py*q0 + 4*pz*q1 - 2.0*q0*ysl - 4*q1*zsl + 2.0*q3*xsl);
j.coeffRef(0,5) = -vzx*(4*px*q2 - 2.0*py*q1 + 2.0*pz*q0 - 2.0*q0*zsl + 2.0*q1*ysl - 4*q2*xsl) - vzy*(-2.0*px*q1 - 2.0*pz*q3 + 2.0*q1*xsl + 2.0*q3*zsl) - vzz*(-2.0*px*q0 - 2.0*py*q3 + 4*pz*q2 + 2.0*q0*xsl - 4*q2*zsl + 2.0*q3*ysl);
j.coeffRef(0,6) = -vzx*(4*px*q3 - 2.0*py*q0 - 2.0*pz*q1 + 2.0*q0*ysl + 2.0*q1*zsl - 4*q3*xsl) - vzy*(2.0*px*q0 + 4*py*q3 - 2.0*pz*q2 - 2.0*q0*xsl + 2.0*q2*zsl - 4*q3*ysl) - vzz*(-2.0*px*q1 - 2.0*py*q2 + 2.0*q1*xsl + 2.0*q2*ysl);
}#endif