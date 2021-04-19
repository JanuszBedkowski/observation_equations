inline void point_to_line_quaternion_wc(Eigen::Matrix<double, 2, 1> &delta, double px, double py, double pz, double q0, double q1, double q2, double q3, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)
{delta.coeffRef(0,0) = -vxx*(px + xsl*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) - xtg + ysl*(-2.0*q0*q3 + 2.0*q1*q2) + zsl*(2.0*q0*q2 + 2.0*q1*q3)) - vxy*(py + xsl*(2.0*q0*q3 + 2.0*q1*q2) + ysl*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) - ytg + zsl*(-2.0*q0*q1 + 2.0*q2*q3)) - vxz*(pz + xsl*(-2.0*q0*q2 + 2.0*q1*q3) + ysl*(2.0*q0*q1 + 2.0*q2*q3) + zsl*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1) - ztg);
delta.coeffRef(1,0) = -vyx*(px + xsl*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) - xtg + ysl*(-2.0*q0*q3 + 2.0*q1*q2) + zsl*(2.0*q0*q2 + 2.0*q1*q3)) - vyy*(py + xsl*(2.0*q0*q3 + 2.0*q1*q2) + ysl*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) - ytg + zsl*(-2.0*q0*q1 + 2.0*q2*q3)) - vyz*(pz + xsl*(-2.0*q0*q2 + 2.0*q1*q3) + ysl*(2.0*q0*q1 + 2.0*q2*q3) + zsl*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1) - ztg);
}
inline void point_to_line_quaternion_wc_jacobian(Eigen::Matrix<double, 2, 7> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)
{j.coeffRef(0,0) = -vxx;
j.coeffRef(0,1) = -vxy;
j.coeffRef(0,2) = -vxz;
j.coeffRef(0,3) = -vxx*(2.0*q2*zsl - 2.0*q3*ysl) - vxy*(-2.0*q1*zsl + 2.0*q3*xsl) - vxz*(2.0*q1*ysl - 2.0*q2*xsl);
j.coeffRef(0,4) = -vxx*(2.0*q2*ysl + 2.0*q3*zsl) - vxy*(-2.0*q0*zsl - 4*q1*ysl + 2.0*q2*xsl) - vxz*(2.0*q0*ysl - 4*q1*zsl + 2.0*q3*xsl);
j.coeffRef(0,5) = -vxx*(2.0*q0*zsl + 2.0*q1*ysl - 4*q2*xsl) - vxy*(2.0*q1*xsl + 2.0*q3*zsl) - vxz*(-2.0*q0*xsl - 4*q2*zsl + 2.0*q3*ysl);
j.coeffRef(0,6) = -vxx*(-2.0*q0*ysl + 2.0*q1*zsl - 4*q3*xsl) - vxy*(2.0*q0*xsl + 2.0*q2*zsl - 4*q3*ysl) - vxz*(2.0*q1*xsl + 2.0*q2*ysl);
j.coeffRef(1,0) = -vyx;
j.coeffRef(1,1) = -vyy;
j.coeffRef(1,2) = -vyz;
j.coeffRef(1,3) = -vyx*(2.0*q2*zsl - 2.0*q3*ysl) - vyy*(-2.0*q1*zsl + 2.0*q3*xsl) - vyz*(2.0*q1*ysl - 2.0*q2*xsl);
j.coeffRef(1,4) = -vyx*(2.0*q2*ysl + 2.0*q3*zsl) - vyy*(-2.0*q0*zsl - 4*q1*ysl + 2.0*q2*xsl) - vyz*(2.0*q0*ysl - 4*q1*zsl + 2.0*q3*xsl);
j.coeffRef(1,5) = -vyx*(2.0*q0*zsl + 2.0*q1*ysl - 4*q2*xsl) - vyy*(2.0*q1*xsl + 2.0*q3*zsl) - vyz*(-2.0*q0*xsl - 4*q2*zsl + 2.0*q3*ysl);
j.coeffRef(1,6) = -vyx*(-2.0*q0*ysl + 2.0*q1*zsl - 4*q3*xsl) - vyy*(2.0*q0*xsl + 2.0*q2*zsl - 4*q3*ysl) - vyz*(2.0*q1*xsl + 2.0*q2*ysl);
}