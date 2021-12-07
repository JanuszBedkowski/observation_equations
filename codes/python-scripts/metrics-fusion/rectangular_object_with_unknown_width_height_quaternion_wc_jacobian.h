inline void rectangular_object_with_unknown_width_height_quaternion_wc(Eigen::Matrix<double, 2, 1> &delta, double px, double py, double pz, double q0, double q1, double q2, double q3, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz, double scale_x, double scale_y)
{delta.coeffRef(0,0) = -vxx*(px + scale_x*xsl*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) + scale_y*ysl*(-2.0*q0*q3 + 2.0*q1*q2) - xtg + zsl*(2.0*q0*q2 + 2.0*q1*q3)) - vxy*(py + scale_x*xsl*(2.0*q0*q3 + 2.0*q1*q2) + scale_y*ysl*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) - ytg + zsl*(-2.0*q0*q1 + 2.0*q2*q3)) - vxz*(pz + scale_x*xsl*(-2.0*q0*q2 + 2.0*q1*q3) + scale_y*ysl*(2.0*q0*q1 + 2.0*q2*q3) + zsl*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1) - ztg);
delta.coeffRef(1,0) = -vyx*(px + scale_x*xsl*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) + scale_y*ysl*(-2.0*q0*q3 + 2.0*q1*q2) - xtg + zsl*(2.0*q0*q2 + 2.0*q1*q3)) - vyy*(py + scale_x*xsl*(2.0*q0*q3 + 2.0*q1*q2) + scale_y*ysl*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) - ytg + zsl*(-2.0*q0*q1 + 2.0*q2*q3)) - vyz*(pz + scale_x*xsl*(-2.0*q0*q2 + 2.0*q1*q3) + scale_y*ysl*(2.0*q0*q1 + 2.0*q2*q3) + zsl*(-2*pow(q1, 2) - 2*pow(q2, 2) + 1) - ztg);
}
inline void rectangular_object_with_unknown_width_height_quaternion_wc_jacobian(Eigen::Matrix<double, 2, 9> &j, double px, double py, double pz, double q0, double q1, double q2, double q3, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz, double scale_x, double scale_y)
{j.coeffRef(0,0) = -vxx;
j.coeffRef(0,1) = -vxy;
j.coeffRef(0,2) = -vxz;
j.coeffRef(0,3) = -vxx*(2.0*q2*zsl - 2.0*q3*scale_y*ysl) - vxy*(-2.0*q1*zsl + 2.0*q3*scale_x*xsl) - vxz*(2.0*q1*scale_y*ysl - 2.0*q2*scale_x*xsl);
j.coeffRef(0,4) = -vxx*(2.0*q2*scale_y*ysl + 2.0*q3*zsl) - vxy*(-2.0*q0*zsl - 4*q1*scale_y*ysl + 2.0*q2*scale_x*xsl) - vxz*(2.0*q0*scale_y*ysl - 4*q1*zsl + 2.0*q3*scale_x*xsl);
j.coeffRef(0,5) = -vxx*(2.0*q0*zsl + 2.0*q1*scale_y*ysl - 4*q2*scale_x*xsl) - vxy*(2.0*q1*scale_x*xsl + 2.0*q3*zsl) - vxz*(-2.0*q0*scale_x*xsl - 4*q2*zsl + 2.0*q3*scale_y*ysl);
j.coeffRef(0,6) = -vxx*(-2.0*q0*scale_y*ysl + 2.0*q1*zsl - 4*q3*scale_x*xsl) - vxy*(2.0*q0*scale_x*xsl + 2.0*q2*zsl - 4*q3*scale_y*ysl) - vxz*(2.0*q1*scale_x*xsl + 2.0*q2*scale_y*ysl);
j.coeffRef(0,7) = -vxx*xsl*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) - vxy*xsl*(2.0*q0*q3 + 2.0*q1*q2) - vxz*xsl*(-2.0*q0*q2 + 2.0*q1*q3);
j.coeffRef(0,8) = -vxx*ysl*(-2.0*q0*q3 + 2.0*q1*q2) - vxy*ysl*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) - vxz*ysl*(2.0*q0*q1 + 2.0*q2*q3);
j.coeffRef(1,0) = -vyx;
j.coeffRef(1,1) = -vyy;
j.coeffRef(1,2) = -vyz;
j.coeffRef(1,3) = -vyx*(2.0*q2*zsl - 2.0*q3*scale_y*ysl) - vyy*(-2.0*q1*zsl + 2.0*q3*scale_x*xsl) - vyz*(2.0*q1*scale_y*ysl - 2.0*q2*scale_x*xsl);
j.coeffRef(1,4) = -vyx*(2.0*q2*scale_y*ysl + 2.0*q3*zsl) - vyy*(-2.0*q0*zsl - 4*q1*scale_y*ysl + 2.0*q2*scale_x*xsl) - vyz*(2.0*q0*scale_y*ysl - 4*q1*zsl + 2.0*q3*scale_x*xsl);
j.coeffRef(1,5) = -vyx*(2.0*q0*zsl + 2.0*q1*scale_y*ysl - 4*q2*scale_x*xsl) - vyy*(2.0*q1*scale_x*xsl + 2.0*q3*zsl) - vyz*(-2.0*q0*scale_x*xsl - 4*q2*zsl + 2.0*q3*scale_y*ysl);
j.coeffRef(1,6) = -vyx*(-2.0*q0*scale_y*ysl + 2.0*q1*zsl - 4*q3*scale_x*xsl) - vyy*(2.0*q0*scale_x*xsl + 2.0*q2*zsl - 4*q3*scale_y*ysl) - vyz*(2.0*q1*scale_x*xsl + 2.0*q2*scale_y*ysl);
j.coeffRef(1,7) = -vyx*xsl*(-2*pow(q2, 2) - 2*pow(q3, 2) + 1) - vyy*xsl*(2.0*q0*q3 + 2.0*q1*q2) - vyz*xsl*(-2.0*q0*q2 + 2.0*q1*q3);
j.coeffRef(1,8) = -vyx*ysl*(-2.0*q0*q3 + 2.0*q1*q2) - vyy*ysl*(-2*pow(q1, 2) - 2*pow(q3, 2) + 1) - vyz*ysl*(2.0*q0*q1 + 2.0*q2*q3);
}