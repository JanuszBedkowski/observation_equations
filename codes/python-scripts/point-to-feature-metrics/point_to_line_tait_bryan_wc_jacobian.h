inline void point_to_line_tait_bryan_wc(Eigen::Matrix<double, 2, 1> &delta, double tx, double ty, double tz, double om, double fi, double ka, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)
{delta.coeffRef(0,0) = -vxx*(tx + xsl*cos(fi)*cos(ka) - xtg - ysl*sin(ka)*cos(fi) + zsl*sin(fi)) - vxy*(ty + xsl*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + ysl*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - ytg - zsl*sin(om)*cos(fi)) - vxz*(tz + xsl*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) + ysl*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + zsl*cos(fi)*cos(om) - ztg);
delta.coeffRef(1,0) = -vyx*(tx + xsl*cos(fi)*cos(ka) - xtg - ysl*sin(ka)*cos(fi) + zsl*sin(fi)) - vyy*(ty + xsl*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + ysl*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - ytg - zsl*sin(om)*cos(fi)) - vyz*(tz + xsl*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) + ysl*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + zsl*cos(fi)*cos(om) - ztg);
}
inline void point_to_line_tait_bryan_wc_jacobian(Eigen::Matrix<double, 2, 6> &j, double tx, double ty, double tz, double om, double fi, double ka, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vxx, double vxy, double vxz, double vyx, double vyy, double vyz)
{j.coeffRef(0,0) = -vxx;
j.coeffRef(0,1) = -vxy;
j.coeffRef(0,2) = -vxz;
j.coeffRef(0,3) = -vxy*(xsl*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + ysl*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) - zsl*cos(fi)*cos(om)) - vxz*(xsl*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + ysl*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - zsl*sin(om)*cos(fi));
j.coeffRef(0,4) = -vxx*(-xsl*sin(fi)*cos(ka) + ysl*sin(fi)*sin(ka) + zsl*cos(fi)) - vxy*(xsl*sin(om)*cos(fi)*cos(ka) - ysl*sin(ka)*sin(om)*cos(fi) + zsl*sin(fi)*sin(om)) - vxz*(-xsl*cos(fi)*cos(ka)*cos(om) + ysl*sin(ka)*cos(fi)*cos(om) - zsl*sin(fi)*cos(om));
j.coeffRef(0,5) = -vxx*(-xsl*sin(ka)*cos(fi) - ysl*cos(fi)*cos(ka)) - vxy*(xsl*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + ysl*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om))) - vxz*(xsl*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + ysl*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)));
j.coeffRef(1,0) = -vyx;
j.coeffRef(1,1) = -vyy;
j.coeffRef(1,2) = -vyz;
j.coeffRef(1,3) = -vyy*(xsl*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + ysl*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) - zsl*cos(fi)*cos(om)) - vyz*(xsl*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + ysl*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - zsl*sin(om)*cos(fi));
j.coeffRef(1,4) = -vyx*(-xsl*sin(fi)*cos(ka) + ysl*sin(fi)*sin(ka) + zsl*cos(fi)) - vyy*(xsl*sin(om)*cos(fi)*cos(ka) - ysl*sin(ka)*sin(om)*cos(fi) + zsl*sin(fi)*sin(om)) - vyz*(-xsl*cos(fi)*cos(ka)*cos(om) + ysl*sin(ka)*cos(fi)*cos(om) - zsl*sin(fi)*cos(om));
j.coeffRef(1,5) = -vyx*(-xsl*sin(ka)*cos(fi) - ysl*cos(fi)*cos(ka)) - vyy*(xsl*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + ysl*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om))) - vyz*(xsl*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + ysl*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)));
}