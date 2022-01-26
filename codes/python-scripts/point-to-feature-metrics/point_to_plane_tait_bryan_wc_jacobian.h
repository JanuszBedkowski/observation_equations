inline void point_to_plane_tait_bryan_wc(Eigen::Matrix<double, 1, 1> &delta, double tx, double ty, double tz, double om, double fi, double ka, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vzx, double vzy, double vzz)
{delta.coeffRef(0,0) = -vzx*(tx + xsl*cos(fi)*cos(ka) - xtg - ysl*sin(ka)*cos(fi) + zsl*sin(fi)) - vzy*(ty + xsl*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + ysl*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - ytg - zsl*sin(om)*cos(fi)) - vzz*(tz + xsl*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) + ysl*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + zsl*cos(fi)*cos(om) - ztg);
}
inline void point_to_plane_tait_bryan_wc_jacobian(Eigen::Matrix<double, 1, 6> &j, double tx, double ty, double tz, double om, double fi, double ka, double xsl, double ysl, double zsl, double xtg, double ytg, double ztg, double vzx, double vzy, double vzz)
{j.coeffRef(0,0) = -vzx;
j.coeffRef(0,1) = -vzy;
j.coeffRef(0,2) = -vzz;
j.coeffRef(0,3) = -vzy*(xsl*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + ysl*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) - zsl*cos(fi)*cos(om)) - vzz*(xsl*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + ysl*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - zsl*sin(om)*cos(fi));
j.coeffRef(0,4) = -vzx*(-xsl*sin(fi)*cos(ka) + ysl*sin(fi)*sin(ka) + zsl*cos(fi)) - vzy*(xsl*sin(om)*cos(fi)*cos(ka) - ysl*sin(ka)*sin(om)*cos(fi) + zsl*sin(fi)*sin(om)) - vzz*(-xsl*cos(fi)*cos(ka)*cos(om) + ysl*sin(ka)*cos(fi)*cos(om) - zsl*sin(fi)*cos(om));
j.coeffRef(0,5) = -vzx*(-xsl*sin(ka)*cos(fi) - ysl*cos(fi)*cos(ka)) - vzy*(xsl*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + ysl*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om))) - vzz*(xsl*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + ysl*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)));
}