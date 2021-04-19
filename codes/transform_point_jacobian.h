inline void transform_point3D_Tait_Bryan_jacobian(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double x, double y, double z, double px, double py, double pz, double om, double fi, double ka)
{j.coeffRef(0,0) = 1;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = 0;
j.coeffRef(0,4) = -x*sin(fi)*cos(ka) + y*sin(fi)*sin(ka) + z*cos(fi);
j.coeffRef(0,5) = -(x*sin(ka) + y*cos(ka))*cos(fi);
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = 1;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = x*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) - y*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) - z*cos(fi)*cos(om);
j.coeffRef(1,4) = (x*cos(fi)*cos(ka) - y*sin(ka)*cos(fi) + z*sin(fi))*sin(om);
j.coeffRef(1,5) = -x*(sin(fi)*sin(ka)*sin(om) - cos(ka)*cos(om)) - y*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om));
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = 1;
j.coeffRef(2,3) = x*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) - y*(sin(fi)*sin(ka)*sin(om) - cos(ka)*cos(om)) - z*sin(om)*cos(fi);
j.coeffRef(2,4) = (-x*cos(fi)*cos(ka) + y*sin(ka)*cos(fi) - z*sin(fi))*cos(om);
j.coeffRef(2,5) = x*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + y*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om));
}