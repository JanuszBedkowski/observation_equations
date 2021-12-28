inline void point_to_projection_onto_line_tait_bryan_wc(Eigen::Matrix<double, 3, 1> &delta, double px, double py, double pz, double om, double fi, double ka, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double x_trg_ln, double y_trg_ln, double z_trg_ln)
{delta.coeffRef(0,0) = -px - x_src_l*cos(fi)*cos(ka) + x_trg_g + x_trg_ln*(x_trg_ln*(px + x_src_l*cos(fi)*cos(ka) - x_trg_g - y_src_l*sin(ka)*cos(fi) + z_src_l*sin(fi)) + y_trg_ln*(py + x_src_l*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + y_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - y_trg_g - z_src_l*sin(om)*cos(fi)) + z_trg_ln*(pz + x_src_l*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) + y_src_l*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + z_src_l*cos(fi)*cos(om) - z_trg_g))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) + y_src_l*sin(ka)*cos(fi) - z_src_l*sin(fi);
delta.coeffRef(1,0) = -py - x_src_l*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) - y_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + y_trg_g + y_trg_ln*(x_trg_ln*(px + x_src_l*cos(fi)*cos(ka) - x_trg_g - y_src_l*sin(ka)*cos(fi) + z_src_l*sin(fi)) + y_trg_ln*(py + x_src_l*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + y_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - y_trg_g - z_src_l*sin(om)*cos(fi)) + z_trg_ln*(pz + x_src_l*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) + y_src_l*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + z_src_l*cos(fi)*cos(om) - z_trg_g))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) + z_src_l*sin(om)*cos(fi);
delta.coeffRef(2,0) = -pz - x_src_l*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) - y_src_l*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) - z_src_l*cos(fi)*cos(om) + z_trg_g + z_trg_ln*(x_trg_ln*(px + x_src_l*cos(fi)*cos(ka) - x_trg_g - y_src_l*sin(ka)*cos(fi) + z_src_l*sin(fi)) + y_trg_ln*(py + x_src_l*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + y_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - y_trg_g - z_src_l*sin(om)*cos(fi)) + z_trg_ln*(pz + x_src_l*(-sin(fi)*cos(ka)*cos(om) + sin(ka)*sin(om)) + y_src_l*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + z_src_l*cos(fi)*cos(om) - z_trg_g))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
}
inline void point_to_projection_onto_line_tait_bryan_wc_jacobian(Eigen::Matrix<double, 3, 6> &j, double px, double py, double pz, double om, double fi, double ka, double x_src_l, double y_src_l, double z_src_l, double x_trg_g, double y_trg_g, double z_trg_g, double x_trg_ln, double y_trg_ln, double z_trg_ln)
{j.coeffRef(0,0) = pow(x_trg_ln, 2)/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) - 1;
j.coeffRef(0,1) = x_trg_ln*y_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(0,2) = x_trg_ln*z_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(0,3) = x_trg_ln*(y_trg_ln*(x_src_l*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + y_src_l*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) - z_src_l*cos(fi)*cos(om)) + z_trg_ln*(x_src_l*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + y_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - z_src_l*sin(om)*cos(fi)))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(0,4) = x_src_l*sin(fi)*cos(ka) + x_trg_ln*(x_trg_ln*(-x_src_l*sin(fi)*cos(ka) + y_src_l*sin(fi)*sin(ka) + z_src_l*cos(fi)) + y_trg_ln*(x_src_l*sin(om)*cos(fi)*cos(ka) - y_src_l*sin(ka)*sin(om)*cos(fi) + z_src_l*sin(fi)*sin(om)) + z_trg_ln*(-x_src_l*cos(fi)*cos(ka)*cos(om) + y_src_l*sin(ka)*cos(fi)*cos(om) - z_src_l*sin(fi)*cos(om)))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) - y_src_l*sin(fi)*sin(ka) - z_src_l*cos(fi);
j.coeffRef(0,5) = x_src_l*sin(ka)*cos(fi) + x_trg_ln*(x_trg_ln*(-x_src_l*sin(ka)*cos(fi) - y_src_l*cos(fi)*cos(ka)) + y_trg_ln*(x_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + y_src_l*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om))) + z_trg_ln*(x_src_l*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + y_src_l*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om))))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) + y_src_l*cos(fi)*cos(ka);
j.coeffRef(1,0) = x_trg_ln*y_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(1,1) = pow(y_trg_ln, 2)/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) - 1;
j.coeffRef(1,2) = y_trg_ln*z_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(1,3) = -x_src_l*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) - y_src_l*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) + y_trg_ln*(y_trg_ln*(x_src_l*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + y_src_l*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) - z_src_l*cos(fi)*cos(om)) + z_trg_ln*(x_src_l*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + y_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - z_src_l*sin(om)*cos(fi)))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) + z_src_l*cos(fi)*cos(om);
j.coeffRef(1,4) = -x_src_l*sin(om)*cos(fi)*cos(ka) + y_src_l*sin(ka)*sin(om)*cos(fi) + y_trg_ln*(x_trg_ln*(-x_src_l*sin(fi)*cos(ka) + y_src_l*sin(fi)*sin(ka) + z_src_l*cos(fi)) + y_trg_ln*(x_src_l*sin(om)*cos(fi)*cos(ka) - y_src_l*sin(ka)*sin(om)*cos(fi) + z_src_l*sin(fi)*sin(om)) + z_trg_ln*(-x_src_l*cos(fi)*cos(ka)*cos(om) + y_src_l*sin(ka)*cos(fi)*cos(om) - z_src_l*sin(fi)*cos(om)))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) - z_src_l*sin(fi)*sin(om);
j.coeffRef(1,5) = -x_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - y_src_l*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om)) + y_trg_ln*(x_trg_ln*(-x_src_l*sin(ka)*cos(fi) - y_src_l*cos(fi)*cos(ka)) + y_trg_ln*(x_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + y_src_l*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om))) + z_trg_ln*(x_src_l*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + y_src_l*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om))))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,0) = x_trg_ln*z_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,1) = y_trg_ln*z_trg_ln/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,2) = pow(z_trg_ln, 2)/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2)) - 1;
j.coeffRef(2,3) = -x_src_l*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) - y_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + z_src_l*sin(om)*cos(fi) + z_trg_ln*(y_trg_ln*(x_src_l*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + y_src_l*(-sin(fi)*sin(ka)*cos(om) - sin(om)*cos(ka)) - z_src_l*cos(fi)*cos(om)) + z_trg_ln*(x_src_l*(sin(fi)*sin(om)*cos(ka) + sin(ka)*cos(om)) + y_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) - z_src_l*sin(om)*cos(fi)))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,4) = x_src_l*cos(fi)*cos(ka)*cos(om) - y_src_l*sin(ka)*cos(fi)*cos(om) + z_src_l*sin(fi)*cos(om) + z_trg_ln*(x_trg_ln*(-x_src_l*sin(fi)*cos(ka) + y_src_l*sin(fi)*sin(ka) + z_src_l*cos(fi)) + y_trg_ln*(x_src_l*sin(om)*cos(fi)*cos(ka) - y_src_l*sin(ka)*sin(om)*cos(fi) + z_src_l*sin(fi)*sin(om)) + z_trg_ln*(-x_src_l*cos(fi)*cos(ka)*cos(om) + y_src_l*sin(ka)*cos(fi)*cos(om) - z_src_l*sin(fi)*cos(om)))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
j.coeffRef(2,5) = -x_src_l*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) - y_src_l*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om)) + z_trg_ln*(x_trg_ln*(-x_src_l*sin(ka)*cos(fi) - y_src_l*cos(fi)*cos(ka)) + y_trg_ln*(x_src_l*(-sin(fi)*sin(ka)*sin(om) + cos(ka)*cos(om)) + y_src_l*(-sin(fi)*sin(om)*cos(ka) - sin(ka)*cos(om))) + z_trg_ln*(x_src_l*(sin(fi)*sin(ka)*cos(om) + sin(om)*cos(ka)) + y_src_l*(sin(fi)*cos(ka)*cos(om) - sin(ka)*sin(om))))/(pow(x_trg_ln, 2) + pow(y_trg_ln, 2) + pow(z_trg_ln, 2));
}