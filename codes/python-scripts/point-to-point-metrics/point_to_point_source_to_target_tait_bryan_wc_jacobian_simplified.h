#ifndef _point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_h_
#define _point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_h_
inline void point_to_point_source_to_target_tait_bryan_wc_simplified(double &delta_x, double &delta_y, double &delta_z, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)
{
double sin_om = sin(om);
double cos_om = cos(om);
double sin_fi = sin(fi);
double cos_fi = cos(fi);
double sin_ka = sin(ka);
double cos_ka = cos(ka);
double x0 = cos_fi*z_s;
double x1 = cos_ka*x_s;
double x2 = sin_ka*y_s;
double x3 = cos_ka*sin_om;
double x4 = cos_om*sin_ka;
double x5 = sin_fi*x4 + x3;
double x6 = sin_ka*sin_om;
double x7 = cos_ka*cos_om;
double x8 = sin_fi*x7 - x6;
double x9 = sin_fi*z_s;
double x10 = cos_fi*x2;
double x11 = cos_fi*x1;
double x12 = sin_fi*x3 + x4;
double x13 = -sin_fi*x6 + x7;
delta_x = -cos_fi*cos_ka*x_s + cos_fi*sin_ka*y_s - sin_fi*z_s - tx + x_t;
delta_y = cos_fi*sin_om*z_s - ty - x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) - y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) + y_t;
delta_z = -cos_fi*cos_om*z_s - tz + x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) - y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + z_t;
}
inline void point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s)
{
double sin_om = sin(om);
double cos_om = cos(om);
double sin_fi = sin(fi);
double cos_fi = cos(fi);
double sin_ka = sin(ka);
double cos_ka = cos(ka);
double x0 = cos_fi*z_s;
double x1 = cos_ka*x_s;
double x2 = sin_ka*y_s;
double x3 = cos_ka*sin_om;
double x4 = cos_om*sin_ka;
double x5 = sin_fi*x4 + x3;
double x6 = sin_ka*sin_om;
double x7 = cos_ka*cos_om;
double x8 = sin_fi*x7 - x6;
double x9 = sin_fi*z_s;
double x10 = cos_fi*x2;
double x11 = cos_fi*x1;
double x12 = sin_fi*x3 + x4;
double x13 = -sin_fi*x6 + x7;
j.coeffRef(0,0) = -1;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = 0;
j.coeffRef(0,4) = -cos_fi*z_s + cos_ka*sin_fi*x_s - sin_fi*sin_ka*y_s;
j.coeffRef(0,5) = cos_fi*(cos_ka*y_s + sin_ka*x_s);
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = -1;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka);
j.coeffRef(1,4) = sin_om*(-cos_fi*cos_ka*x_s + cos_fi*sin_ka*y_s - sin_fi*z_s);
j.coeffRef(1,5) = -x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) + y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka);
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = -1;
j.coeffRef(2,3) = cos_fi*sin_om*z_s - x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) - y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om);
j.coeffRef(2,4) = cos_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s);
j.coeffRef(2,5) = -x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) - y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om);
}
inline void point_to_point_source_to_target_tait_bryan_wc_AtPA_simplified(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &AtPA, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double p11, double p12, double p13, double p21, double p22, double p23, double p31, double p32, double p33)
{
double sin_om = sin(om);
double cos_om = cos(om);
double sin_fi = sin(fi);
double cos_fi = cos(fi);
double sin_ka = sin(ka);
double cos_ka = cos(ka);
double x0 = cos_fi*z_s;
double x1 = cos_ka*x_s;
double x2 = sin_ka*y_s;
double x3 = cos_ka*sin_om;
double x4 = cos_om*sin_ka;
double x5 = sin_fi*x4 + x3;
double x6 = sin_ka*sin_om;
double x7 = cos_ka*cos_om;
double x8 = sin_fi*x7 - x6;
double x9 = sin_fi*z_s;
double x10 = cos_fi*x2;
double x11 = cos_fi*x1;
double x12 = sin_fi*x3 + x4;
double x13 = -sin_fi*x6 + x7;
AtPA.coeffRef(0,0) = p11;
AtPA.coeffRef(0,1) = p12;
AtPA.coeffRef(0,2) = p13;
AtPA.coeffRef(0,3) = -p12*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) + p13*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om));
AtPA.coeffRef(0,4) = -cos_om*p13*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p11*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p12*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s);
AtPA.coeffRef(0,5) = -cos_fi*p11*(cos_ka*y_s + sin_ka*x_s) + p12*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p13*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om));
AtPA.coeffRef(1,0) = p21;
AtPA.coeffRef(1,1) = p22;
AtPA.coeffRef(1,2) = p23;
AtPA.coeffRef(1,3) = -p22*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) + p23*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om));
AtPA.coeffRef(1,4) = -cos_om*p23*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p21*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p22*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s);
AtPA.coeffRef(1,5) = -cos_fi*p21*(cos_ka*y_s + sin_ka*x_s) + p22*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p23*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om));
AtPA.coeffRef(2,0) = p31;
AtPA.coeffRef(2,1) = p32;
AtPA.coeffRef(2,2) = p33;
AtPA.coeffRef(2,3) = -p32*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) + p33*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om));
AtPA.coeffRef(2,4) = -cos_om*p33*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p31*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p32*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s);
AtPA.coeffRef(2,5) = -cos_fi*p31*(cos_ka*y_s + sin_ka*x_s) + p32*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p33*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om));
AtPA.coeffRef(3,0) = -p21*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) + p31*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om));
AtPA.coeffRef(3,1) = -p22*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) + p32*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om));
AtPA.coeffRef(3,2) = -p23*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) + p33*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om));
AtPA.coeffRef(3,3) = (p22*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p32*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om)))*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - (p23*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p33*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om)))*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om));
AtPA.coeffRef(3,4) = cos_om*(p23*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p33*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om)))*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) - sin_om*(p22*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p32*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om)))*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) - (p21*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p31*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om)))*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s);
AtPA.coeffRef(3,5) = cos_fi*(cos_ka*y_s + sin_ka*x_s)*(p21*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p31*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om))) - (p22*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p32*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om)))*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) - (p23*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p33*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om)))*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om));
AtPA.coeffRef(4,0) = -cos_om*p31*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p11*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p21*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s);
AtPA.coeffRef(4,1) = -cos_om*p32*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p12*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p22*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s);
AtPA.coeffRef(4,2) = -cos_om*p33*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p13*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p23*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s);
AtPA.coeffRef(4,3) = -(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka))*(-cos_om*p32*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p12*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p22*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s)) + (-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om))*(-cos_om*p33*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p13*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p23*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s));
AtPA.coeffRef(4,4) = -cos_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s)*(-cos_om*p33*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p13*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p23*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s)) + sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s)*(-cos_om*p32*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p12*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p22*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s)) + (cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s)*(-cos_om*p31*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p11*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p21*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s));
AtPA.coeffRef(4,5) = -cos_fi*(cos_ka*y_s + sin_ka*x_s)*(-cos_om*p31*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p11*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p21*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s)) + (x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka))*(-cos_om*p32*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p12*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p22*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s)) + (x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om))*(-cos_om*p33*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p13*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p23*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s));
AtPA.coeffRef(5,0) = -cos_fi*p11*(cos_ka*y_s + sin_ka*x_s) + p21*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p31*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om));
AtPA.coeffRef(5,1) = -cos_fi*p12*(cos_ka*y_s + sin_ka*x_s) + p22*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p32*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om));
AtPA.coeffRef(5,2) = -cos_fi*p13*(cos_ka*y_s + sin_ka*x_s) + p23*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p33*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om));
AtPA.coeffRef(5,3) = -(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka))*(-cos_fi*p12*(cos_ka*y_s + sin_ka*x_s) + p22*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p32*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om))) + (-cos_fi*p13*(cos_ka*y_s + sin_ka*x_s) + p23*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p33*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om)))*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om));
AtPA.coeffRef(5,4) = -cos_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s)*(-cos_fi*p13*(cos_ka*y_s + sin_ka*x_s) + p23*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p33*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om))) + sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s)*(-cos_fi*p12*(cos_ka*y_s + sin_ka*x_s) + p22*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p32*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om))) + (cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s)*(-cos_fi*p11*(cos_ka*y_s + sin_ka*x_s) + p21*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p31*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om)));
AtPA.coeffRef(5,5) = -cos_fi*(cos_ka*y_s + sin_ka*x_s)*(-cos_fi*p11*(cos_ka*y_s + sin_ka*x_s) + p21*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p31*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om))) + (x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka))*(-cos_fi*p12*(cos_ka*y_s + sin_ka*x_s) + p22*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p32*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om))) + (x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om))*(-cos_fi*p13*(cos_ka*y_s + sin_ka*x_s) + p23*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p33*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om)));
}
inline void point_to_point_source_to_target_tait_bryan_wc_AtPB_simplified(Eigen::Matrix<double, 6, 1> &AtPB, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double p11, double p12, double p13, double p21, double p22, double p23, double p31, double p32, double p33, double x_t, double y_t, double z_t)
{
double sin_om = sin(om);
double cos_om = cos(om);
double sin_fi = sin(fi);
double cos_fi = cos(fi);
double sin_ka = sin(ka);
double cos_ka = cos(ka);
double x0 = cos_fi*z_s;
double x1 = cos_ka*x_s;
double x2 = sin_ka*y_s;
double x3 = cos_ka*sin_om;
double x4 = cos_om*sin_ka;
double x5 = sin_fi*x4 + x3;
double x6 = sin_ka*sin_om;
double x7 = cos_ka*cos_om;
double x8 = sin_fi*x7 - x6;
double x9 = sin_fi*z_s;
double x10 = cos_fi*x2;
double x11 = cos_fi*x1;
double x12 = sin_fi*x3 + x4;
double x13 = -sin_fi*x6 + x7;
AtPB.coeffRef(0) = p11*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s + tx - x_t) + p12*(-cos_fi*sin_om*z_s + ty + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_t) + p13*(cos_fi*cos_om*z_s + tz - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) - z_t);
AtPB.coeffRef(1) = p21*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s + tx - x_t) + p22*(-cos_fi*sin_om*z_s + ty + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_t) + p23*(cos_fi*cos_om*z_s + tz - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) - z_t);
AtPB.coeffRef(2) = p31*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s + tx - x_t) + p32*(-cos_fi*sin_om*z_s + ty + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_t) + p33*(cos_fi*cos_om*z_s + tz - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) - z_t);
AtPB.coeffRef(3) = -(p21*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p31*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om)))*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s + tx - x_t) - (p22*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p32*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om)))*(-cos_fi*sin_om*z_s + ty + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_t) - (p23*(cos_fi*cos_om*z_s - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka)) - p33*(-cos_fi*sin_om*z_s + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om)))*(cos_fi*cos_om*z_s + tz - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) - z_t);
AtPB.coeffRef(4) = (-cos_om*p31*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p11*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p21*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s))*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s + tx - x_t) + (-cos_om*p32*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p12*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p22*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s))*(-cos_fi*sin_om*z_s + ty + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_t) + (-cos_om*p33*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s) + p13*(cos_fi*z_s - cos_ka*sin_fi*x_s + sin_fi*sin_ka*y_s) + p23*sin_om*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s))*(cos_fi*cos_om*z_s + tz - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) - z_t);
AtPB.coeffRef(5) = (-cos_fi*p11*(cos_ka*y_s + sin_ka*x_s) + p21*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p31*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om)))*(cos_fi*cos_ka*x_s - cos_fi*sin_ka*y_s + sin_fi*z_s + tx - x_t) + (-cos_fi*p12*(cos_ka*y_s + sin_ka*x_s) + p22*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p32*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om)))*(-cos_fi*sin_om*z_s + ty + x_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka) + y_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_t) + (-cos_fi*p13*(cos_ka*y_s + sin_ka*x_s) + p23*(x_s*(cos_ka*cos_om - sin_fi*sin_ka*sin_om) - y_s*(cos_ka*sin_fi*sin_om + cos_om*sin_ka)) + p33*(x_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) + y_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om)))*(cos_fi*cos_om*z_s + tz - x_s*(cos_ka*cos_om*sin_fi - sin_ka*sin_om) + y_s*(cos_ka*sin_om + cos_om*sin_fi*sin_ka) - z_t);
}
#endif
