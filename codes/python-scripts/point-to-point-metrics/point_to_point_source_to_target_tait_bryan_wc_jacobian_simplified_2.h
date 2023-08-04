#ifndef _point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_2_h_
#define _point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_2_h_
inline void point_to_point_source_to_target_tait_bryan_wc_simplified_2(double &delta_x, double &delta_y, double &delta_z, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double x_t, double y_t, double z_t)
{
delta_x = -tx - x_s + x_t;
delta_y = -ty - y_s + y_t;
delta_z = -tz - z_s + z_t;
}
inline void point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_2(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s)
{
j.coeffRef(0,0) = -1;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = 0;
j.coeffRef(0,4) = -z_s;
j.coeffRef(0,5) = y_s;
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = -1;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = z_s;
j.coeffRef(1,4) = 0;
j.coeffRef(1,5) = -x_s;
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = -1;
j.coeffRef(2,3) = -y_s;
j.coeffRef(2,4) = x_s;
j.coeffRef(2,5) = 0;
}
inline void point_to_point_source_to_target_tait_bryan_wc_AtPA_simplified_2(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &AtPA, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double p11, double p12, double p13, double p21, double p22, double p23, double p31, double p32, double p33)
{
double x0 = p12*z_s;
double x1 = p13*y_s;
double x2 = p11*z_s;
double x3 = p11*y_s;
double x4 = p22*z_s;
double x5 = p21*z_s;
double x6 = -p23*x_s;
double x7 = -p22*x_s;
double x8 = -p33*y_s;
double x9 = -p33*x_s;
double x10 = p31*y_s;
double x11 = -p32*x_s;
double x12 = -x10 + x5;
double x13 = -p32*y_s + x4;
double x14 = p23*z_s + x8;
double x15 = -p31*x_s + x2;
double x16 = x0 + x11;
double x17 = p13*z_s + x9;
double x18 = -p21*x_s + x3;
double x19 = p12*y_s + x7;
double x20 = x1 + x6;
AtPA.coeffRef(0,0) = p11;
AtPA.coeffRef(0,1) = p12;
AtPA.coeffRef(0,2) = p13;
AtPA.coeffRef(0,3) = -x0 + x1;
AtPA.coeffRef(0,4) = -p13*x_s + x2;
AtPA.coeffRef(0,5) = p12*x_s - x3;
AtPA.coeffRef(1,0) = p21;
AtPA.coeffRef(1,1) = p22;
AtPA.coeffRef(1,2) = p23;
AtPA.coeffRef(1,3) = p23*y_s - x4;
AtPA.coeffRef(1,4) = x5 + x6;
AtPA.coeffRef(1,5) = -p21*y_s - x7;
AtPA.coeffRef(2,0) = p31;
AtPA.coeffRef(2,1) = p32;
AtPA.coeffRef(2,2) = p33;
AtPA.coeffRef(2,3) = -p32*z_s - x8;
AtPA.coeffRef(2,4) = p31*z_s + x9;
AtPA.coeffRef(2,5) = -x10 - x11;
AtPA.coeffRef(3,0) = -x12;
AtPA.coeffRef(3,1) = -x13;
AtPA.coeffRef(3,2) = -x14;
AtPA.coeffRef(3,3) = x13*z_s - x14*y_s;
AtPA.coeffRef(3,4) = -x12*z_s + x14*x_s;
AtPA.coeffRef(3,5) = x12*y_s - x13*x_s;
AtPA.coeffRef(4,0) = x15;
AtPA.coeffRef(4,1) = x16;
AtPA.coeffRef(4,2) = x17;
AtPA.coeffRef(4,3) = -x16*z_s + x17*y_s;
AtPA.coeffRef(4,4) = x15*z_s - x17*x_s;
AtPA.coeffRef(4,5) = -x15*y_s + x16*x_s;
AtPA.coeffRef(5,0) = -x18;
AtPA.coeffRef(5,1) = -x19;
AtPA.coeffRef(5,2) = -x20;
AtPA.coeffRef(5,3) = x19*z_s - x20*y_s;
AtPA.coeffRef(5,4) = -x18*z_s + x20*x_s;
AtPA.coeffRef(5,5) = x18*y_s - x19*x_s;
}
inline void point_to_point_source_to_target_tait_bryan_wc_AtPB_simplified_2(Eigen::Matrix<double, 6, 1> &AtPB, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double p11, double p12, double p13, double p21, double p22, double p23, double p31, double p32, double p33, double x_t, double y_t, double z_t)
{
double x0 = tx + x_s - x_t;
double x1 = ty + y_s - y_t;
double x2 = tz + z_s - z_t;
AtPB.coeffRef(0) = p11*x0 + p12*x1 + p13*x2;
AtPB.coeffRef(1) = p21*x0 + p22*x1 + p23*x2;
AtPB.coeffRef(2) = p31*x0 + p32*x1 + p33*x2;
AtPB.coeffRef(3) = -x0*(p21*z_s - p31*y_s) - x1*(p22*z_s - p32*y_s) - x2*(p23*z_s - p33*y_s);
AtPB.coeffRef(4) = x0*(p11*z_s - p31*x_s) + x1*(p12*z_s - p32*x_s) + x2*(p13*z_s - p33*x_s);
AtPB.coeffRef(5) = -x0*(p11*y_s - p21*x_s) - x1*(p12*y_s - p22*x_s) - x2*(p13*y_s - p23*x_s);
}
#endif
