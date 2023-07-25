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
AtPA.coeffRef(0,0) = p11;
AtPA.coeffRef(0,1) = p12;
AtPA.coeffRef(0,2) = p13;
AtPA.coeffRef(0,3) = -p12*z_s + p13*y_s;
AtPA.coeffRef(0,4) = p11*z_s - p13*x_s;
AtPA.coeffRef(0,5) = -p11*y_s + p12*x_s;
AtPA.coeffRef(1,0) = p21;
AtPA.coeffRef(1,1) = p22;
AtPA.coeffRef(1,2) = p23;
AtPA.coeffRef(1,3) = -p22*z_s + p23*y_s;
AtPA.coeffRef(1,4) = p21*z_s - p23*x_s;
AtPA.coeffRef(1,5) = -p21*y_s + p22*x_s;
AtPA.coeffRef(2,0) = p31;
AtPA.coeffRef(2,1) = p32;
AtPA.coeffRef(2,2) = p33;
AtPA.coeffRef(2,3) = -p32*z_s + p33*y_s;
AtPA.coeffRef(2,4) = p31*z_s - p33*x_s;
AtPA.coeffRef(2,5) = -p31*y_s + p32*x_s;
AtPA.coeffRef(3,0) = -p21*z_s + p31*y_s;
AtPA.coeffRef(3,1) = -p22*z_s + p32*y_s;
AtPA.coeffRef(3,2) = -p23*z_s + p33*y_s;
AtPA.coeffRef(3,3) = -y_s*(p23*z_s - p33*y_s) + z_s*(p22*z_s - p32*y_s);
AtPA.coeffRef(3,4) = x_s*(p23*z_s - p33*y_s) - z_s*(p21*z_s - p31*y_s);
AtPA.coeffRef(3,5) = -x_s*(p22*z_s - p32*y_s) + y_s*(p21*z_s - p31*y_s);
AtPA.coeffRef(4,0) = p11*z_s - p31*x_s;
AtPA.coeffRef(4,1) = p12*z_s - p32*x_s;
AtPA.coeffRef(4,2) = p13*z_s - p33*x_s;
AtPA.coeffRef(4,3) = y_s*(p13*z_s - p33*x_s) - z_s*(p12*z_s - p32*x_s);
AtPA.coeffRef(4,4) = -x_s*(p13*z_s - p33*x_s) + z_s*(p11*z_s - p31*x_s);
AtPA.coeffRef(4,5) = x_s*(p12*z_s - p32*x_s) - y_s*(p11*z_s - p31*x_s);
AtPA.coeffRef(5,0) = -p11*y_s + p21*x_s;
AtPA.coeffRef(5,1) = -p12*y_s + p22*x_s;
AtPA.coeffRef(5,2) = -p13*y_s + p23*x_s;
AtPA.coeffRef(5,3) = -y_s*(p13*y_s - p23*x_s) + z_s*(p12*y_s - p22*x_s);
AtPA.coeffRef(5,4) = x_s*(p13*y_s - p23*x_s) - z_s*(p11*y_s - p21*x_s);
AtPA.coeffRef(5,5) = -x_s*(p12*y_s - p22*x_s) + y_s*(p11*y_s - p21*x_s);
}
inline void point_to_point_source_to_target_tait_bryan_wc_AtPB_simplified_2(Eigen::Matrix<double, 6, 1> &AtPB, double tx, double ty, double tz, double om, double fi, double ka, double x_s, double y_s, double z_s, double p11, double p12, double p13, double p21, double p22, double p23, double p31, double p32, double p33, double x_t, double y_t, double z_t)
{
AtPB.coeffRef(0) = p11*(tx + x_s - x_t) + p12*(ty + y_s - y_t) + p13*(tz + z_s - z_t);
AtPB.coeffRef(1) = p21*(tx + x_s - x_t) + p22*(ty + y_s - y_t) + p23*(tz + z_s - z_t);
AtPB.coeffRef(2) = p31*(tx + x_s - x_t) + p32*(ty + y_s - y_t) + p33*(tz + z_s - z_t);
AtPB.coeffRef(3) = -(p21*z_s - p31*y_s)*(tx + x_s - x_t) - (p22*z_s - p32*y_s)*(ty + y_s - y_t) - (p23*z_s - p33*y_s)*(tz + z_s - z_t);
AtPB.coeffRef(4) = (p11*z_s - p31*x_s)*(tx + x_s - x_t) + (p12*z_s - p32*x_s)*(ty + y_s - y_t) + (p13*z_s - p33*x_s)*(tz + z_s - z_t);
AtPB.coeffRef(5) = -(p11*y_s - p21*x_s)*(tx + x_s - x_t) - (p12*y_s - p22*x_s)*(ty + y_s - y_t) - (p13*y_s - p23*x_s)*(tz + z_s - z_t);
}
#endif
