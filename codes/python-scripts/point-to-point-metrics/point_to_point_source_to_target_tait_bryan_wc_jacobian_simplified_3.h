#ifndef _point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_3_h_
#define _point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_3_h_
inline void point_to_point_source_to_target_tait_bryan_wc_simplified_3(double &delta_x, double &delta_y, double &delta_z, const double &tx, const double &ty, const double &tz, const double &om, const double &fi, const double &ka, const double &x_s, const double &y_s, const double &z_s, const double &x_t, const double &y_t, const double &z_t)
{
delta_x = -tx - x_s + x_t;
delta_y = -ty - y_s + y_t;
delta_z = -tz - z_s + z_t;
}
inline void point_to_point_source_to_target_tait_bryan_wc_jacobian_simplified_3(Eigen::Matrix<double, 3, 6, Eigen::RowMajor> &j, const double &tx, const double &ty, const double &tz, const double &om, const double &fi, const double &ka, const double &x_s, const double &y_s, const double &z_s)
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
inline void point_to_point_source_to_target_tait_bryan_wc_AtPA_simplified_3(Eigen::Matrix<double, 6, 6, Eigen::RowMajor> &AtPA, const double &tx, const double &ty, const double &tz, const double &om, const double &fi, const double &ka, const double &x_s, const double &y_s, const double &z_s)
{
double x0 = -y_s;
double x1 = -z_s;
double x2 = -x_s;
double x3 = pow(y_s, 2);
double x4 = pow(z_s, 2);
double x5 = -x_s*y_s;
double x6 = -x_s*z_s;
double x7 = pow(x_s, 2);
double x8 = -y_s*z_s;
AtPA.coeffRef(0,0) = 1;
AtPA.coeffRef(0,1) = 0;
AtPA.coeffRef(0,2) = 0;
AtPA.coeffRef(0,3) = 0;
AtPA.coeffRef(0,4) = z_s;
AtPA.coeffRef(0,5) = x0;
AtPA.coeffRef(1,0) = 0;
AtPA.coeffRef(1,1) = 1;
AtPA.coeffRef(1,2) = 0;
AtPA.coeffRef(1,3) = x1;
AtPA.coeffRef(1,4) = 0;
AtPA.coeffRef(1,5) = x_s;
AtPA.coeffRef(2,0) = 0;
AtPA.coeffRef(2,1) = 0;
AtPA.coeffRef(2,2) = 1;
AtPA.coeffRef(2,3) = y_s;
AtPA.coeffRef(2,4) = x2;
AtPA.coeffRef(2,5) = 0;
AtPA.coeffRef(3,0) = 0;
AtPA.coeffRef(3,1) = x1;
AtPA.coeffRef(3,2) = y_s;
AtPA.coeffRef(3,3) = x3 + x4;
AtPA.coeffRef(3,4) = x5;
AtPA.coeffRef(3,5) = x6;
AtPA.coeffRef(4,0) = z_s;
AtPA.coeffRef(4,1) = 0;
AtPA.coeffRef(4,2) = x2;
AtPA.coeffRef(4,3) = x5;
AtPA.coeffRef(4,4) = x4 + x7;
AtPA.coeffRef(4,5) = x8;
AtPA.coeffRef(5,0) = x0;
AtPA.coeffRef(5,1) = x_s;
AtPA.coeffRef(5,2) = 0;
AtPA.coeffRef(5,3) = x6;
AtPA.coeffRef(5,4) = x8;
AtPA.coeffRef(5,5) = x3 + x7;
}
inline void point_to_point_source_to_target_tait_bryan_wc_AtPB_simplified_3(Eigen::Matrix<double, 6, 1> &AtPB, const double &tx, const double &ty, const double &tz, const double &om, const double &fi, const double &ka, const double &x_s, const double &y_s, const double &z_s, const double &x_t, const double &y_t, const double &z_t)
{
double x0 = tx + x_s - x_t;
double x1 = ty + y_s - y_t;
double x2 = tz + z_s - z_t;
AtPB.coeffRef(0) = x0;
AtPB.coeffRef(1) = x1;
AtPB.coeffRef(2) = x2;
AtPB.coeffRef(3) = -x1*z_s + x2*y_s;
AtPB.coeffRef(4) = x0*z_s - x2*x_s;
AtPB.coeffRef(5) = -x0*y_s + x1*x_s;
}
#endif
