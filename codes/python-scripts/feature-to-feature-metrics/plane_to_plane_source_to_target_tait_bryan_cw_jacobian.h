#ifndef _plane_to_plane_source_to_target_tait_bryan_cw_jacobian_h_
#define _plane_to_plane_source_to_target_tait_bryan_cw_jacobian_h_
inline void plane_to_plane_source_to_target_tait_bryan_cw(Eigen::Matrix<double, 4, 1> &delta, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double a_1, double b_1, double c_1, double d_1, double a_2, double b_2, double c_2, double d_2)
{delta.coeffRef(0,0) = -a_1*cos(fi_1)*cos(ka_1) + a_2 - b_1*(sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1)) - c_1*(-sin(fi_1)*cos(ka_1)*cos(om_1) + sin(ka_1)*sin(om_1));
delta.coeffRef(1,0) = a_1*sin(ka_1)*cos(fi_1) - b_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1)) + b_2 - c_1*(sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1));
delta.coeffRef(2,0) = -a_1*sin(fi_1) + b_1*sin(om_1)*cos(fi_1) - c_1*cos(fi_1)*cos(om_1) + c_2;
delta.coeffRef(3,0) = -a_1*tx_1 - b_1*ty_1 - c_1*tz_1 - d_1 + d_2;
}
inline void plane_to_plane_source_to_target_tait_bryan_cw_jacobian(Eigen::Matrix<double, 4, 6, Eigen::RowMajor> &j, double tx_1, double ty_1, double tz_1, double om_1, double fi_1, double ka_1, double a_1, double b_1, double c_1, double d_1, double a_2, double b_2, double c_2, double d_2)
{j.coeffRef(0,0) = 0;
j.coeffRef(0,1) = 0;
j.coeffRef(0,2) = 0;
j.coeffRef(0,3) = -b_1*(sin(fi_1)*cos(ka_1)*cos(om_1) - sin(ka_1)*sin(om_1)) - c_1*(sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1));
j.coeffRef(0,4) = a_1*sin(fi_1)*cos(ka_1) - b_1*sin(om_1)*cos(fi_1)*cos(ka_1) + c_1*cos(fi_1)*cos(ka_1)*cos(om_1);
j.coeffRef(0,5) = a_1*sin(ka_1)*cos(fi_1) - b_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1)) - c_1*(sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1));
j.coeffRef(1,0) = 0;
j.coeffRef(1,1) = 0;
j.coeffRef(1,2) = 0;
j.coeffRef(1,3) = -b_1*(-sin(fi_1)*sin(ka_1)*cos(om_1) - sin(om_1)*cos(ka_1)) - c_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1));
j.coeffRef(1,4) = -a_1*sin(fi_1)*sin(ka_1) + b_1*sin(ka_1)*sin(om_1)*cos(fi_1) - c_1*sin(ka_1)*cos(fi_1)*cos(om_1);
j.coeffRef(1,5) = a_1*cos(fi_1)*cos(ka_1) - b_1*(-sin(fi_1)*sin(om_1)*cos(ka_1) - sin(ka_1)*cos(om_1)) - c_1*(sin(fi_1)*cos(ka_1)*cos(om_1) - sin(ka_1)*sin(om_1));
j.coeffRef(2,0) = 0;
j.coeffRef(2,1) = 0;
j.coeffRef(2,2) = 0;
j.coeffRef(2,3) = b_1*cos(fi_1)*cos(om_1) + c_1*sin(om_1)*cos(fi_1);
j.coeffRef(2,4) = -a_1*cos(fi_1) - b_1*sin(fi_1)*sin(om_1) + c_1*sin(fi_1)*cos(om_1);
j.coeffRef(2,5) = 0;
j.coeffRef(3,0) = -a_1;
j.coeffRef(3,1) = -b_1;
j.coeffRef(3,2) = -c_1;
j.coeffRef(3,3) = 0;
j.coeffRef(3,4) = 0;
j.coeffRef(3,5) = 0;
}
#endif