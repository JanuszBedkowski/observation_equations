inline void point_to_point_tait_bryan_cw(double &delta_x, double &delta_y, double &delta_z, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2)
{delta_x = px_1*cos(fi_1)*cos(ka_1) - px_2*cos(fi_2)*cos(ka_2) - py_1*(-sin(fi_1)*sin(om_1)*cos(ka_1) - sin(ka_1)*cos(om_1)) + py_2*(-sin(fi_2)*sin(om_2)*cos(ka_2) - sin(ka_2)*cos(om_2)) - pz_1*(sin(fi_1)*cos(ka_1)*cos(om_1) - sin(ka_1)*sin(om_1)) + pz_2*(sin(fi_2)*cos(ka_2)*cos(om_2) - sin(ka_2)*sin(om_2)) - x_1*cos(fi_1)*cos(ka_1) + x_2*cos(fi_2)*cos(ka_2) - y_1*(sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1)) + y_2*(sin(fi_2)*sin(om_2)*cos(ka_2) + sin(ka_2)*cos(om_2)) - z_1*(-sin(fi_1)*cos(ka_1)*cos(om_1) + sin(ka_1)*sin(om_1)) + z_2*(-sin(fi_2)*cos(ka_2)*cos(om_2) + sin(ka_2)*sin(om_2));
delta_y = -px_1*sin(ka_1)*cos(fi_1) + px_2*sin(ka_2)*cos(fi_2) - py_1*(sin(fi_1)*sin(ka_1)*sin(om_1) - cos(ka_1)*cos(om_1)) + py_2*(sin(fi_2)*sin(ka_2)*sin(om_2) - cos(ka_2)*cos(om_2)) - pz_1*(-sin(fi_1)*sin(ka_1)*cos(om_1) - sin(om_1)*cos(ka_1)) + pz_2*(-sin(fi_2)*sin(ka_2)*cos(om_2) - sin(om_2)*cos(ka_2)) + x_1*sin(ka_1)*cos(fi_1) - x_2*sin(ka_2)*cos(fi_2) - y_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1)) + y_2*(-sin(fi_2)*sin(ka_2)*sin(om_2) + cos(ka_2)*cos(om_2)) - z_1*(sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1)) + z_2*(sin(fi_2)*sin(ka_2)*cos(om_2) + sin(om_2)*cos(ka_2));
delta_z = px_1*sin(fi_1) - px_2*sin(fi_2) - py_1*sin(om_1)*cos(fi_1) + py_2*sin(om_2)*cos(fi_2) + pz_1*cos(fi_1)*cos(om_1) - pz_2*cos(fi_2)*cos(om_2) - x_1*sin(fi_1) + x_2*sin(fi_2) + y_1*sin(om_1)*cos(fi_1) - y_2*sin(om_2)*cos(fi_2) - z_1*cos(fi_1)*cos(om_1) + z_2*cos(fi_2)*cos(om_2);
}
inline void point_to_point_tait_bryan_cw_jacobian(Eigen::Matrix<double, 3, 12, Eigen::RowMajor> &j, double px_1, double py_1, double pz_1, double om_1, double fi_1, double ka_1, double px_2, double py_2, double pz_2, double om_2, double fi_2, double ka_2, double x_1, double y_1, double z_1, double x_2, double y_2, double z_2)
{j.coeffRef(0,0) = cos(fi_1)*cos(ka_1);
j.coeffRef(0,1) = sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1);
j.coeffRef(0,2) = -sin(fi_1)*cos(ka_1)*cos(om_1) + sin(ka_1)*sin(om_1);
j.coeffRef(0,3) = -py_1*(-sin(fi_1)*cos(ka_1)*cos(om_1) + sin(ka_1)*sin(om_1)) - pz_1*(-sin(fi_1)*sin(om_1)*cos(ka_1) - sin(ka_1)*cos(om_1)) - y_1*(sin(fi_1)*cos(ka_1)*cos(om_1) - sin(ka_1)*sin(om_1)) - z_1*(sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1));
j.coeffRef(0,4) = -px_1*sin(fi_1)*cos(ka_1) + py_1*sin(om_1)*cos(fi_1)*cos(ka_1) - pz_1*cos(fi_1)*cos(ka_1)*cos(om_1) + x_1*sin(fi_1)*cos(ka_1) - y_1*sin(om_1)*cos(fi_1)*cos(ka_1) + z_1*cos(fi_1)*cos(ka_1)*cos(om_1);
j.coeffRef(0,5) = -px_1*sin(ka_1)*cos(fi_1) - py_1*(sin(fi_1)*sin(ka_1)*sin(om_1) - cos(ka_1)*cos(om_1)) - pz_1*(-sin(fi_1)*sin(ka_1)*cos(om_1) - sin(om_1)*cos(ka_1)) + x_1*sin(ka_1)*cos(fi_1) - y_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1)) - z_1*(sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1));
j.coeffRef(0,6) = -cos(fi_2)*cos(ka_2);
j.coeffRef(0,7) = -sin(fi_2)*sin(om_2)*cos(ka_2) - sin(ka_2)*cos(om_2);
j.coeffRef(0,8) = sin(fi_2)*cos(ka_2)*cos(om_2) - sin(ka_2)*sin(om_2);
j.coeffRef(0,9) = py_2*(-sin(fi_2)*cos(ka_2)*cos(om_2) + sin(ka_2)*sin(om_2)) + pz_2*(-sin(fi_2)*sin(om_2)*cos(ka_2) - sin(ka_2)*cos(om_2)) + y_2*(sin(fi_2)*cos(ka_2)*cos(om_2) - sin(ka_2)*sin(om_2)) + z_2*(sin(fi_2)*sin(om_2)*cos(ka_2) + sin(ka_2)*cos(om_2));
j.coeffRef(0,10) = px_2*sin(fi_2)*cos(ka_2) - py_2*sin(om_2)*cos(fi_2)*cos(ka_2) + pz_2*cos(fi_2)*cos(ka_2)*cos(om_2) - x_2*sin(fi_2)*cos(ka_2) + y_2*sin(om_2)*cos(fi_2)*cos(ka_2) - z_2*cos(fi_2)*cos(ka_2)*cos(om_2);
j.coeffRef(0,11) = px_2*sin(ka_2)*cos(fi_2) + py_2*(sin(fi_2)*sin(ka_2)*sin(om_2) - cos(ka_2)*cos(om_2)) + pz_2*(-sin(fi_2)*sin(ka_2)*cos(om_2) - sin(om_2)*cos(ka_2)) - x_2*sin(ka_2)*cos(fi_2) + y_2*(-sin(fi_2)*sin(ka_2)*sin(om_2) + cos(ka_2)*cos(om_2)) + z_2*(sin(fi_2)*sin(ka_2)*cos(om_2) + sin(om_2)*cos(ka_2));
j.coeffRef(1,0) = -sin(ka_1)*cos(fi_1);
j.coeffRef(1,1) = -sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1);
j.coeffRef(1,2) = sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1);
j.coeffRef(1,3) = -py_1*(sin(fi_1)*sin(ka_1)*cos(om_1) + sin(om_1)*cos(ka_1)) - pz_1*(sin(fi_1)*sin(ka_1)*sin(om_1) - cos(ka_1)*cos(om_1)) - y_1*(-sin(fi_1)*sin(ka_1)*cos(om_1) - sin(om_1)*cos(ka_1)) - z_1*(-sin(fi_1)*sin(ka_1)*sin(om_1) + cos(ka_1)*cos(om_1));
j.coeffRef(1,4) = px_1*sin(fi_1)*sin(ka_1) - py_1*sin(ka_1)*sin(om_1)*cos(fi_1) + pz_1*sin(ka_1)*cos(fi_1)*cos(om_1) - x_1*sin(fi_1)*sin(ka_1) + y_1*sin(ka_1)*sin(om_1)*cos(fi_1) - z_1*sin(ka_1)*cos(fi_1)*cos(om_1);
j.coeffRef(1,5) = -px_1*cos(fi_1)*cos(ka_1) - py_1*(sin(fi_1)*sin(om_1)*cos(ka_1) + sin(ka_1)*cos(om_1)) - pz_1*(-sin(fi_1)*cos(ka_1)*cos(om_1) + sin(ka_1)*sin(om_1)) + x_1*cos(fi_1)*cos(ka_1) - y_1*(-sin(fi_1)*sin(om_1)*cos(ka_1) - sin(ka_1)*cos(om_1)) - z_1*(sin(fi_1)*cos(ka_1)*cos(om_1) - sin(ka_1)*sin(om_1));
j.coeffRef(1,6) = sin(ka_2)*cos(fi_2);
j.coeffRef(1,7) = sin(fi_2)*sin(ka_2)*sin(om_2) - cos(ka_2)*cos(om_2);
j.coeffRef(1,8) = -sin(fi_2)*sin(ka_2)*cos(om_2) - sin(om_2)*cos(ka_2);
j.coeffRef(1,9) = py_2*(sin(fi_2)*sin(ka_2)*cos(om_2) + sin(om_2)*cos(ka_2)) + pz_2*(sin(fi_2)*sin(ka_2)*sin(om_2) - cos(ka_2)*cos(om_2)) + y_2*(-sin(fi_2)*sin(ka_2)*cos(om_2) - sin(om_2)*cos(ka_2)) + z_2*(-sin(fi_2)*sin(ka_2)*sin(om_2) + cos(ka_2)*cos(om_2));
j.coeffRef(1,10) = -px_2*sin(fi_2)*sin(ka_2) + py_2*sin(ka_2)*sin(om_2)*cos(fi_2) - pz_2*sin(ka_2)*cos(fi_2)*cos(om_2) + x_2*sin(fi_2)*sin(ka_2) - y_2*sin(ka_2)*sin(om_2)*cos(fi_2) + z_2*sin(ka_2)*cos(fi_2)*cos(om_2);
j.coeffRef(1,11) = px_2*cos(fi_2)*cos(ka_2) + py_2*(sin(fi_2)*sin(om_2)*cos(ka_2) + sin(ka_2)*cos(om_2)) + pz_2*(-sin(fi_2)*cos(ka_2)*cos(om_2) + sin(ka_2)*sin(om_2)) - x_2*cos(fi_2)*cos(ka_2) + y_2*(-sin(fi_2)*sin(om_2)*cos(ka_2) - sin(ka_2)*cos(om_2)) + z_2*(sin(fi_2)*cos(ka_2)*cos(om_2) - sin(ka_2)*sin(om_2));
j.coeffRef(2,0) = sin(fi_1);
j.coeffRef(2,1) = -sin(om_1)*cos(fi_1);
j.coeffRef(2,2) = cos(fi_1)*cos(om_1);
j.coeffRef(2,3) = -py_1*cos(fi_1)*cos(om_1) - pz_1*sin(om_1)*cos(fi_1) + y_1*cos(fi_1)*cos(om_1) + z_1*sin(om_1)*cos(fi_1);
j.coeffRef(2,4) = px_1*cos(fi_1) + py_1*sin(fi_1)*sin(om_1) - pz_1*sin(fi_1)*cos(om_1) - x_1*cos(fi_1) - y_1*sin(fi_1)*sin(om_1) + z_1*sin(fi_1)*cos(om_1);
j.coeffRef(2,5) = 0;
j.coeffRef(2,6) = -sin(fi_2);
j.coeffRef(2,7) = sin(om_2)*cos(fi_2);
j.coeffRef(2,8) = -cos(fi_2)*cos(om_2);
j.coeffRef(2,9) = py_2*cos(fi_2)*cos(om_2) + pz_2*sin(om_2)*cos(fi_2) - y_2*cos(fi_2)*cos(om_2) - z_2*sin(om_2)*cos(fi_2);
j.coeffRef(2,10) = -px_2*cos(fi_2) - py_2*sin(fi_2)*sin(om_2) + pz_2*sin(fi_2)*cos(om_2) + x_2*cos(fi_2) + y_2*sin(fi_2)*sin(om_2) - z_2*sin(fi_2)*cos(om_2);
j.coeffRef(2,11) = 0;
}