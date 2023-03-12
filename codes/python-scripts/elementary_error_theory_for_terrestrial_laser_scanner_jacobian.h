#ifndef __ELEMENTARY_ERROR_THEORY_FOR_TERRESTRIAL_LASER_SCANNER_JACOBIAN_H__
#define __ELEMENTARY_ERROR_THEORY_FOR_TERRESTRIAL_LASER_SCANNER_JACOBIAN_H__
inline void elementary_error_theory_for_terrestrial_laser_scanner_jacobian(Eigen::Matrix<double, 3, 3> &j, double r, double alpha, double theta)
{
j(0,0) = sin(alpha)*cos(theta);
j(0,1) = r*cos(alpha)*cos(theta);
j(0,2) = -r*sin(alpha)*sin(theta);
j(1,0) = sin(alpha)*sin(theta);
j(1,1) = r*sin(theta)*cos(alpha);
j(1,2) = r*sin(alpha)*cos(theta);
j(2,0) = cos(theta);
j(2,1) = 0;
j(2,2) = -r*sin(theta);
}
#endif