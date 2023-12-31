#include <Eigen/Eigen>
#include "structures.h"
#include "transformations.h"
#include <uncertainty_calculator.h>

#include <iostream>

int main(int argc, char *argv[]){
    Eigen::Matrix<double, 6, 6> j_tb;
    TaitBryanPose pose_tb;
    pose_tb.px = 10.01;
    pose_tb.py = 1000.02;
    pose_tb.pz = 110.03;
    pose_tb.om = 3.0 / 180.0 * M_PI;
    pose_tb.fi = 3.0 / 180.0 * M_PI;
    pose_tb.ka = 90.0 / 180.0 * M_PI;

    Eigen::Matrix<double, 6, 6> cov_tb;
    for(int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++)
        {
            cov_tb(i, j) = 0.0;
        }
    }
    cov_tb(0, 0) = 0.01;
    cov_tb(1, 1) = 0.02;
    cov_tb(2, 2) = 0.03;
    cov_tb(3, 3) = pow(3.0 / 180.0 * M_PI, 2);
    cov_tb(4, 4) = pow(4.0 / 180.0 * M_PI, 2);
    cov_tb(5, 5) = pow(5.0 / 180.0 * M_PI, 2);

    std::cout << "cov_tb: " << std::endl << cov_tb << std::endl;

    uncertainty_pose_tait_bryan_to_rodrigues(j_tb, pose_tb.px, pose_tb.py, pose_tb.pz, pose_tb.om, pose_tb.fi, pose_tb.ka);

    std::cout << "j_tb: " << std::endl << j_tb << std::endl;

    Eigen::Matrix<double, 6, 6> cov_rodrigues;
    cov_rodrigues = j_tb * cov_tb * j_tb.transpose();

    std::cout << "cov_rodrigues: " << std::endl
              << cov_rodrigues << std::endl;

    //check
    Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(pose_tb);
    RodriguesPose pose_rodrigues = pose_rodrigues_from_affine_matrix(m);

    Eigen::Matrix<double, 6, 6> j_cov_rodrigues_check;

    uncertainty_pose_rodrigues_to_tait_bryan(j_cov_rodrigues_check, pose_rodrigues.px, pose_rodrigues.py, pose_rodrigues.pz,
                                             pose_rodrigues.sx, pose_rodrigues.sy, pose_rodrigues.sz);

    Eigen::Matrix<double, 6, 6> cov_tb_check;
    cov_tb_check = j_cov_rodrigues_check * cov_rodrigues * j_cov_rodrigues_check.transpose();

    std::cout << "cov_tb_check: " << std::endl
              << cov_tb_check << std::endl;

    std::cout << "check diff (cov_tb_check - cov_tb)" << std::endl << cov_tb_check - cov_tb << std::endl;

    //////////////////////////
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            cov_rodrigues(i, j) = 0.0;
        }
    }
    cov_rodrigues(0, 0) = 0.0001;
    cov_rodrigues(1, 1) = 0.0001;
    cov_rodrigues(2, 2) = 0.0001;
    cov_rodrigues(3, 3) = 0.0001;
    cov_rodrigues(4, 4) = 0.0001;
    cov_rodrigues(5, 5) = 0.0001;

    cov_tb_check;
    cov_tb_check = j_cov_rodrigues_check * cov_rodrigues * j_cov_rodrigues_check.transpose();

    std::cout << "cov_tb_check diag(0.01): " << std::endl
              << cov_tb_check << std::endl;

    return 0;
}