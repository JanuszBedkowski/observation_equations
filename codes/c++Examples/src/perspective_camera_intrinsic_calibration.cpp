#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"
#include "perspective_camera_with_intrinsics_tait_bryan_cw_jacobian.h"
#include "perspective_camera_with_intrinsics_rodrigues_cw_jacobian.h"
#include "perspective_camera_with_intrinsics_quaternion_cw_jacobian.h"
#include "quaternion_constraint_jacobian.h"
#include "perspective_camera_tait_bryan_wc_jacobian.h"

double k1 = 0.1686;
double k2 = -0.515827;
double k3 = 0.446983;
double k4 = 0.0;
double k5 = 0.0;
double k6 = 0.0;
double p1 = -0.00138148;
double p2 = 0.000127175;

double k1tmp;
double k2tmp;
double k3tmp;
double k4tmp;
double k5tmp;
double k6tmp;
double p1tmp;
double p2tmp;

struct KeyPoint{
	double u;
	double v;
	std::pair<int,int> index_to_tie_point;
};

struct Camera{
	Eigen::Affine3d pose;
	std::vector<std::vector<KeyPoint>> key_points;
};

std::vector<std::vector<Eigen::Vector3d>> tie_points;
std::vector<Camera> cameras;
PerspectiveCameraParams cam_params;

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = -62, rotate_y = 21;
float translate_z = -12.4502;
float translate_x = -0.899999, translate_y = -1.49993;
bool show_only_one = true;

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

int main(int argc, char *argv[]){
	k1tmp = k1;
	k2tmp = k2;
	k3tmp = k3;
	p1tmp = p1;
	p2tmp = p2;

	if (false == initGL(&argc, argv)) {
		return 4;
	}

	cam_params.fx = 1727;
	cam_params.fy = 1727;
	cam_params.cx = 985;
	cam_params.cy = 522;

	for(size_t i = 0 ; i < 8; i++){
		std::vector<Eigen::Vector3d> tps;
		for(size_t j = 0 ; j < 8; j++){
		Eigen::Vector3d p;
		p.x() = i;
		p.y() = j;
		p.z() = 0;
		tps.push_back(p);
		}
		tie_points.push_back(tps);
	}

	for(float i = 3 ; i < 5; i+=0.5){
		Camera c;
		c.pose = Eigen::Affine3d::Identity();
		c.pose(0,3) = 4;
		c.pose(1,3) = i;
		c.pose(2,3) = -10;
		cameras.push_back(c);
	}

	for(size_t i = 0 ; i < cameras.size(); i++){
		for(size_t j = 0; j < tie_points.size(); j++){
			std::vector<KeyPoint> kps;
			for(size_t k = 0; k < tie_points[j].size(); k++){
				KeyPoint kp;
				kp.index_to_tie_point = std::pair<int,int>(j,k);

				TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose.inverse());
				projection_perspective_camera_with_intrinsics_tait_bryan_cw(kp.u, kp.v, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
						tie_points[j][k].x(), tie_points[j][k].y(), tie_points[j][k].z(),  k1,  k2,  k3, p1, p2);
				kps.push_back(kp);
			}
			cameras[i].key_points.push_back(kps);
		}
	}

	k1 = 0.0;
	k2 = 0.0;
	p1 = 0.0;
	p2 = 0.0;
	k3 = 0.0;

	printHelp();
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutMainLoop();

	return 0;
}

bool initGL(int *argc, char **argv) {
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow("perspective_camera_intrinsic_calibration");
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMotionFunc(motion);

	// default initialization
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glEnable(GL_DEPTH_TEST);

	// viewport
	glViewport(0, 0, window_width, window_height);

	// projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (GLfloat) window_width / (GLfloat) window_height, 0.01,
			10000.0);
	glutReshapeFunc(reshape);

	return true;
}

Eigen::Vector3d get_intersection(float u, float v, Eigen::Affine3d camera_pose_wc){
	cv::Mat mat(1, 2, CV_32F);
	mat.at<float>(0, 0) = u;
	mat.at<float>(0, 1) = v;

	cv::Mat cv_cam_matrix = (cv::Mat_<float>(3, 3) << cam_params.fx, 0, cam_params.cx, 0, cam_params.fy, cam_params.cy, 0, 0, 1);
	cv::Mat cv_dist_params = (cv::Mat_<float>(5, 1) << k1, k2, p1, p2, k3);

	mat = mat.reshape(2);
	cv::undistortPoints(mat, mat, cv_cam_matrix, cv_dist_params, cv::Mat(), cv_cam_matrix,
						cv::TermCriteria(cv::TermCriteria::EPS | cv::TermCriteria::MAX_ITER, 20, 1e-6));
	mat = mat.reshape(1);

	double undistorted_u = mat.at<float>(0, 0);
	double undistorted_v = mat.at<float>(0, 1);

	Eigen::Matrix3d K;
	K(0,0) = cam_params.fx;
	K(0,1) = 0;
	K(0,2) = cam_params.cx;
	K(1,0) = 0;
	K(1,1) = cam_params.fy;
	K(1,2) = cam_params.cy;
	K(2,0) = 0;
	K(2,1) = 0;
	K(2,2) = 1;

	Eigen::Vector3d p(undistorted_u, undistorted_v, 1);
	Eigen::Vector3d r =  K.inverse() * p;
	Eigen::Matrix3d R = camera_pose_wc.rotation();
	Eigen::Vector3d T = camera_pose_wc.translation();
	Eigen::Vector3d rt = R*r + T;

	Eigen::Vector3d n(0,0,-1);
	Eigen::Vector3d a(camera_pose_wc(0,3), camera_pose_wc(1,3), camera_pose_wc(2,3));
	Eigen::Vector3d b = rt;

	Eigen::Vector3d ba = b-a;
	Eigen::Vector3d intersection = a + (((0 - n.dot(a))/n.dot(ba)) * ba);
	return intersection;
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(translate_x, translate_y, translate_z);
	glRotatef(rotate_x, 1.0, 0.0, 0.0);
	glRotatef(rotate_y, 0.0, 0.0, 1.0);

	glLineWidth(2);
	glColor3f(1,0,0);
	glBegin(GL_LINES);
	for(size_t i = 0 ; i < tie_points.size(); i++){
		for(size_t j = 0 ; j < tie_points[i].size(); j++){
			if(j+1 < tie_points[i].size()){
				glVertex3f(tie_points[i][j].x(), tie_points[i][j].y(), tie_points[i][j].z());
				glVertex3f(tie_points[i][j+1].x(), tie_points[i][j+1].y(), tie_points[i][j+1].z());
			}

			if(i+1 < tie_points.size()){
				glVertex3f(tie_points[i][j].x(), tie_points[i][j].y(), tie_points[i][j].z());
				glVertex3f(tie_points[i+1][j].x(), tie_points[i+1][j].y(), tie_points[i+1][j].z());
			}
		}
	}
	glEnd();

	for(size_t i = 0 ; i < cameras.size(); i++){
		Eigen::Affine3d m = cameras[i].pose;

		glBegin(GL_LINES);
			glColor3f(1.0f, 0.0f, 0.0f);
			glVertex3f(m(0,3), m(1,3), m(2,3));
			glVertex3f(m(0,3) + m(0,0), m(1,3) + m(1,0), m(2,3) + m(2,0));

			glColor3f(0.0f, 1.0f, 0.0f);
			glVertex3f(m(0,3), m(1,3), m(2,3));
			glVertex3f(m(0,3) + m(0,1), m(1,3) + m(1,1), m(2,3) + m(2,1));

			glColor3f(0.0f, 0.0f, 1.0f);
			glVertex3f(m(0,3), m(1,3), m(2,3));
			glVertex3f(m(0,3) + m(0,2), m(1,3) + m(1,2), m(2,3) + m(2,2));
		glEnd();
	}

	int number_cameras = 1;
	if(!show_only_one){
		number_cameras = cameras.size();
	}
	glColor3f(0,1,0);
	glBegin(GL_LINES);
	for(size_t i = 0; i < number_cameras; i++){
		for(size_t j = 0; j < cameras[i].key_points.size(); j++){
			for(size_t k = 0; k < cameras[i].key_points[j].size(); k++){
				if(k+1 < tie_points[j].size()){
					Eigen::Vector3d intersection = get_intersection(cameras[i].key_points[j][k].u, cameras[i].key_points[j][k].v, cameras[i].pose);
					glVertex3f(intersection.x(), intersection.y(), intersection.z());

					intersection = get_intersection(cameras[i].key_points[j][k+1].u, cameras[i].key_points[j][k+1].v, cameras[i].pose);
					glVertex3f(intersection.x(), intersection.y(), intersection.z());
				}

				if(j+1 < tie_points.size()){
					Eigen::Vector3d intersection = get_intersection(cameras[i].key_points[j][k].u, cameras[i].key_points[j][k].v, cameras[i].pose);
					glVertex3f(intersection.x(), intersection.y(), intersection.z());

					intersection = get_intersection(cameras[i].key_points[j+1][k].u, cameras[i].key_points[j+1][k].v, cameras[i].pose);
					glVertex3f(intersection.x(), intersection.y(), intersection.z());
				}
			}
		}
	}
	glEnd();
	glPointSize(1);

	glutSwapBuffers();
}

void keyboard(unsigned char key, int /*x*/, int /*y*/) {
	switch (key) {
		case (27): {
			glutDestroyWindow(glutGetWindow());
			return;
		}
		case 'c':{
			for(size_t i = 0; i < cameras.size(); i++){
				TaitBryanPose pose;
				pose.px = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1;
				pose.py = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1;
				pose.pz = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1;
				pose.om = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.1;
				pose.fi = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.1;
				pose.ka = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.1;

				Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(pose);
				cameras[i].pose = cameras[i].pose * m;
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			int number_intrinsic_paramters = 5;


			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras[i].key_points.size(); j++){
					for(size_t k = 0; k < cameras[i].key_points[j].size(); k++){
						TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose.inverse());

						Eigen::Matrix<double, 2, 1> delta;
						observation_equation_perspective_camera_with_intrinsics_tait_bryan_cw(delta, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy,
								pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].x(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].y(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].z(),
								std::floor(cameras[i].key_points[j][k].u), std::floor(cameras[i].key_points[j][k].v), k1, k2, k3, p1, p2);

						Eigen::Matrix<double, 2, 14> jacobian;
						observation_equation_perspective_camera_with_intrinsics_tait_bryan_cw_jacobian(jacobian, cam_params.fx, cam_params.fy,
								pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].x(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].y(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].z(),
								k1, k2, k3, p1, p2);

						int ir = tripletListB.size();
						int ic_camera = i * 6;

						for(size_t l = 0; l < number_intrinsic_paramters; l++){
							tripletListA.emplace_back(ir     , l, -jacobian(0,l));
						}
						for(size_t l = 0; l < number_intrinsic_paramters; l++){
							tripletListA.emplace_back(ir + 1  , l, -jacobian(1,l));
						}

						tripletListA.emplace_back(ir     , ic_camera     + number_intrinsic_paramters, -jacobian(0,0 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 1 + number_intrinsic_paramters, -jacobian(0,1 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 2 + number_intrinsic_paramters, -jacobian(0,2 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 3 + number_intrinsic_paramters, -jacobian(0,3 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 4 + number_intrinsic_paramters, -jacobian(0,4 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 5 + number_intrinsic_paramters, -jacobian(0,5 + number_intrinsic_paramters));

						tripletListA.emplace_back(ir + 1 , ic_camera     + number_intrinsic_paramters, -jacobian(1,0 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 1 + number_intrinsic_paramters, -jacobian(1,1 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 2 + number_intrinsic_paramters, -jacobian(1,2 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 3 + number_intrinsic_paramters, -jacobian(1,3 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 4 + number_intrinsic_paramters, -jacobian(1,4 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 5 + number_intrinsic_paramters, -jacobian(1,5 + number_intrinsic_paramters));

						tripletListP.emplace_back(ir    , ir    , cauchy(delta(0,0), 1));
						tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta(1,0), 1));

						tripletListB.emplace_back(ir    , 0,  delta(0,0));
						tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
					}
				}
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 6 + number_intrinsic_paramters);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 6 + number_intrinsic_paramters, cameras.size() * 6 + number_intrinsic_paramters);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 6 + number_intrinsic_paramters, 1);

			{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP) * matA;
			AtPB = (AtP) * matB;
			}

			tripletListA.clear();
			tripletListP.clear();
			tripletListB.clear();


			std::cout << "AtPA.size: " << AtPA.size() << std::endl;
			std::cout << "AtPB.size: " << AtPB.size() << std::endl;

			std::cout << "start solving AtPA=AtPB" << std::endl;
			Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(AtPA);

			std::cout << "x = solver.solve(AtPB)" << std::endl;
			Eigen::SparseMatrix<double> x = solver.solve(AtPB);

			std::vector<double> h_x;

			for (int k=0; k<x.outerSize(); ++k){
				for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it){
					h_x.push_back(it.value());
				}
			}

			if(h_x.size() == cameras.size() * 6 + number_intrinsic_paramters){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				k1 += h_x[counter++];
				k2 += h_x[counter++];
				k3 += h_x[counter++];
				p1 += h_x[counter++];
				p2 += h_x[counter++];

				for(size_t i = 0; i < cameras.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose.inverse());
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_tait_bryan(pose).inverse();
				}

				std::cout << "desired intrinsic parameters:" << std::endl;
				std::cout << "k1: " << k1tmp << " k2: " << k2tmp << " k3: " << k3tmp << " k4: " << k4tmp << " k5: " << k5tmp << " k6: " << k6tmp << " p1: " << p1tmp << " p2: " << p2tmp << std::endl;
				std::cout << "computer intrinsic parameters:" << std::endl;
				std::cout << "k1: " << k1 << " k2: " << k2 << " k3: " << k3 << " k4: " << k4 << " k5: " << k5 << " k6: " << k6 << " p1: " << p1 << " p2: " << p2 << std::endl;

			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'r':{
			for(size_t i = 0; i < cameras.size(); i++){
				TaitBryanPose posetb = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
				posetb.px += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.py += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.pz += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.om += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.fi += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.ka += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				cameras[i].pose = affine_matrix_from_pose_tait_bryan(posetb);
			}

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			int number_intrinsic_paramters = 5;


			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras[i].key_points.size(); j++){
					for(size_t k = 0; k < cameras[i].key_points[j].size(); k++){

						RodriguesPose pose = pose_rodrigues_from_affine_matrix(cameras[i].pose.inverse());

						Eigen::Matrix<double, 2, 1> delta;
						observation_equation_perspective_camera_with_intrinsics_rodrigues_cw(delta, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy,
								pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].x(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].y(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].z(),
								std::floor(cameras[i].key_points[j][k].u), std::floor(cameras[i].key_points[j][k].v), k1, k2, k3, p1, p2);

						Eigen::Matrix<double, 2, 14> jacobian;
						observation_equation_perspective_camera_with_intrinsics_rodrigues_cw_jacobian(jacobian, cam_params.fx, cam_params.fy,
								pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].x(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].y(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].z(),
								k1, k2, k3, p1, p2);

						int ir = tripletListB.size();
						int ic_camera = i * 6;

						tripletListA.emplace_back(ir     , 0, -jacobian(0,0));
						tripletListA.emplace_back(ir     , 1, -jacobian(0,1));
						tripletListA.emplace_back(ir     , 2, -jacobian(0,2));
						tripletListA.emplace_back(ir     , 3, -jacobian(0,3));
						tripletListA.emplace_back(ir     , 4, -jacobian(0,4));

						tripletListA.emplace_back(ir + 1 , 0, -jacobian(1,0));
						tripletListA.emplace_back(ir + 1 , 1, -jacobian(1,1));
						tripletListA.emplace_back(ir + 1 , 2, -jacobian(1,2));
						tripletListA.emplace_back(ir + 1 , 3, -jacobian(1,3));
						tripletListA.emplace_back(ir + 1 , 4, -jacobian(1,4));

						tripletListA.emplace_back(ir     , ic_camera     + number_intrinsic_paramters, -jacobian(0,0 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 1 + number_intrinsic_paramters, -jacobian(0,1 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 2 + number_intrinsic_paramters, -jacobian(0,2 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 3 + number_intrinsic_paramters, -jacobian(0,3 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 4 + number_intrinsic_paramters, -jacobian(0,4 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 5 + number_intrinsic_paramters, -jacobian(0,5 + number_intrinsic_paramters));

						tripletListA.emplace_back(ir + 1 , ic_camera     + number_intrinsic_paramters, -jacobian(1,0 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 1 + number_intrinsic_paramters, -jacobian(1,1 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 2 + number_intrinsic_paramters, -jacobian(1,2 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 3 + number_intrinsic_paramters, -jacobian(1,3 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 4 + number_intrinsic_paramters, -jacobian(1,4 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 5 + number_intrinsic_paramters, -jacobian(1,5 + number_intrinsic_paramters));

						tripletListP.emplace_back(ir    , ir    , cauchy(delta(0,0), 1));
						tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta(1,0), 1));

						tripletListB.emplace_back(ir    , 0,  delta(0,0));
						tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
					}
				}
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 6 + number_intrinsic_paramters);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 6 + number_intrinsic_paramters, cameras.size() * 6 + number_intrinsic_paramters);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 6 + number_intrinsic_paramters, 1);

			{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP) * matA;
			AtPB = (AtP) * matB;
			}

			tripletListA.clear();
			tripletListP.clear();
			tripletListB.clear();


			std::cout << "AtPA.size: " << AtPA.size() << std::endl;
			std::cout << "AtPB.size: " << AtPB.size() << std::endl;

			std::cout << "start solving AtPA=AtPB" << std::endl;
			Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(AtPA);

			std::cout << "x = solver.solve(AtPB)" << std::endl;
			Eigen::SparseMatrix<double> x = solver.solve(AtPB);

			std::vector<double> h_x;

			for (int k=0; k<x.outerSize(); ++k){
				for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it){
					h_x.push_back(it.value());
				}
			}

			if(h_x.size() == cameras.size() * 6 + number_intrinsic_paramters){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				k1 += h_x[counter++];
				k2 += h_x[counter++];
				k3 += h_x[counter++];
				p1 += h_x[counter++];
				p2 += h_x[counter++];

				for(size_t i = 0; i < cameras.size(); i++){
					RodriguesPose pose = pose_rodrigues_from_affine_matrix(cameras[i].pose.inverse());
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.sx += h_x[counter++];
					pose.sy += h_x[counter++];
					pose.sz += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_rodrigues(pose).inverse();
				}

				std::cout << "desired intrinsic parameters:" << std::endl;
				std::cout << "k1: " << k1tmp << " k2: " << k2tmp << " k3: " << k3tmp << " p1: " << p1tmp << " p2: " << p2tmp << std::endl;
				std::cout << "computer intrinsic parameters:" << std::endl;
				std::cout << "k1: " << k1 << " k2: " << k2 << " k3: " << k3 << " p1: " << p1 << " p2: " << p2 << std::endl;

			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'q':{
			for(size_t i = 0; i < cameras.size(); i++){
				TaitBryanPose posetb = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
				posetb.px += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;
				posetb.py += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;
				posetb.pz += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;
				posetb.om += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;
				posetb.fi += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;
				posetb.ka += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;
				cameras[i].pose = affine_matrix_from_pose_tait_bryan(posetb);
			}
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			int number_intrinsic_paramters = 5;

			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras[i].key_points.size(); j++){
					for(size_t k = 0; k < cameras[i].key_points[j].size(); k++){
						QuaternionPose pose = pose_quaternion_from_affine_matrix(cameras[i].pose.inverse());

						Eigen::Matrix<double, 2, 1> delta;
						observation_equation_perspective_camera_with_intrinsics_quaternion_cw(delta, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy,
								pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].x(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].y(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].z(),
								std::floor(cameras[i].key_points[j][k].u), std::floor(cameras[i].key_points[j][k].v), k1, k2, k3, p1, p2);

						Eigen::Matrix<double, 2, 15> jacobian;
						observation_equation_perspective_camera_with_intrinsics_quaternion_cw_jacobian(jacobian, cam_params.fx, cam_params.fy,
								pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].x(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].y(),
								tie_points[cameras[i].key_points[j][k].index_to_tie_point.first][cameras[i].key_points[j][k].index_to_tie_point.second].z(),
								k1, k2, k3, p1, p2);

						int ir = tripletListB.size();
						int ic_camera = i * 7;

						tripletListA.emplace_back(ir     , 0, -jacobian(0,0));
						tripletListA.emplace_back(ir     , 1, -jacobian(0,1));
						tripletListA.emplace_back(ir     , 2, -jacobian(0,2));
						tripletListA.emplace_back(ir     , 3, -jacobian(0,3));
						tripletListA.emplace_back(ir     , 4, -jacobian(0,4));

						tripletListA.emplace_back(ir + 1 , 0, -jacobian(1,0));
						tripletListA.emplace_back(ir + 1 , 1, -jacobian(1,1));
						tripletListA.emplace_back(ir + 1 , 2, -jacobian(1,2));
						tripletListA.emplace_back(ir + 1 , 3, -jacobian(1,3));
						tripletListA.emplace_back(ir + 1 , 4, -jacobian(1,4));

						tripletListA.emplace_back(ir     , ic_camera     + number_intrinsic_paramters, -jacobian(0,0 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 1 + number_intrinsic_paramters, -jacobian(0,1 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 2 + number_intrinsic_paramters, -jacobian(0,2 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 3 + number_intrinsic_paramters, -jacobian(0,3 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 4 + number_intrinsic_paramters, -jacobian(0,4 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 5 + number_intrinsic_paramters, -jacobian(0,5 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir     , ic_camera + 6 + number_intrinsic_paramters, -jacobian(0,6 + number_intrinsic_paramters));

						tripletListA.emplace_back(ir + 1 , ic_camera     + number_intrinsic_paramters, -jacobian(1,0 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 1 + number_intrinsic_paramters, -jacobian(1,1 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 2 + number_intrinsic_paramters, -jacobian(1,2 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 3 + number_intrinsic_paramters, -jacobian(1,3 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 4 + number_intrinsic_paramters, -jacobian(1,4 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 5 + number_intrinsic_paramters, -jacobian(1,5 + number_intrinsic_paramters));
						tripletListA.emplace_back(ir + 1 , ic_camera + 6 + number_intrinsic_paramters, -jacobian(1,6 + number_intrinsic_paramters));

						tripletListP.emplace_back(ir    , ir    , cauchy(delta(0,0), 1));
						tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta(1,0), 1));

						tripletListB.emplace_back(ir    , 0,  delta(0,0));
						tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
					}
				}
			}
			std::cout << "joo2" << std::endl;
			for(size_t i = 0 ; i < cameras.size(); i++){
				int ic = i * 7 + number_intrinsic_paramters;
				int ir = tripletListB.size();
				QuaternionPose pose = pose_quaternion_from_affine_matrix(cameras[i].pose);

				double delta;
				quaternion_constraint(delta, pose.q0, pose.q1, pose.q2, pose.q3);

				Eigen::Matrix<double, 1, 4> jacobian;
				quaternion_constraint_jacobian(jacobian, pose.q0, pose.q1, pose.q2, pose.q3);

				tripletListA.emplace_back(ir, ic + 3 , -jacobian(0,0));
				tripletListA.emplace_back(ir, ic + 4 , -jacobian(0,1));
				tripletListA.emplace_back(ir, ic + 5 , -jacobian(0,2));
				tripletListA.emplace_back(ir, ic + 6 , -jacobian(0,3));

				tripletListP.emplace_back(ir, ir, 1000000.0);

				tripletListB.emplace_back(ir, 0, delta);
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 7 + number_intrinsic_paramters);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 7 + number_intrinsic_paramters, cameras.size() * 7 + number_intrinsic_paramters);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 7 + number_intrinsic_paramters, 1);

			{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP) * matA;
			AtPB = (AtP) * matB;
			}

			tripletListA.clear();
			tripletListP.clear();
			tripletListB.clear();


			std::cout << "AtPA.size: " << AtPA.size() << std::endl;
			std::cout << "AtPB.size: " << AtPB.size() << std::endl;

			std::cout << "start solving AtPA=AtPB" << std::endl;
			Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(AtPA);

			std::cout << "x = solver.solve(AtPB)" << std::endl;
			Eigen::SparseMatrix<double> x = solver.solve(AtPB);

			std::vector<double> h_x;

			for (int k=0; k<x.outerSize(); ++k){
				for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it){
					h_x.push_back(it.value());
				}
			}

			if(h_x.size() == cameras.size() * 7 + number_intrinsic_paramters){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				k1 += h_x[counter++];
				k2 += h_x[counter++];
				k3 += h_x[counter++];
				p1 += h_x[counter++];
				p2 += h_x[counter++];

				for(size_t i = 0; i < cameras.size(); i++){
					QuaternionPose pose = pose_quaternion_from_affine_matrix(cameras[i].pose.inverse());
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.q0 += h_x[counter++];
					pose.q1 += h_x[counter++];
					pose.q2 += h_x[counter++];
					pose.q3 += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_quaternion(pose).inverse();
				}

				std::cout << "desired intrinsic parameters:" << std::endl;
				std::cout << "k1: " << k1tmp << " k2: " << k2tmp << " k3: " << k3tmp << " p1: " << p1tmp << " p2: " << p2tmp << std::endl;
				std::cout << "computer intrinsic parameters:" << std::endl;
				std::cout << "k1: " << k1 << " k2: " << k2 << " k3: " << k3 << " p1: " << p1 << " p2: " << p2 << std::endl;

			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case '1':{
			k1 -= 0.01;
			std::cout << "k1: " << k1 << std::endl;
			break;
		}
		case '2':{
			k1 += 0.01;
			std::cout << "k1: " << k1 << std::endl;
			break;
		}
		case '3':{
			k2 -= 0.01;
			std::cout << "k2: " << k2 << std::endl;
			break;
		}
		case '4':{
			k2 += 0.01;
			std::cout << "k2: " << k2 << std::endl;
			break;
		}
		case '5':{
			k3 -= 0.01;
			std::cout << "k3: " << k3 << std::endl;
			break;
		}
		case '6':{
			k3 += 0.01;
			std::cout << "k3: " << k3 << std::endl;
			break;
		}
		case '7':{
			p1 -= 0.01;
			std::cout << "p1: " << p1 << std::endl;
			break;
		}
		case '8':{
			p1 += 0.01;
			std::cout << "p1: " << p1 << std::endl;
			break;
		}
		case '9':{
			p2 -= 0.01;
			std::cout << "p2: " << p2 << std::endl;
			break;
		}
		case '0':{
			p2 += 0.01;
			std::cout << "p2: " << p2 << std::endl;
			break;
		}
		case 'e':{
			show_only_one =! show_only_one;
			break;
		}
	}

	printHelp();
	glutPostRedisplay();
}


void mouse(int button, int state, int x, int y) {
	if (state == GLUT_DOWN) {
		mouse_buttons |= 1 << button;
	} else if (state == GLUT_UP) {
		mouse_buttons = 0;
	}

	mouse_old_x = x;
	mouse_old_y = y;
}

void motion(int x, int y) {
	float dx, dy;
	dx = (float) (x - mouse_old_x);
	dy = (float) (y - mouse_old_y);

	if (mouse_buttons & 1) {
		rotate_x += dy * 0.2f;
		rotate_y += dx * 0.2f;

	} else if (mouse_buttons & 4) {
		translate_z += dy * 0.05f;
	} else if (mouse_buttons & 3) {
		translate_x += dx * 0.05f;
		translate_y -= dy * 0.05f;
	}

	mouse_old_x = x;
	mouse_old_y = y;

	glutPostRedisplay();
}

void reshape(int w, int h) {
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (GLfloat) w / (GLfloat) h, 0.01, 10000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void printHelp() {
	std::cout << "-------help-------" << std::endl;
	std::cout << "t: optimize (Tait-Bryan wc)" << std::endl;
	std::cout << "r: optimize (Rodriguez cw)" << std::endl;
	std::cout << "q: optimize (Quaternion cw)" << std::endl;
	std::cout << "c: add noise to cameras" << std::endl;
	std::cout << "1: k1 -= 0.01" << std::endl;
	std::cout << "2: k1 += 0.01" << std::endl;
	std::cout << "3: k2 -= 0.01" << std::endl;
	std::cout << "4: k2 += 0.01" << std::endl;
	std::cout << "5: k3 -= 0.01" << std::endl;
	std::cout << "6: k3 += 0.01" << std::endl;
	std::cout << "7: p1 -= 0.01" << std::endl;
	std::cout << "8: p1 += 0.01" << std::endl;
	std::cout << "9: p2 -= 0.01" << std::endl;
	std::cout << "0: p2 += 0.01" << std::endl;
	std::cout << "e: show_only_one=!show_only_one" << std::endl;
}




