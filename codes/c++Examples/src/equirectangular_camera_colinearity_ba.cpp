#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>

#include "equirectangular_camera_colinearity_tait_bryan_wc_jacobian.h"
#include "equirectangular_camera_colinearity_rodrigues_wc_jacobian.h"
#include "equirectangular_camera_colinearity_quaternion_wc_jacobian.h"
#include "quaternion_constraint_jacobian.h"
#include "structures.h"
#include "transformations.h"
#include "cauchy.h"

struct KeyPoint{
	double u;
	double v;
	int index_to_tie_point;
};

struct Camera{
	Eigen::Affine3d pose;
	std::vector<KeyPoint> key_points;
};

std::vector<Eigen::Vector3d> tie_points;
std::vector<Camera> cameras;
int cols = 4096;
int rows = 2048;

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = -215, rotate_y = 273;
float translate_z = -68.0;
float translate_x = 4, translate_y = 30.0;

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

int main(int argc, char *argv[]){

	if (false == initGL(&argc, argv)) {
		return 4;
	}

	for(size_t i = 0 ; i < 100; i++){
		Eigen::Vector3d p;
		p.x() = ((rand()%1000000)/1000000.0 - 0.5) * 2.0 * 100.0;
		p.y() = ((rand()%1000000)/1000000.0 - 0.5) * 2.0 * 200.0;
		p.z() = ((rand()%1000000)/1000000.0 - 0.5) * 2.0 * 100.0;
		tie_points.push_back(p);
	}

	for(int i = -50 ; i < 50; i+=20){
		Camera c;
		c.pose = Eigen::Affine3d::Identity();
		c.pose(0,3) = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1.0;
		c.pose(1,3) = i;
		c.pose(2,3) = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1.0;
		cameras.push_back(c);
	}

	for(size_t i = 0 ; i < cameras.size(); i++){
		for(size_t j = 0; j < tie_points.size(); j++){
			KeyPoint kp;
			kp.index_to_tie_point = j;

			TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
			equrectangular_camera_colinearity_tait_bryan_wc(kp.u, kp.v, rows, cols, M_PI, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka, tie_points[j].x(), tie_points[j].y(), tie_points[j].z());

			cameras[i].key_points.push_back(kp);
		}
	}

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
	glutCreateWindow("equirectangular_camera_colinearity_ba");
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

void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(translate_x, translate_y, translate_z);
	glRotatef(rotate_x, 1.0, 0.0, 0.0);
	glRotatef(rotate_y, 0.0, 0.0, 1.0);

	glColor3f(1,0,0);
	glPointSize(5);
	glBegin(GL_POINTS);
	for(size_t i = 0 ; i < tie_points.size(); i++){
		glVertex3f(tie_points[i].x(), tie_points[i].y(), tie_points[i].z());
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

	glColor3f(0.8,0.8,0.8);
	glBegin(GL_LINE_STRIP);
	for(size_t i = 0 ; i < cameras.size(); i++){
		Eigen::Affine3d m = cameras[i].pose;
		glVertex3f(m(0,3), m(1,3), m(2,3));
	}
	glEnd();

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
				pose.px = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.1;
				pose.py = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.1;
				pose.pz = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.1;
				pose.om = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.01;
				pose.fi = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.01;
				pose.ka = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.01;

				Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(pose);
				cameras[i].pose = cameras[i].pose * m;
			}
			break;
		}
		case 'p':{
			for(size_t i = 0; i < tie_points.size(); i++){
				tie_points[i].x() += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1.0;
				tie_points[i].y() += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1.0;
				tie_points[i].z() += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1.0;
			}
			break;
		}
		case 't':{
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

			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras[i].key_points.size(); j++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
					double u;
					double v;
					equrectangular_camera_colinearity_tait_bryan_wc(u, v, rows, cols, M_PI, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka, tie_points[cameras[i].key_points[j].index_to_tie_point].x(), tie_points[cameras[i].key_points[j].index_to_tie_point].y(), tie_points[cameras[i].key_points[j].index_to_tie_point].z());

					if(fabs(std::floor(cameras[i].key_points[j].u) - u) > 100 || std::floor(int(cameras[i].key_points[j].v) - v) > 100){
						continue;
					}

					Eigen::Matrix<double, 2, 1> delta;
					observation_equation_equrectangular_camera_colinearity_tait_bryan_wc(delta, rows, cols, M_PI, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
							tie_points[cameras[i].key_points[j].index_to_tie_point].x(), tie_points[cameras[i].key_points[j].index_to_tie_point].y(), tie_points[cameras[i].key_points[j].index_to_tie_point].z(),
							std::floor(cameras[i].key_points[j].u), std::floor(cameras[i].key_points[j].v));

					Eigen::Matrix<double, 2, 9, Eigen::RowMajor> jacobian;
					observation_equation_equrectangular_camera_colinearity_tait_bryan_wc_jacobian(jacobian, rows, cols, M_PI, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
							tie_points[cameras[i].key_points[j].index_to_tie_point].x(), tie_points[cameras[i].key_points[j].index_to_tie_point].y(), tie_points[cameras[i].key_points[j].index_to_tie_point].z(),
							std::floor(cameras[i].key_points[j].u), std::floor(cameras[i].key_points[j].v));

					int ir = tripletListB.size();
					int ic_camera = i * 6;
					int ic_tie_point = cameras.size() * 6 + cameras[i].key_points[j].index_to_tie_point * 3;

					tripletListA.emplace_back(ir     , ic_camera    , -jacobian(0,0));
					tripletListA.emplace_back(ir     , ic_camera + 1, -jacobian(0,1));
					tripletListA.emplace_back(ir     , ic_camera + 2, -jacobian(0,2));
					tripletListA.emplace_back(ir     , ic_camera + 3, -jacobian(0,3));
					tripletListA.emplace_back(ir     , ic_camera + 4, -jacobian(0,4));
					tripletListA.emplace_back(ir     , ic_camera + 5, -jacobian(0,5));

					tripletListA.emplace_back(ir     , ic_tie_point,     -jacobian(0,6));
					tripletListA.emplace_back(ir     , ic_tie_point + 1, -jacobian(0,7));
					tripletListA.emplace_back(ir     , ic_tie_point + 2, -jacobian(0,8));

					tripletListA.emplace_back(ir + 1 , ic_camera    , -jacobian(1,0));
					tripletListA.emplace_back(ir + 1 , ic_camera + 1, -jacobian(1,1));
					tripletListA.emplace_back(ir + 1 , ic_camera + 2, -jacobian(1,2));
					tripletListA.emplace_back(ir + 1 , ic_camera + 3, -jacobian(1,3));
					tripletListA.emplace_back(ir + 1 , ic_camera + 4, -jacobian(1,4));
					tripletListA.emplace_back(ir + 1 , ic_camera + 5, -jacobian(1,5));

					tripletListA.emplace_back(ir + 1 , ic_tie_point,     -jacobian(1,6));
					tripletListA.emplace_back(ir + 1 , ic_tie_point + 1, -jacobian(1,7));
					tripletListA.emplace_back(ir + 1 , ic_tie_point + 2, -jacobian(1,8));

					tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
					tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));

					tripletListB.emplace_back(ir    , 0,  delta(0,0));
					tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
				}
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);
			tripletListA.emplace_back(ir + 3 , 3, 1);
			tripletListA.emplace_back(ir + 4 , 4, 1);
			tripletListA.emplace_back(ir + 5 , 5, 1);

			tripletListP.emplace_back(ir     , ir,     1000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 1000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 1000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 1000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 1000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 1000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 6 + tie_points.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 6 + tie_points.size() * 3, cameras.size() * 6 + tie_points.size() * 3);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 6 + tie_points.size() * 3, 1);

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

			if(h_x.size() == cameras.size() * 6 + tie_points.size() * 3){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < cameras.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_tait_bryan(pose);
				}

				for(size_t i = 0; i < tie_points.size(); i++){
					tie_points[i].x() += h_x[counter++];
					tie_points[i].y() += h_x[counter++];
					tie_points[i].z() += h_x[counter++];
				}
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

			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras[i].key_points.size(); j++){
					TaitBryanPose pose_check = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
					double u;
					double v;
					equrectangular_camera_colinearity_tait_bryan_wc(u, v, rows, cols, M_PI, pose_check.px, pose_check.py, pose_check.pz, pose_check.om, pose_check.fi, pose_check.ka, tie_points[cameras[i].key_points[j].index_to_tie_point].x(), tie_points[cameras[i].key_points[j].index_to_tie_point].y(), tie_points[cameras[i].key_points[j].index_to_tie_point].z());

					if(fabs(std::floor(cameras[i].key_points[j].u) - u) > 100 || std::floor(int(cameras[i].key_points[j].v) - v) > 100){
						continue;
					}

					RodriguesPose pose = pose_rodrigues_from_affine_matrix(cameras[i].pose);

					Eigen::Matrix<double, 2, 1> delta;
					observation_equation_equrectangular_camera_colinearity_rodrigues_wc(delta, rows, cols, M_PI, pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
							tie_points[cameras[i].key_points[j].index_to_tie_point].x(), tie_points[cameras[i].key_points[j].index_to_tie_point].y(), tie_points[cameras[i].key_points[j].index_to_tie_point].z(),
							std::floor(cameras[i].key_points[j].u), std::floor(cameras[i].key_points[j].v));


					Eigen::Matrix<double, 2, 9, Eigen::RowMajor> jacobian;
					observation_equation_equrectangular_camera_colinearity_rodrigues_wc_jacobian(jacobian, rows, cols, M_PI, pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
							tie_points[cameras[i].key_points[j].index_to_tie_point].x(), tie_points[cameras[i].key_points[j].index_to_tie_point].y(), tie_points[cameras[i].key_points[j].index_to_tie_point].z(),
							std::floor(cameras[i].key_points[j].u), std::floor(cameras[i].key_points[j].v));

					int ir = tripletListB.size();
					int ic_camera = i * 6;
					int ic_tie_point = cameras.size() * 6 + cameras[i].key_points[j].index_to_tie_point * 3;

					tripletListA.emplace_back(ir     , ic_camera    , -jacobian(0,0));
					tripletListA.emplace_back(ir     , ic_camera + 1, -jacobian(0,1));
					tripletListA.emplace_back(ir     , ic_camera + 2, -jacobian(0,2));
					tripletListA.emplace_back(ir     , ic_camera + 3, -jacobian(0,3));
					tripletListA.emplace_back(ir     , ic_camera + 4, -jacobian(0,4));
					tripletListA.emplace_back(ir     , ic_camera + 5, -jacobian(0,5));

					tripletListA.emplace_back(ir     , ic_tie_point,     -jacobian(0,6));
					tripletListA.emplace_back(ir     , ic_tie_point + 1, -jacobian(0,7));
					tripletListA.emplace_back(ir     , ic_tie_point + 2, -jacobian(0,8));

					tripletListA.emplace_back(ir + 1 , ic_camera    , -jacobian(1,0));
					tripletListA.emplace_back(ir + 1 , ic_camera + 1, -jacobian(1,1));
					tripletListA.emplace_back(ir + 1 , ic_camera + 2, -jacobian(1,2));
					tripletListA.emplace_back(ir + 1 , ic_camera + 3, -jacobian(1,3));
					tripletListA.emplace_back(ir + 1 , ic_camera + 4, -jacobian(1,4));
					tripletListA.emplace_back(ir + 1 , ic_camera + 5, -jacobian(1,5));

					tripletListA.emplace_back(ir + 1 , ic_tie_point,     -jacobian(1,6));
					tripletListA.emplace_back(ir + 1 , ic_tie_point + 1, -jacobian(1,7));
					tripletListA.emplace_back(ir + 1 , ic_tie_point + 2, -jacobian(1,8));

					tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
					tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));

					tripletListB.emplace_back(ir    , 0,  delta(0,0));
					tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
				}
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);
			tripletListA.emplace_back(ir + 3 , 3, 1);
			tripletListA.emplace_back(ir + 4 , 4, 1);
			tripletListA.emplace_back(ir + 5 , 5, 1);

			tripletListP.emplace_back(ir     , ir,     1000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 1000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 1000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 1000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 1000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 1000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 6 + tie_points.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 6 + tie_points.size() * 3, cameras.size() * 6 + tie_points.size() * 3);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 6 + tie_points.size() * 3, 1);

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

			if(h_x.size() == cameras.size() * 6 + tie_points.size() * 3){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < cameras.size(); i++){
					RodriguesPose pose = pose_rodrigues_from_affine_matrix(cameras[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.sx += h_x[counter++];
					pose.sy += h_x[counter++];
					pose.sz += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_rodrigues(pose);
				}

				for(size_t i = 0; i < tie_points.size(); i++){
					tie_points[i].x() += h_x[counter++];
					tie_points[i].y() += h_x[counter++];
					tie_points[i].z() += h_x[counter++];
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'q':{
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


			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras[i].key_points.size(); j++){
					TaitBryanPose pose_check = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
					double u;
					double v;
					equrectangular_camera_colinearity_tait_bryan_wc(u, v, rows, cols, M_PI, pose_check.px, pose_check.py, pose_check.pz, pose_check.om, pose_check.fi, pose_check.ka, tie_points[cameras[i].key_points[j].index_to_tie_point].x(), tie_points[cameras[i].key_points[j].index_to_tie_point].y(), tie_points[cameras[i].key_points[j].index_to_tie_point].z());

					if(fabs(std::floor(cameras[i].key_points[j].u) - u) > 100 || std::floor(int(cameras[i].key_points[j].v) - v) > 100){
						continue;
					}

					QuaternionPose pose = pose_quaternion_from_affine_matrix(cameras[i].pose);

					Eigen::Matrix<double, 2, 1> delta;
					observation_equation_equrectangular_camera_colinearity_quaternion_wc(delta, rows, cols, M_PI, pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
							tie_points[cameras[i].key_points[j].index_to_tie_point].x(), tie_points[cameras[i].key_points[j].index_to_tie_point].y(), tie_points[cameras[i].key_points[j].index_to_tie_point].z(),
							std::floor(cameras[i].key_points[j].u), std::floor(cameras[i].key_points[j].v));

					Eigen::Matrix<double, 2, 10, Eigen::RowMajor> jacobian;
					observation_equation_equrectangular_camera_colinearity_quaternion_wc_jacobian(jacobian, rows, cols, M_PI, pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
							tie_points[cameras[i].key_points[j].index_to_tie_point].x(), tie_points[cameras[i].key_points[j].index_to_tie_point].y(), tie_points[cameras[i].key_points[j].index_to_tie_point].z(),
							std::floor(cameras[i].key_points[j].u), std::floor(cameras[i].key_points[j].v));

					int ir = tripletListB.size();
					int ic_camera = i * 7;
					int ic_tie_point = cameras.size() * 7 + cameras[i].key_points[j].index_to_tie_point * 3;

					tripletListA.emplace_back(ir     , ic_camera    , -jacobian(0,0));
					tripletListA.emplace_back(ir     , ic_camera + 1, -jacobian(0,1));
					tripletListA.emplace_back(ir     , ic_camera + 2, -jacobian(0,2));
					tripletListA.emplace_back(ir     , ic_camera + 3, -jacobian(0,3));
					tripletListA.emplace_back(ir     , ic_camera + 4, -jacobian(0,4));
					tripletListA.emplace_back(ir     , ic_camera + 5, -jacobian(0,5));
					tripletListA.emplace_back(ir     , ic_camera + 6, -jacobian(0,6));

					tripletListA.emplace_back(ir     , ic_tie_point,     -jacobian(0,7));
					tripletListA.emplace_back(ir     , ic_tie_point + 1, -jacobian(0,8));
					tripletListA.emplace_back(ir     , ic_tie_point + 2, -jacobian(0,9));

					tripletListA.emplace_back(ir + 1 , ic_camera    , -jacobian(1,0));
					tripletListA.emplace_back(ir + 1 , ic_camera + 1, -jacobian(1,1));
					tripletListA.emplace_back(ir + 1 , ic_camera + 2, -jacobian(1,2));
					tripletListA.emplace_back(ir + 1 , ic_camera + 3, -jacobian(1,3));
					tripletListA.emplace_back(ir + 1 , ic_camera + 4, -jacobian(1,4));
					tripletListA.emplace_back(ir + 1 , ic_camera + 5, -jacobian(1,5));
					tripletListA.emplace_back(ir + 1 , ic_camera + 6, -jacobian(1,6));

					tripletListA.emplace_back(ir + 1 , ic_tie_point,     -jacobian(1,7));
					tripletListA.emplace_back(ir + 1 , ic_tie_point + 1, -jacobian(1,8));
					tripletListA.emplace_back(ir + 1 , ic_tie_point + 2, -jacobian(1,9));

					tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
					tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));

					tripletListB.emplace_back(ir    , 0,  delta(0,0));
					tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
				}
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);
			tripletListA.emplace_back(ir + 3 , 3, 1);
			tripletListA.emplace_back(ir + 4 , 4, 1);
			tripletListA.emplace_back(ir + 5 , 5, 1);
			tripletListA.emplace_back(ir + 6 , 6, 1);

			tripletListP.emplace_back(ir     , ir,     1000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 1000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 1000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 1000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 1000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 1000000);
			tripletListP.emplace_back(ir + 6 , ir + 6, 1000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);
			tripletListB.emplace_back(ir + 6 , 0, 0);


			for(size_t i = 0 ; i < cameras.size(); i++){
				int ic = i * 7;
				ir = tripletListB.size();
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

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 7 + tie_points.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 7 + tie_points.size() * 3, cameras.size() * 7 + tie_points.size() * 3);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 7 + tie_points.size() * 3, 1);

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

			if(h_x.size() == cameras.size() * 7 + tie_points.size() * 3){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < cameras.size(); i++){
					QuaternionPose pose = pose_quaternion_from_affine_matrix(cameras[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.q0 += h_x[counter++];
					pose.q1 += h_x[counter++];
					pose.q2 += h_x[counter++];
					pose.q3 += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_quaternion(pose);
				}

				for(size_t i = 0; i < tie_points.size(); i++){
					tie_points[i].x() += h_x[counter++];
					tie_points[i].y() += h_x[counter++];
					tie_points[i].z() += h_x[counter++];
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
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
	std::cout << "t: optimize (Tait-Bryan)" << std::endl;
	std::cout << "r: optimize (Rodriguez)" << std::endl;
	std::cout << "q: optimize (Quaternion)" << std::endl;
	std::cout << "c: add noise to cameras" << std::endl;
	std::cout << "p: add noise to tie points" << std::endl;
}




