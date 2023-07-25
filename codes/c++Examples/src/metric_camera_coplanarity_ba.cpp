#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"
#include "metric_camera_colinearity_tait_bryan_wc_jacobian.h"
#include "metric_camera_coplanarity_tait_bryan_wc_jacobian.h"
#include "metric_camera_coplanarity_rodrigues_wc_jacobian.h"
#include "metric_camera_coplanarity_quaternion_wc_jacobian.h"
#include "quaternion_constraint_jacobian.h"


struct KeyPoint{
	double ksi;
	double eta;
	int index_to_tie_point;
};

struct Camera{
	Eigen::Affine3d pose;
	std::vector<KeyPoint> key_points;
};

std::vector<Eigen::Vector3d> tie_points;
std::vector<Camera> cameras;
MetricCameraParams cam_params;

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

	cam_params.c = 5000;
	cam_params.eta_0 = 5000;
	cam_params.ksi_0 = 5000;

	for(size_t i = 0 ; i < 100; i++){
		Eigen::Vector3d p;
		p.x() = random(-100.0, 100.0);
		p.y() = random(-100.0, 100.0);
		p.z() = 100 + random(-5.0, 5.0);
		tie_points.push_back(p);
	}

	for(int i = -50 ; i < 50; i+=20){
		TaitBryanPose pose;
		pose.px = random(-1.0, 1.0);
		pose.py = i;
		pose.pz = random(-1.0, 1.0);
		pose.om = random(-0.1, 0.1);
		pose.fi = random(-0.1, 0.1);
		pose.ka = random(-0.1, 0.1);

		Camera c;
		c.pose = affine_matrix_from_pose_tait_bryan(pose);
		cameras.push_back(c);
	}


	for(size_t i = 0 ; i < cameras.size(); i++){
		for(size_t j = 0; j < tie_points.size(); j++){
			KeyPoint kp;
			kp.index_to_tie_point = j;

			TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
			metric_camera_colinearity_tait_bryan_wc(kp.ksi, kp.eta, cam_params.ksi_0, cam_params.eta_0, cam_params.c, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka, tie_points[j].x(), tie_points[j].y(), tie_points[j].z());

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
	glutCreateWindow("metric_camera_coplanarity_ba");
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
				pose.om = random(-0.1, 0.1);
				pose.fi = random(-0.1, 0.1);
				pose.ka = random(-0.1, 0.1);
				Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(pose);
				cameras[i].pose = cameras[i].pose * m;
			}
			break;
		}
		case 'p':{
			for(size_t i = 0; i < tie_points.size(); i++){
				tie_points[i].x() += random(-1.0, 1.0);
				tie_points[i].y() += random(-1.0, 1.0);
				tie_points[i].z() += random(-1.0, 1.0);
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras.size(); j++){
					for(size_t k = 0; k < cameras[i].key_points.size(); k++){
						if(i != j){
							TaitBryanPose pose_1 = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
							TaitBryanPose pose_2 = pose_tait_bryan_from_affine_matrix(cameras[j].pose);

							double px_1 = pose_1.px;
							double py_1 = pose_1.py;
							double pz_1 = pose_1.pz;
							double om_1 = pose_1.om;
							double fi_1 = pose_1.fi;
							double ka_1 = pose_1.ka;

							double px_2 = pose_2.px;
							double py_2 = pose_2.py;
							double pz_2 = pose_2.pz;
							double om_2 = pose_2.om;
							double fi_2 = pose_2.fi;
							double ka_2 = pose_2.ka;

							double ksi_1 = (cameras[i].key_points[k].ksi);
							double eta_1 = (cameras[i].key_points[k].eta);
							double ksi_2 = (cameras[j].key_points[k].ksi);
							double eta_2 = (cameras[j].key_points[k].eta);


							double delta;
							observation_equation_metric_camera_coplanarity_tait_bryan_wc(delta, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_1, eta_1, px_1, py_1, pz_1, om_1, fi_1, ka_1, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_2, eta_2, px_2, py_2, pz_2, om_2, fi_2, ka_2);


							Eigen::Matrix<double, 1, 12, Eigen::RowMajor> jacobian;
							observation_equation_metric_camera_coplanarity_tait_bryan_wc_jacobian(jacobian, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_1, eta_1, px_1, py_1, pz_1, om_1, fi_1, ka_1, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_2, eta_2, px_2, py_2, pz_2, om_2, fi_2, ka_2);

							int ir = tripletListB.size();
							int ic_1 = i * 3;
							int ic_2 = j * 3;

							tripletListA.emplace_back(ir     , ic_1 + 0, -jacobian(0,3));
							tripletListA.emplace_back(ir     , ic_1 + 1, -jacobian(0,4));
							tripletListA.emplace_back(ir     , ic_1 + 2, -jacobian(0,5));

							tripletListA.emplace_back(ir     , ic_2 + 0, -jacobian(0,9));
							tripletListA.emplace_back(ir     , ic_2 + 1, -jacobian(0,10));
							tripletListA.emplace_back(ir     , ic_2 + 2, -jacobian(0,11));

							tripletListP.emplace_back(ir    , ir    ,   cauchy(delta, 1));

							tripletListB.emplace_back(ir    , 0,  delta);
						}
					}
				}
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 3, cameras.size() * 3);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 3, 1);

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

			if(h_x.size() == cameras.size() * 3){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < cameras.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_tait_bryan(pose);
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}

			break;
		}
		case 'r':{
			for(size_t i = 0; i < cameras.size(); i++){
				TaitBryanPose posetb = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
				posetb.om += random(-0.000001, 0.000001);
				posetb.fi += random(-0.000001, 0.000001);
				posetb.ka += random(-0.000001, 0.000001);
				cameras[i].pose = affine_matrix_from_pose_tait_bryan(posetb);
			}

			std::vector<Eigen::Triplet<double>> tripletListA;
				std::vector<Eigen::Triplet<double>> tripletListP;
				std::vector<Eigen::Triplet<double>> tripletListB;

				for(size_t i = 0; i < cameras.size(); i++){
					for(size_t j = 0; j < cameras.size(); j++){
						for(size_t k = 0; k < cameras[i].key_points.size(); k++){
							if(i != j){
								RodriguesPose pose_1 = pose_rodrigues_from_affine_matrix(cameras[i].pose);
								RodriguesPose pose_2 = pose_rodrigues_from_affine_matrix(cameras[j].pose);

								double px_1 = pose_1.px;
								double py_1 = pose_1.py;
								double pz_1 = pose_1.pz;
								double sx_1 = pose_1.sx;
								double sy_1 = pose_1.sy;
								double sz_1 = pose_1.sz;

								double px_2 = pose_2.px;
								double py_2 = pose_2.py;
								double pz_2 = pose_2.pz;
								double sx_2 = pose_2.sx;
								double sy_2 = pose_2.sy;
								double sz_2 = pose_2.sz;

								double ksi_1 = (cameras[i].key_points[k].ksi);
								double eta_1 = (cameras[i].key_points[k].eta);
								double ksi_2 = (cameras[j].key_points[k].ksi);
								double eta_2 = (cameras[j].key_points[k].eta);


								double delta;
								observation_equation_metric_camera_coplanarity_rodrigues_wc(delta, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_1, eta_1, px_1, py_1, pz_1, sx_1, sy_1, sz_1, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_2, eta_2, px_2, py_2, pz_2, sx_2, sy_2, sz_2);

								Eigen::Matrix<double, 1, 12, Eigen::RowMajor> jacobian;
								observation_equation_metric_camera_coplanarity_rodrigues_wc_jacobian(jacobian, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_1, eta_1, px_1, py_1, pz_1, sx_1, sy_1, sz_1, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_2, eta_2, px_2, py_2, pz_2, sx_2, sy_2, sz_2);

								int ir = tripletListB.size();
								int ic_1 = i * 3;
								int ic_2 = j * 3;

								tripletListA.emplace_back(ir     , ic_1 + 0, -jacobian(0,3));
								tripletListA.emplace_back(ir     , ic_1 + 1, -jacobian(0,4));
								tripletListA.emplace_back(ir     , ic_1 + 2, -jacobian(0,5));

								tripletListA.emplace_back(ir     , ic_2 + 0, -jacobian(0,9));
								tripletListA.emplace_back(ir     , ic_2 + 1, -jacobian(0,10));
								tripletListA.emplace_back(ir     , ic_2 + 2, -jacobian(0,11));

								tripletListP.emplace_back(ir    , ir    ,   cauchy(delta, 1));

								tripletListB.emplace_back(ir    , 0,  delta);
							}
						}
					}
				}

				Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 3);
				Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
				Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

				matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
				matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
				matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

				Eigen::SparseMatrix<double> AtPA(cameras.size() * 3, cameras.size() * 3);
				Eigen::SparseMatrix<double> AtPB(cameras.size() * 3, 1);

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

				if(h_x.size() == cameras.size() * 3){
					for(size_t i = 0 ; i < h_x.size(); i++){
						std::cout << h_x[i] << std::endl;
					}
					std::cout << "AtPA=AtPB SOLVED" << std::endl;
					std::cout << "update" << std::endl;

					int counter = 0;

					for(size_t i = 0; i < cameras.size(); i++){
						RodriguesPose pose = pose_rodrigues_from_affine_matrix(cameras[i].pose);
						pose.sx += h_x[counter++];
						pose.sy += h_x[counter++];
						pose.sz += h_x[counter++];

						cameras[i].pose = affine_matrix_from_pose_rodrigues(pose);
					}
				}else{
					std::cout << "AtPA=AtPB FAILED" << std::endl;
				}
			break;
		}
		case 'q':{
			for(size_t i = 0; i < cameras.size(); i++){
				TaitBryanPose posetb = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
				posetb.om += random(-0.000001, 0.000001);
				posetb.fi += random(-0.000001, 0.000001);
				posetb.ka += random(-0.000001, 0.000001);
				cameras[i].pose = affine_matrix_from_pose_tait_bryan(posetb);
			}

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras.size(); j++){
					for(size_t k = 0; k < cameras[i].key_points.size(); k++){
						if(i != j){
							QuaternionPose pose_1 = pose_quaternion_from_affine_matrix(cameras[i].pose);
							QuaternionPose pose_2 = pose_quaternion_from_affine_matrix(cameras[j].pose);

							double px_1 = pose_1.px;
							double py_1 = pose_1.py;
							double pz_1 = pose_1.pz;
							double q0_1 = pose_1.q0;
							double q1_1 = pose_1.q1;
							double q2_1 = pose_1.q2;
							double q3_1 = pose_1.q3;

							double px_2 = pose_2.px;
							double py_2 = pose_2.py;
							double pz_2 = pose_2.pz;
							double q0_2 = pose_2.q0;
							double q1_2 = pose_2.q1;
							double q2_2 = pose_2.q2;
							double q3_2 = pose_2.q3;

							double ksi_1 = (cameras[i].key_points[k].ksi);
							double eta_1 = (cameras[i].key_points[k].eta);
							double ksi_2 = (cameras[j].key_points[k].ksi);
							double eta_2 = (cameras[j].key_points[k].eta);


							double delta;
							observation_equation_metric_camera_coplanarity_quaternion_wc(delta, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_1, eta_1, px_1, py_1, pz_1, q0_1, q1_1, q2_1,q3_1, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_2, eta_2, px_2, py_2, pz_2, q0_2, q1_2, q2_2, q3_2);

							Eigen::Matrix<double, 1, 14, Eigen::RowMajor> jacobian;
							observation_equation_metric_camera_coplanarity_quaternion_wc_jacobian(jacobian, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_1, eta_1, px_1, py_1, pz_1, q0_1, q1_1, q2_1, q3_1, cam_params.ksi_0, cam_params.eta_0, cam_params.c, ksi_2, eta_2, px_2, py_2, pz_2, q0_2, q1_2, q2_2, q3_2);

							int ir = tripletListB.size();
							int ic_1 = i * 4;
							int ic_2 = j * 4;

							tripletListA.emplace_back(ir     , ic_1 + 0, -jacobian(0,3));
							tripletListA.emplace_back(ir     , ic_1 + 1, -jacobian(0,4));
							tripletListA.emplace_back(ir     , ic_1 + 2, -jacobian(0,5));
							tripletListA.emplace_back(ir     , ic_1 + 3, -jacobian(0,6));

							tripletListA.emplace_back(ir     , ic_2 + 0, -jacobian(0,10));
							tripletListA.emplace_back(ir     , ic_2 + 1, -jacobian(0,11));
							tripletListA.emplace_back(ir     , ic_2 + 2, -jacobian(0,12));
							tripletListA.emplace_back(ir     , ic_2 + 3, -jacobian(0,13));

							tripletListP.emplace_back(ir    , ir    ,   cauchy(delta, 1));

							tripletListB.emplace_back(ir    , 0,  delta);
						}
					}
				}
			}

			for(size_t i = 0 ; i < cameras.size(); i++){
				int ic = i * 4;
				int ir = tripletListB.size();
				QuaternionPose pose = pose_quaternion_from_affine_matrix(cameras[i].pose);

				double delta;
				quaternion_constraint(delta, pose.q0, pose.q1, pose.q2, pose.q3);

				Eigen::Matrix<double, 1, 4> jacobian;
				quaternion_constraint_jacobian(jacobian, pose.q0, pose.q1, pose.q2, pose.q3);

				tripletListA.emplace_back(ir, ic + 0 , -jacobian(0,0));
				tripletListA.emplace_back(ir, ic + 1 , -jacobian(0,1));
				tripletListA.emplace_back(ir, ic + 2 , -jacobian(0,2));
				tripletListA.emplace_back(ir, ic + 3 , -jacobian(0,3));

				tripletListP.emplace_back(ir, ir, 1000000.0);

				tripletListB.emplace_back(ir, 0, delta);
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 4);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 4, cameras.size() * 4);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 4, 1);

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

			if(h_x.size() == cameras.size() * 4){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < cameras.size(); i++){
					QuaternionPose pose = pose_quaternion_from_affine_matrix(cameras[i].pose);
					pose.q0 += h_x[counter++];
					pose.q1 += h_x[counter++];
					pose.q2 += h_x[counter++];
					pose.q3 += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_quaternion(pose);
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




