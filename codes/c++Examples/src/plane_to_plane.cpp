#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"
#include "plane_to_plane_tait_bryan_wc_jacobian.h"
#include "plane_to_plane_rodrigues_wc_jacobian.h"
#include "plane_to_plane_quaternion_wc_jacobian.h"
#include "quaternion_constraint_jacobian.h"

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -10.0;
float translate_x, translate_y = 0.0;

std::vector<Eigen::Affine3d> trajectory;
std::vector<Eigen::Affine3d> planes_global;
std::vector<std::vector<Eigen::Affine3d>> planes_local;

struct Plane{
	double a;
	double b;
	double c;
	double d;
};

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

int main(int argc, char *argv[]){
	TaitBryanPose pose;
	pose.px = -40;
	pose.py = 0;
	pose.pz = 0;
	pose.om = 0;
	pose.fi = 0;
	pose.ka = 0;
	trajectory.push_back(affine_matrix_from_pose_tait_bryan(pose));

	pose.px = -20;
	pose.py = 0;
	pose.pz = 0;
	pose.om = 0;
	pose.fi = 0;
	pose.ka = 0;
	trajectory.push_back(affine_matrix_from_pose_tait_bryan(pose));

	pose.px = -0;
	pose.py = 0;
	pose.pz = 0;
	pose.om = 0;
	pose.fi = 0;
	pose.ka = 0;
	trajectory.push_back(affine_matrix_from_pose_tait_bryan(pose));

	pose.px = 20;
	pose.py = 0;
	pose.pz = 0;
	pose.om = 0;
	pose.fi = 0;
	pose.ka = 0;
	trajectory.push_back(affine_matrix_from_pose_tait_bryan(pose));

	pose.px = 40;
	pose.py = 0;
	pose.pz = 0;
	pose.om = 0;
	pose.fi = 0;
	pose.ka = 0;
	trajectory.push_back(affine_matrix_from_pose_tait_bryan(pose));

	for(size_t i = 0 ; i < 100; i++){
		pose.px = random(-100.0, 100.0);
		pose.py = random(-100.0, 100.0);
		pose.pz = random(-100.0, 100.0);

		pose.om = random(-10.0, 10.0);
		pose.fi = random(-10.0, 10.0);
		pose.ka = random(-10.0, 10.0);
		planes_global.push_back(affine_matrix_from_pose_tait_bryan(pose));
	}

	planes_local.resize(trajectory.size());
	for(size_t i = 0 ; i < trajectory.size(); i++){
		Eigen::Affine3d m_inv = trajectory[i].inverse();

		for(size_t j = 0 ; j < planes_global.size(); j++){
			Eigen::Affine3d m_plane_local = m_inv * planes_global[j];
			planes_local[i].push_back(m_plane_local);
		}
	}

	if (false == initGL(&argc, argv)) {
		return 4;
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
	glutCreateWindow("plane_to_plane");
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
	gluPerspective(60.0, (GLfloat) window_width / (GLfloat) window_height, 0.01, 10000.0);
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

	for(size_t i = 0 ; i < trajectory.size(); i++){
		Eigen::Affine3d &m = trajectory[i];

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
	for(size_t i = 0 ; i < trajectory.size(); i++){
		Eigen::Affine3d &m = trajectory[i];
		glVertex3f(m(0,3), m(1,3), m(2,3));
	}
	glEnd();


	for(size_t i = 0 ; i < planes_local.size(); i++){
		Eigen::Affine3d &m = trajectory[i];

		for(size_t j = 0 ; j < planes_local[i].size(); j++){
			Eigen::Affine3d m_plane_global = m * planes_local[i][j];

			glBegin(GL_LINES);
				glColor3f(1.0f, 0.0f, 0.0f);
				glVertex3f(m_plane_global(0,3), m_plane_global(1,3), m_plane_global(2,3));
				glVertex3f(m_plane_global(0,3) + m_plane_global(0,0), m_plane_global(1,3) + m_plane_global(1,0), m_plane_global(2,3) + m_plane_global(2,0));

				glColor3f(0.0f, 1.0f, 0.0f);
				glVertex3f(m_plane_global(0,3), m_plane_global(1,3), m_plane_global(2,3));
				glVertex3f(m_plane_global(0,3) + m_plane_global(0,1), m_plane_global(1,3) + m_plane_global(1,1), m_plane_global(2,3) + m_plane_global(2,1));

				glColor3f(0.0f, 0.0f, 1.0f);
				glVertex3f(m_plane_global(0,3), m_plane_global(1,3), m_plane_global(2,3));
				glVertex3f(m_plane_global(0,3) + m_plane_global(0,2), m_plane_global(1,3) + m_plane_global(1,2), m_plane_global(2,3) + m_plane_global(2,2));
			glEnd();
		}
	}

	glutSwapBuffers();
}


void keyboard(unsigned char key, int /*x*/, int /*y*/) {
	switch (key) {
		case (27): {
			glutDestroyWindow(glutGetWindow());
			return;
		}
		case 'n':{
			for(size_t i = 0 ; i  < trajectory.size(); i++){
				TaitBryanPose pose;
				pose.px = random(-1.0, 1.0);
				pose.py = random(-1.0, 1.0);
				pose.pz = random(-1.0, 1.0);
				pose.om = random(-0.1, 0.1);
				pose.fi = random(-0.1, 0.1);
				pose.ka = random(-0.1, 0.1);
				trajectory[i] = trajectory[i] * affine_matrix_from_pose_tait_bryan(pose);
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < planes_local.size(); i++){
				for(size_t j = i+1 ; j < planes_local.size(); j++){
					if(i != j){
						int index_pose_from = i;
						int index_pose_to = j;

						TaitBryanPose pose_from = pose_tait_bryan_from_affine_matrix(trajectory[i]);
						TaitBryanPose pose_to = pose_tait_bryan_from_affine_matrix(trajectory[j]);

						for(size_t k = 0 ; k < planes_local[i].size(); k++){
							Plane plane_from;
							plane_from.a = planes_local[i][k](0,2);
							plane_from.b = planes_local[i][k](1,2);
							plane_from.c = planes_local[i][k](2,2);
							plane_from.d = -plane_from.a * planes_local[i][k](0,3) - plane_from.b * planes_local[i][k](1,3) - plane_from.c * planes_local[i][k](2,3);

							Plane plane_to;
							plane_to.a = planes_local[j][k](0,2);
							plane_to.b = planes_local[j][k](1,2);
							plane_to.c = planes_local[j][k](2,2);
							plane_to.d = -plane_from.a * planes_local[j][k](0,3) - plane_from.b * planes_local[j][k](1,3) - plane_from.c * planes_local[j][k](2,3);

							Eigen::Matrix<double, 4, 1> delta;
							plane_to_plane_tait_bryan_wc(delta,
									pose_from.px, pose_from.py, pose_from.pz, pose_from.om, pose_from.fi, pose_from.ka,
									pose_to.px, pose_to.py, pose_to.pz, pose_to.om, pose_to.fi, pose_to.ka,
									plane_from.a, plane_from.b, plane_from.c, plane_from.d,
									plane_to.a, plane_to.b, plane_to.c, plane_to.d);

							Eigen::Matrix<double, 4, 12, Eigen::RowMajor> delta_jacobian;
							plane_to_plane_tait_bryan_wc_jacobian(delta_jacobian,
									pose_from.px, pose_from.py, pose_from.pz, pose_from.om, pose_from.fi, pose_from.ka,
									pose_to.px, pose_to.py, pose_to.pz, pose_to.om, pose_to.fi, pose_to.ka,
									plane_from.a, plane_from.b, plane_from.c, plane_from.d,
									plane_to.a, plane_to.b, plane_to.c, plane_to.d);

							int ir = tripletListB.size();

							for(int ii = 0; ii < 4; ii++){
								for(int jj = 0; jj < 6; jj++){
									int ic = i * 6;
									if(delta_jacobian(ii,jj) != 0.0){
										tripletListA.emplace_back(ir + ii, ic + jj , -delta_jacobian(ii,jj));
									}
									ic = j * 6;
									if(delta_jacobian(ii,jj + 6) != 0.0){
										tripletListA.emplace_back(ir + ii, ic + jj , -delta_jacobian(ii,jj + 6));
									}
								}
							}

							tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
							tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));
							tripletListP.emplace_back(ir + 2, ir + 2,  cauchy(delta(2,0), 1));
							tripletListP.emplace_back(ir + 3, ir + 3,  cauchy(delta(3,0), 1));

							tripletListB.emplace_back(ir    , 0,  delta(0,0));
							tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
							tripletListB.emplace_back(ir + 2, 0,  delta(2,0));
							tripletListB.emplace_back(ir + 3, 0,  delta(3,0));
						}
					}
				}
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);
			tripletListA.emplace_back(ir + 3 , 3, 1);
			tripletListA.emplace_back(ir + 4 , 4, 1);
			tripletListA.emplace_back(ir + 5 , 5, 1);

			tripletListP.emplace_back(ir     , ir,     10000000000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 10000000000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 10000000000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 10000000000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 10000000000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 10000000000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);

			Eigen::SparseMatrix<double> matA(tripletListB.size(), trajectory.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(trajectory.size() * 6, trajectory.size() * 6);
			Eigen::SparseMatrix<double> AtPB(trajectory.size() * 6, 1);

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

			if(h_x.size() == trajectory.size() * 6){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < trajectory.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(trajectory[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];

					trajectory[i] = affine_matrix_from_pose_tait_bryan(pose);
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}

			break;
		}
		case 'r':{
			for(size_t i = 0; i < trajectory.size(); i++){
				TaitBryanPose posetb = pose_tait_bryan_from_affine_matrix(trajectory[i]);
				posetb.om += random(-0.000001, 0.000001);
				posetb.fi += random(-0.000001, 0.000001);
				posetb.ka += random(-0.000001, 0.000001);
				trajectory[i] = affine_matrix_from_pose_tait_bryan(posetb);
			}

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < planes_local.size(); i++){
				for(size_t j = i+1 ; j < planes_local.size(); j++){
					if(i != j){
						int index_pose_from = i;
						int index_pose_to = j;

						RodriguesPose pose_from = pose_rodrigues_from_affine_matrix(trajectory[i]);
						RodriguesPose pose_to = pose_rodrigues_from_affine_matrix(trajectory[j]);

						for(size_t k = 0 ; k < planes_local[i].size(); k++){
							Plane plane_from;
							plane_from.a = planes_local[i][k](0,2);
							plane_from.b = planes_local[i][k](1,2);
							plane_from.c = planes_local[i][k](2,2);
							plane_from.d = -plane_from.a * planes_local[i][k](0,3) - plane_from.b * planes_local[i][k](1,3) - plane_from.c * planes_local[i][k](2,3);

							Plane plane_to;
							plane_to.a = planes_local[j][k](0,2);
							plane_to.b = planes_local[j][k](1,2);
							plane_to.c = planes_local[j][k](2,2);
							plane_to.d = -plane_from.a * planes_local[j][k](0,3) - plane_from.b * planes_local[j][k](1,3) - plane_from.c * planes_local[j][k](2,3);

							Eigen::Matrix<double, 4, 1> delta;
							plane_to_plane_rodrigues_wc(delta,
									pose_from.px, pose_from.py, pose_from.pz, pose_from.sx, pose_from.sy, pose_from.sz,
									pose_to.px, pose_to.py, pose_to.pz, pose_to.sx, pose_to.sy, pose_to.sz,
									plane_from.a, plane_from.b, plane_from.c, plane_from.d,
									plane_to.a, plane_to.b, plane_to.c, plane_to.d);

							Eigen::Matrix<double, 4, 12, Eigen::RowMajor> delta_jacobian;
							plane_to_plane_rodrigues_wc_jacobian(delta_jacobian,
									pose_from.px, pose_from.py, pose_from.pz, pose_from.sx, pose_from.sy, pose_from.sz,
									pose_to.px, pose_to.py, pose_to.pz, pose_to.sx, pose_to.sy, pose_to.sz,
									plane_from.a, plane_from.b, plane_from.c, plane_from.d,
									plane_to.a, plane_to.b, plane_to.c, plane_to.d);

							int ir = tripletListB.size();

							for(int ii = 0; ii < 4; ii++){
								for(int jj = 0; jj < 6; jj++){
									int ic = i * 6;
									if(delta_jacobian(ii,jj) != 0.0){
										tripletListA.emplace_back(ir + ii, ic + jj , -delta_jacobian(ii,jj));
									}
									ic = j * 6;
									if(delta_jacobian(ii,jj + 6) != 0.0){
										tripletListA.emplace_back(ir + ii, ic + jj , -delta_jacobian(ii,jj + 6));
									}
								}
							}

							tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
							tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));
							tripletListP.emplace_back(ir + 2, ir + 2,  cauchy(delta(2,0), 1));
							tripletListP.emplace_back(ir + 3, ir + 3,  cauchy(delta(3,0), 1));

							tripletListB.emplace_back(ir    , 0,  delta(0,0));
							tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
							tripletListB.emplace_back(ir + 2, 0,  delta(2,0));
							tripletListB.emplace_back(ir + 3, 0,  delta(3,0));
						}
					}
				}
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);
			tripletListA.emplace_back(ir + 3 , 3, 1);
			tripletListA.emplace_back(ir + 4 , 4, 1);
			tripletListA.emplace_back(ir + 5 , 5, 1);

			tripletListP.emplace_back(ir     , ir,     10000000000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 10000000000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 10000000000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 10000000000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 10000000000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 10000000000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);

			Eigen::SparseMatrix<double> matA(tripletListB.size(), trajectory.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(trajectory.size() * 6, trajectory.size() * 6);
			Eigen::SparseMatrix<double> AtPB(trajectory.size() * 6, 1);

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

			if(h_x.size() == trajectory.size() * 6){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < trajectory.size(); i++){
					RodriguesPose pose = pose_rodrigues_from_affine_matrix(trajectory[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.sx += h_x[counter++];
					pose.sy += h_x[counter++];
					pose.sz += h_x[counter++];

					trajectory[i] = affine_matrix_from_pose_rodrigues(pose);
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'q':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < planes_local.size(); i++){
				for(size_t j = i+1 ; j < planes_local.size(); j++){
					if(i != j){
						int index_pose_from = i;
						int index_pose_to = j;

						QuaternionPose pose_from = pose_quaternion_from_affine_matrix(trajectory[i]);
						QuaternionPose pose_to = pose_quaternion_from_affine_matrix(trajectory[j]);

						for(size_t k = 0 ; k < planes_local[i].size(); k++){
							Plane plane_from;
							plane_from.a = planes_local[i][k](0,2);
							plane_from.b = planes_local[i][k](1,2);
							plane_from.c = planes_local[i][k](2,2);
							plane_from.d = -plane_from.a * planes_local[i][k](0,3) - plane_from.b * planes_local[i][k](1,3) - plane_from.c * planes_local[i][k](2,3);

							Plane plane_to;
							plane_to.a = planes_local[j][k](0,2);
							plane_to.b = planes_local[j][k](1,2);
							plane_to.c = planes_local[j][k](2,2);
							plane_to.d = -plane_from.a * planes_local[j][k](0,3) - plane_from.b * planes_local[j][k](1,3) - plane_from.c * planes_local[j][k](2,3);

							Eigen::Matrix<double, 4, 1> delta;
							plane_to_plane_quaternion_wc(delta,
									pose_from.px, pose_from.py, pose_from.pz, pose_from.q0, pose_from.q1, pose_from.q2, pose_from.q3,
									pose_to.px, pose_to.py, pose_to.pz, pose_to.q0, pose_to.q1, pose_to.q2, pose_to.q3,
									plane_from.a, plane_from.b, plane_from.c, plane_from.d,
									plane_to.a, plane_to.b, plane_to.c, plane_to.d);

							Eigen::Matrix<double, 4, 14, Eigen::RowMajor> delta_jacobian;
							plane_to_plane_quaternion_wc_jacobian(delta_jacobian,
									pose_from.px, pose_from.py, pose_from.pz, pose_from.q0, pose_from.q1, pose_from.q2, pose_from.q3,
									pose_to.px, pose_to.py, pose_to.pz, pose_to.q0, pose_to.q1, pose_to.q2, pose_to.q3,
									plane_from.a, plane_from.b, plane_from.c, plane_from.d,
									plane_to.a, plane_to.b, plane_to.c, plane_to.d);

							int ir = tripletListB.size();

							for(int ii = 0; ii < 4; ii++){
								for(int jj = 0; jj < 7; jj++){
									int ic = i * 7;
									if(delta_jacobian(ii,jj) != 0.0){
										tripletListA.emplace_back(ir + ii, ic + jj , -delta_jacobian(ii,jj));
									}
									ic = j * 7;
									if(delta_jacobian(ii,jj + 7) != 0.0){
										tripletListA.emplace_back(ir + ii, ic + jj , -delta_jacobian(ii,jj + 7));
									}
								}
							}

							tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
							tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));
							tripletListP.emplace_back(ir + 2, ir + 2,  cauchy(delta(2,0), 1));
							tripletListP.emplace_back(ir + 3, ir + 3,  cauchy(delta(3,0), 1));

							tripletListB.emplace_back(ir    , 0,  delta(0,0));
							tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
							tripletListB.emplace_back(ir + 2, 0,  delta(2,0));
							tripletListB.emplace_back(ir + 3, 0,  delta(3,0));
						}
					}
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

			tripletListP.emplace_back(ir     , ir,     10000000000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 10000000000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 10000000000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 10000000000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 10000000000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 10000000000000);
			tripletListP.emplace_back(ir + 6 , ir + 6, 10000000000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);
			tripletListB.emplace_back(ir + 6 , 0, 0);

			for(size_t i = 0 ; i < trajectory.size(); i++){
				int ic = i * 7;
				ir = tripletListB.size();
				QuaternionPose pose = pose_quaternion_from_affine_matrix(trajectory[i]);

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

			Eigen::SparseMatrix<double> matA(tripletListB.size(), trajectory.size() * 7);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(trajectory.size() * 7, trajectory.size() * 7);
			Eigen::SparseMatrix<double> AtPB(trajectory.size() * 7, 1);

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

			if(h_x.size() == trajectory.size() * 7){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < trajectory.size(); i++){
					QuaternionPose pose = pose_quaternion_from_affine_matrix(trajectory[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.q0 += h_x[counter++];
					pose.q1 += h_x[counter++];
					pose.q2 += h_x[counter++];
					pose.q3 += h_x[counter++];
					trajectory[i] = affine_matrix_from_pose_quaternion(pose);
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
	std::cout << "n: add noise to trajectory" << std::endl;
	std::cout << "t: optimize Tait-Bryan" << std::endl;
	std::cout << "r: optimize Rodrigues" << std::endl;
	std::cout << "q: optimize Quaternion" << std::endl;
}
