#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"
#include "point_to_line_tait_bryan_wc_jacobian.h"
#include "point_to_line_rodrigues_wc_jacobian.h"
#include "point_to_line_quaternion_wc_jacobian.h"
#include "quaternion_constraint_jacobian.h"

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -10.0;
float translate_x, translate_y = 0.0;

std::vector<Eigen::Affine3d> trajectory;
std::vector<Eigen::Affine3d> lines_global;
std::vector<std::vector<Eigen::Affine3d>> lines_local;


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
		pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 100;
		pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 100;
		pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 100;

		pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 10;
		pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 10;
		pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 10;
		lines_global.push_back(affine_matrix_from_pose_tait_bryan(pose));
	}

	lines_local.resize(trajectory.size());
	for(size_t i = 0 ; i < trajectory.size(); i++){
		Eigen::Affine3d m_inv = trajectory[i].inverse();

		for(size_t j = 0 ; j < lines_global.size(); j++){
			Eigen::Affine3d m_line_local = m_inv * lines_global[j];
			lines_local[i].push_back(m_line_local);
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
	glutCreateWindow("point_to_line");
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


	for(size_t i = 0 ; i < lines_local.size(); i++){
		Eigen::Affine3d &m = trajectory[i];

		for(size_t j = 0 ; j < lines_local[i].size(); j++){
			Eigen::Affine3d m_line_global = m * lines_local[i][j];

			glBegin(GL_LINES);
				glColor3f(1.0f, 0.0f, 0.0f);
				glVertex3f(m_line_global(0,3), m_line_global(1,3), m_line_global(2,3));
				glVertex3f(m_line_global(0,3) + m_line_global(0,0), m_line_global(1,3) + m_line_global(1,0), m_line_global(2,3) + m_line_global(2,0));

				glColor3f(0.0f, 1.0f, 0.0f);
				glVertex3f(m_line_global(0,3), m_line_global(1,3), m_line_global(2,3));
				glVertex3f(m_line_global(0,3) + m_line_global(0,1), m_line_global(1,3) + m_line_global(1,1), m_line_global(2,3) + m_line_global(2,1));

				glColor3f(0.0f, 0.0f, 1.0f);
				glVertex3f(m_line_global(0,3), m_line_global(1,3), m_line_global(2,3));
				glVertex3f(m_line_global(0,3) + m_line_global(0,2), m_line_global(1,3) + m_line_global(1,2), m_line_global(2,3) + m_line_global(2,2));
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
				pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.0;
				pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.0;
				pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.0;
				pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				trajectory[i] = trajectory[i] * affine_matrix_from_pose_tait_bryan(pose);
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < lines_local.size(); i++){
				for(size_t j = 0 ; j < lines_local.size(); j++){
					if(i != j){
						TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(trajectory[i]);

						for(size_t k = 0 ; k < lines_local[i].size(); k++){
							Eigen::Affine3d line_target_global = trajectory[j] * lines_local[j][k];

							Eigen::Vector3d vx(line_target_global(0,0), line_target_global(1,0), line_target_global(2,0));
							Eigen::Vector3d vy(line_target_global(0,1), line_target_global(1,1), line_target_global(2,1));
							Eigen::Vector3d point_on_target_line(line_target_global(0,3), line_target_global(1,3), line_target_global(2,3));

							Eigen::Vector3d point_source_local(lines_local[i][k](0,3), lines_local[i][k](1,3), lines_local[i][k](2,3));

							Eigen::Matrix<double, 2, 1> delta;
							point_to_line_tait_bryan_wc(delta,
									pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
									point_source_local.x(), point_source_local.y(), point_source_local.z(),
									point_on_target_line.x(), point_on_target_line.y(), point_on_target_line.z(),
									vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z());

							Eigen::Matrix<double, 2, 6> delta_jacobian;
							point_to_line_tait_bryan_wc_jacobian(delta_jacobian,
									pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
									point_source_local.x(), point_source_local.y(), point_source_local.z(),
									point_on_target_line.x(), point_on_target_line.y(), point_on_target_line.z(),
									vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z());

							int ir = tripletListB.size();

							for(int ii = 0; ii < 2; ii++){
								for(int jj = 0; jj < 6; jj++){
									int ic = i * 6;
									if(delta_jacobian(ii,jj) != 0.0){
										tripletListA.emplace_back(ir + ii, ic + jj , -delta_jacobian(ii,jj));
									}
								}
							}

							tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
							tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));

							tripletListB.emplace_back(ir    , 0,  delta(0,0));
							tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
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
				posetb.om += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.fi += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.ka += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				trajectory[i] = affine_matrix_from_pose_tait_bryan(posetb);
			}

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < lines_local.size(); i++){
				for(size_t j = 0 ; j < lines_local.size(); j++){
					if(i != j){
						RodriguesPose pose = pose_rodrigues_from_affine_matrix(trajectory[i]);

						for(size_t k = 0 ; k < lines_local[i].size(); k++){
							Eigen::Affine3d line_target_global = trajectory[j] * lines_local[j][k];

							Eigen::Vector3d vx(line_target_global(0,0), line_target_global(1,0), line_target_global(2,0));
							Eigen::Vector3d vy(line_target_global(0,1), line_target_global(1,1), line_target_global(2,1));
							Eigen::Vector3d point_on_target_line(line_target_global(0,3), line_target_global(1,3), line_target_global(2,3));

							Eigen::Vector3d point_source_local(lines_local[i][k](0,3), lines_local[i][k](1,3), lines_local[i][k](2,3));

							Eigen::Matrix<double, 2, 1> delta;
							point_to_line_rodrigues_wc(delta,
									pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
									point_source_local.x(), point_source_local.y(), point_source_local.z(),
									point_on_target_line.x(), point_on_target_line.y(), point_on_target_line.z(),
									vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z());

							Eigen::Matrix<double, 2, 6> delta_jacobian;
							point_to_line_rodrigues_wc_jacobian(delta_jacobian,
									pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
									point_source_local.x(), point_source_local.y(), point_source_local.z(),
									point_on_target_line.x(), point_on_target_line.y(), point_on_target_line.z(),
									vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z());

							int ir = tripletListB.size();

							for(int ii = 0; ii < 2; ii++){
								for(int jj = 0; jj < 6; jj++){
									int ic = i * 6;
									if(delta_jacobian(ii,jj) != 0.0){
										tripletListA.emplace_back(ir + ii, ic + jj , -delta_jacobian(ii,jj));
									}
								}
							}


							tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
							tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));;

							tripletListB.emplace_back(ir    , 0,  delta(0,0));
							tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
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

			for(size_t i = 0 ; i < lines_local.size(); i++){
				for(size_t j = 0 ; j < lines_local.size(); j++){
					if(i != j){
						QuaternionPose pose = pose_quaternion_from_affine_matrix(trajectory[i]);

						for(size_t k = 0 ; k < lines_local[i].size(); k++){
							Eigen::Affine3d line_target_global = trajectory[j] * lines_local[j][k];

							Eigen::Vector3d vx(line_target_global(0,0), line_target_global(1,0), line_target_global(2,0));
							Eigen::Vector3d vy(line_target_global(0,1), line_target_global(1,1), line_target_global(2,1));
							Eigen::Vector3d point_on_target_line(line_target_global(0,3), line_target_global(1,3), line_target_global(2,3));

							Eigen::Vector3d point_source_local(lines_local[i][k](0,3), lines_local[i][k](1,3), lines_local[i][k](2,3));

							Eigen::Matrix<double, 2, 1> delta;
							point_to_line_quaternion_wc(delta,
									pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
									point_source_local.x(), point_source_local.y(), point_source_local.z(),
									point_on_target_line.x(), point_on_target_line.y(), point_on_target_line.z(),
									vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z());

							Eigen::Matrix<double, 2, 7> delta_jacobian;
							point_to_line_quaternion_wc_jacobian(delta_jacobian,
									pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
									point_source_local.x(), point_source_local.y(), point_source_local.z(),
									point_on_target_line.x(), point_on_target_line.y(), point_on_target_line.z(),
									vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z());

							int ir = tripletListB.size();

							for(int ii = 0; ii < 2; ii++){
								for(int jj = 0; jj < 7; jj++){
									int ic = i * 7;
									if(delta_jacobian(ii,jj) != 0.0){
										tripletListA.emplace_back(ir + ii, ic + jj , -delta_jacobian(ii,jj));
									}
								}
							}

							tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
							tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));

							tripletListB.emplace_back(ir    , 0,  delta(0,0));
							tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
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
