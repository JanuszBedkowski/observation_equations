#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "relative_pose_tait_bryan_wc_jacobian.h"
#include "relative_pose_rodrigues_wc_jacobian.h"
#include "relative_pose_quaternion_wc_jacobian.h"
#include "quaternion_constraint_jacobian.h"
#include "relative_pose_wc_jacobian.h"
#include "relative_pose_2_tait_bryan_wc_jacobian.h"
#include "relative_pose_2_rodrigues_wc_jacobian.h"
#include "relative_pose_2_quaternion_wc_jacobian.h"

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -10.0;
float translate_x, translate_y = 0.0;

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

std::vector<Eigen::Affine3d> m_poses;
std::vector<Eigen::Affine3d> m_poses_desired;

std::vector<std::pair<int, int>> odo_edges;
std::vector<std::pair<int, int>> loop_edges;

int main(int argc, char *argv[]){

	for(size_t i = 0 ; i < 100; i++){
		TaitBryanPose p;
		p.px = i;
		p.py = -1;
		p.pz = 0.0;
		p.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
		p.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
		p.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;

		Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(p);
		m_poses.push_back(m);
	}
	for(size_t i = 0 ; i < 100; i++){
		TaitBryanPose p;
		p.px = i;
		p.py = 1;
		p.pz = 0.0;
		p.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
		p.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
		p.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;

		Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(p);
		m_poses.push_back(m);
	}
	m_poses_desired = m_poses;

	for(size_t i = 1; i < 100; i++){
		odo_edges.emplace_back(i-1,i);
	}

	for(size_t i = 101; i < 200; i++){
		odo_edges.emplace_back(i-1,i);
	}

	for(size_t i = 0; i < 100; i+=10){
		loop_edges.emplace_back(i,i+100);
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
	glutCreateWindow("relative_pose");
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

	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(1.0f, 0.0f, 0.0f);

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 1.0f, 0.0f);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 1.0f);
	glEnd();

	glColor3f(1,0,0);
	glBegin(GL_LINES);
	for(size_t i = 0; i < odo_edges.size(); i++){
		glVertex3f(m_poses[odo_edges[i].first](0,3), m_poses[odo_edges[i].first](1,3), m_poses[odo_edges[i].first](2,3) );
		glVertex3f(m_poses[odo_edges[i].second](0,3), m_poses[odo_edges[i].second](1,3), m_poses[odo_edges[i].second](2,3) );
	}
	glEnd();

	glColor3f(0,1,0);
	glBegin(GL_LINES);
	for(size_t i = 0; i < loop_edges.size(); i++){
		glVertex3f(m_poses[loop_edges[i].first](0,3), m_poses[loop_edges[i].first](1,3), m_poses[loop_edges[i].first](2,3) );
		glVertex3f(m_poses[loop_edges[i].second](0,3), m_poses[loop_edges[i].second](1,3), m_poses[loop_edges[i].second](2,3) );
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
		case 'n':{
			for(size_t i = 0 ; i < m_poses.size(); i++){
				TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m_poses[i]);
				pose.px += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.py += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.pz += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.om += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
				pose.fi += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
				pose.ka += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
				m_poses[i] = affine_matrix_from_pose_tait_bryan(pose);
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<TaitBryanPose> poses;
			std::vector<TaitBryanPose> poses_desired;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_tait_bryan_from_affine_matrix(m_poses[i]));
			}
			for(size_t i = 0 ; i < m_poses_desired.size(); i++){
				poses_desired.push_back(pose_tait_bryan_from_affine_matrix(m_poses_desired[i]));
			}

			for(size_t i = 0 ; i < odo_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_odo;
				relative_pose_tait_bryan_wc_case1(relative_pose_measurement_odo,
						poses_desired[odo_edges[i].first].px,
						poses_desired[odo_edges[i].first].py,
						poses_desired[odo_edges[i].first].pz,
						poses_desired[odo_edges[i].first].om,
						poses_desired[odo_edges[i].first].fi,
						poses_desired[odo_edges[i].first].ka,
						poses_desired[odo_edges[i].second].px,
						poses_desired[odo_edges[i].second].py,
						poses_desired[odo_edges[i].second].pz,
						poses_desired[odo_edges[i].second].om,
						poses_desired[odo_edges[i].second].fi,
						poses_desired[odo_edges[i].second].ka);

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_tait_bryan_wc_case1(
						delta,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].om,
						poses[odo_edges[i].first].fi,
						poses[odo_edges[i].first].ka,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].om,
						poses[odo_edges[i].second].fi,
						poses[odo_edges[i].second].ka,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_tait_bryan_wc_case1_jacobian(jacobian,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].om,
						poses[odo_edges[i].first].fi,
						poses[odo_edges[i].first].ka,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].om,
						poses[odo_edges[i].second].fi,
						poses[odo_edges[i].second].ka);

				int ir = tripletListB.size();

				int ic_1 = odo_edges[i].first * 6;
				int ic_2 = odo_edges[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
			}

			for(size_t i = 0 ; i < loop_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_loop;
				relative_pose_tait_bryan_wc_case1(relative_pose_measurement_loop,
						poses_desired[loop_edges[i].first].px,
						poses_desired[loop_edges[i].first].py,
						poses_desired[loop_edges[i].first].pz,
						poses_desired[loop_edges[i].first].om,
						poses_desired[loop_edges[i].first].fi,
						poses_desired[loop_edges[i].first].ka,
						poses_desired[loop_edges[i].second].px,
						poses_desired[loop_edges[i].second].py,
						poses_desired[loop_edges[i].second].pz,
						poses_desired[loop_edges[i].second].om,
						poses_desired[loop_edges[i].second].fi,
						poses_desired[loop_edges[i].second].ka);

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_tait_bryan_wc_case1(
						delta,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].om,
						poses[loop_edges[i].first].fi,
						poses[loop_edges[i].first].ka,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].om,
						poses[loop_edges[i].second].fi,
						poses[loop_edges[i].second].ka,
						relative_pose_measurement_loop(0,0),
						relative_pose_measurement_loop(1,0),
						relative_pose_measurement_loop(2,0),
						relative_pose_measurement_loop(3,0),
						relative_pose_measurement_loop(4,0),
						relative_pose_measurement_loop(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_tait_bryan_wc_case1_jacobian(jacobian,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].om,
						poses[loop_edges[i].first].fi,
						poses[loop_edges[i].first].ka,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].om,
						poses[loop_edges[i].second].fi,
						poses[loop_edges[i].second].ka);

				int ir = tripletListB.size();

				int ic_1 = loop_edges[i].first * 6;
				int ic_2 = loop_edges[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
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


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 6 , m_poses.size() * 6);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 6 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 6 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];
					m_poses[i] = affine_matrix_from_pose_tait_bryan(pose);
				}
				std::cout << "optimizing with tait bryan finished" << std::endl;
			}else{
				std::cout << "optimizing with tait bryan FAILED" << std::endl;
			}

			break;
		}

		case 'r':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<RodriguesPose> poses;
			std::vector<RodriguesPose> poses_desired;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_rodrigues_from_affine_matrix(m_poses[i]));
			}
			for(size_t i = 0 ; i < m_poses_desired.size(); i++){
				poses_desired.push_back(pose_rodrigues_from_affine_matrix(m_poses_desired[i]));
			}

			for(size_t i = 0 ; i < odo_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_odo;
				relative_pose_rodrigues_wc(relative_pose_measurement_odo,
						poses_desired[odo_edges[i].first].px,
						poses_desired[odo_edges[i].first].py,
						poses_desired[odo_edges[i].first].pz,
						poses_desired[odo_edges[i].first].sx,
						poses_desired[odo_edges[i].first].sy,
						poses_desired[odo_edges[i].first].sz,
						poses_desired[odo_edges[i].second].px,
						poses_desired[odo_edges[i].second].py,
						poses_desired[odo_edges[i].second].pz,
						poses_desired[odo_edges[i].second].sx,
						poses_desired[odo_edges[i].second].sy,
						poses_desired[odo_edges[i].second].sz);

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_rodrigues_wc(
						delta,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].sx,
						poses[odo_edges[i].first].sy,
						poses[odo_edges[i].first].sz,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].sx,
						poses[odo_edges[i].second].sy,
						poses[odo_edges[i].second].sz,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_rodrigues_wc_jacobian(jacobian,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].sx,
						poses[odo_edges[i].first].sy,
						poses[odo_edges[i].first].sz,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].sx,
						poses[odo_edges[i].second].sy,
						poses[odo_edges[i].second].sz);

				int ir = tripletListB.size();

				int ic_1 = odo_edges[i].first * 6;
				int ic_2 = odo_edges[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
			}

			for(size_t i = 0 ; i < loop_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_loop;
				relative_pose_rodrigues_wc(relative_pose_measurement_loop,
						poses_desired[loop_edges[i].first].px,
						poses_desired[loop_edges[i].first].py,
						poses_desired[loop_edges[i].first].pz,
						poses_desired[loop_edges[i].first].sx,
						poses_desired[loop_edges[i].first].sy,
						poses_desired[loop_edges[i].first].sz,
						poses_desired[loop_edges[i].second].px,
						poses_desired[loop_edges[i].second].py,
						poses_desired[loop_edges[i].second].pz,
						poses_desired[loop_edges[i].second].sx,
						poses_desired[loop_edges[i].second].sy,
						poses_desired[loop_edges[i].second].sz);

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_rodrigues_wc(
						delta,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].sx,
						poses[loop_edges[i].first].sy,
						poses[loop_edges[i].first].sz,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].sx,
						poses[loop_edges[i].second].sy,
						poses[loop_edges[i].second].sz,
						relative_pose_measurement_loop(0,0),
						relative_pose_measurement_loop(1,0),
						relative_pose_measurement_loop(2,0),
						relative_pose_measurement_loop(3,0),
						relative_pose_measurement_loop(4,0),
						relative_pose_measurement_loop(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_rodrigues_wc_jacobian(jacobian,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].sx,
						poses[loop_edges[i].first].sy,
						poses[loop_edges[i].first].sz,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].sx,
						poses[loop_edges[i].second].sy,
						poses[loop_edges[i].second].sz);

				int ir = tripletListB.size();

				int ic_1 = loop_edges[i].first * 6;
				int ic_2 = loop_edges[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
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


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 6 , m_poses.size() * 6);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 6 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 6 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					RodriguesPose pose = pose_rodrigues_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.sx += h_x[counter++];
					pose.sy += h_x[counter++];
					pose.sz += h_x[counter++];

					m_poses[i] = affine_matrix_from_pose_rodrigues(pose);

					Eigen::Vector3d vx(m_poses[i](0,0), m_poses[i](1,0), m_poses[i](2,0));
					Eigen::Vector3d vy(m_poses[i](0,1), m_poses[i](1,1), m_poses[i](2,1));
					Eigen::Vector3d vz(m_poses[i](0,2), m_poses[i](1,2), m_poses[i](2,2));

					std::cout << std::setprecision(15);
					std::cout << "norm: "<< vx.norm() << " " << vy.norm() << " " << vz.norm() << " " <<
							vx.dot(vy) << " " << vy.dot(vz) << " " << vx.dot(vz) << std::endl;
				}
				std::cout << "optimizing with rodrigues finished" << std::endl;
			}else{
				std::cout << "optimizing with rodrigues FAILED" << std::endl;
			}

			break;
		}
		case 'q':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<QuaternionPose> poses;
			std::vector<QuaternionPose> poses_desired;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_quaternion_from_affine_matrix(m_poses[i]));
			}
			for(size_t i = 0 ; i < m_poses_desired.size(); i++){
				poses_desired.push_back(pose_quaternion_from_affine_matrix(m_poses_desired[i]));
			}

			for(size_t i = 0 ; i < odo_edges.size(); i++){
				Eigen::Matrix<double, 7, 1> relative_pose_measurement_odo;
				relative_pose_quaternion_wc(relative_pose_measurement_odo,
						poses_desired[odo_edges[i].first].px,
						poses_desired[odo_edges[i].first].py,
						poses_desired[odo_edges[i].first].pz,
						poses_desired[odo_edges[i].first].q0,
						poses_desired[odo_edges[i].first].q1,
						poses_desired[odo_edges[i].first].q2,
						poses_desired[odo_edges[i].first].q3,
						poses_desired[odo_edges[i].second].px,
						poses_desired[odo_edges[i].second].py,
						poses_desired[odo_edges[i].second].pz,
						poses_desired[odo_edges[i].second].q0,
						poses_desired[odo_edges[i].second].q1,
						poses_desired[odo_edges[i].second].q2,
						poses_desired[odo_edges[i].second].q3);

				Eigen::Matrix<double, 7, 1> delta;
				relative_pose_obs_eq_quaternion_wc(
						delta,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].q0,
						poses[odo_edges[i].first].q1,
						poses[odo_edges[i].first].q2,
						poses[odo_edges[i].first].q3,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].q0,
						poses[odo_edges[i].second].q1,
						poses[odo_edges[i].second].q2,
						poses[odo_edges[i].second].q3,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0),
						relative_pose_measurement_odo(6,0));

				Eigen::Matrix<double, 7, 14, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_quaternion_wc_jacobian(jacobian,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].q0,
						poses[odo_edges[i].first].q1,
						poses[odo_edges[i].first].q2,
						poses[odo_edges[i].first].q3,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].q0,
						poses[odo_edges[i].second].q1,
						poses[odo_edges[i].second].q2,
						poses[odo_edges[i].second].q3);

				int ir = tripletListB.size();

				int ic_1 = odo_edges[i].first * 7;
				int ic_2 = odo_edges[i].second * 7;

				for(size_t row = 0 ; row < 7; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));
					tripletListA.emplace_back(ir + row, ic_1 + 6, -jacobian(row,6));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,11));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_2 + 6, -jacobian(row,13));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));
				tripletListB.emplace_back(ir + 6, 0, delta(6,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
				tripletListP.emplace_back(ir + 6, ir + 6, 1);
			}

			for(size_t i = 0 ; i < loop_edges.size(); i++){
				Eigen::Matrix<double, 7, 1> relative_pose_measurement_loop;
				relative_pose_quaternion_wc(relative_pose_measurement_loop,
						poses_desired[loop_edges[i].first].px,
						poses_desired[loop_edges[i].first].py,
						poses_desired[loop_edges[i].first].pz,
						poses_desired[loop_edges[i].first].q0,
						poses_desired[loop_edges[i].first].q1,
						poses_desired[loop_edges[i].first].q2,
						poses_desired[loop_edges[i].first].q3,
						poses_desired[loop_edges[i].second].px,
						poses_desired[loop_edges[i].second].py,
						poses_desired[loop_edges[i].second].pz,
						poses_desired[loop_edges[i].second].q0,
						poses_desired[loop_edges[i].second].q1,
						poses_desired[loop_edges[i].second].q2,
						poses_desired[loop_edges[i].second].q3);

				Eigen::Matrix<double, 7, 1> delta;
				relative_pose_obs_eq_quaternion_wc(
						delta,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].q0,
						poses[loop_edges[i].first].q1,
						poses[loop_edges[i].first].q2,
						poses[loop_edges[i].first].q3,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].q0,
						poses[loop_edges[i].second].q1,
						poses[loop_edges[i].second].q2,
						poses[loop_edges[i].second].q3,
						relative_pose_measurement_loop(0,0),
						relative_pose_measurement_loop(1,0),
						relative_pose_measurement_loop(2,0),
						relative_pose_measurement_loop(3,0),
						relative_pose_measurement_loop(4,0),
						relative_pose_measurement_loop(5,0),
						relative_pose_measurement_loop(6,0));

				Eigen::Matrix<double, 7, 14, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_quaternion_wc_jacobian(jacobian,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].q0,
						poses[loop_edges[i].first].q1,
						poses[loop_edges[i].first].q2,
						poses[loop_edges[i].first].q3,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].q0,
						poses[loop_edges[i].second].q1,
						poses[loop_edges[i].second].q2,
						poses[loop_edges[i].second].q3);

				int ir = tripletListB.size();

				int ic_1 = loop_edges[i].first * 7;
				int ic_2 = loop_edges[i].second * 7;

				for(size_t row = 0 ; row < 7; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));
					tripletListA.emplace_back(ir + row, ic_1 + 6, -jacobian(row,6));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,11));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_2 + 6, -jacobian(row,13));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));
				tripletListB.emplace_back(ir + 6, 0, delta(6,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
				tripletListP.emplace_back(ir + 6, ir + 6, 1);
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

			for(size_t i = 0 ; i < m_poses.size(); i++){
				int ic = i * 7;
				ir = tripletListB.size();
				QuaternionPose pose = pose_quaternion_from_affine_matrix(m_poses[i]);

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


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 7);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 7 , m_poses.size() * 7);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 7 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 7 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					QuaternionPose pose = pose_quaternion_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.q0 += h_x[counter++];
					pose.q1 += h_x[counter++];
					pose.q2 += h_x[counter++];
					pose.q3 += h_x[counter++];
					m_poses[i] = affine_matrix_from_pose_quaternion(pose);
				}
				std::cout << "optimizing with quaternions finished" << std::endl;
			}else{
				std::cout << "optimizing with quaternions FAILED" << std::endl;
			}

			break;
		}
		case 'x':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < odo_edges.size(); i++){
				Eigen::Matrix<double, 12, 1> relative_pose_measurement_odo;
				relative_pose_wc(relative_pose_measurement_odo,
						m_poses_desired[odo_edges[i].first](0,3),
						m_poses_desired[odo_edges[i].first](1,3),
						m_poses_desired[odo_edges[i].first](2,3),
						m_poses_desired[odo_edges[i].first](0,0),
						m_poses_desired[odo_edges[i].first](0,1),
						m_poses_desired[odo_edges[i].first](0,2),
						m_poses_desired[odo_edges[i].first](1,0),
						m_poses_desired[odo_edges[i].first](1,1),
						m_poses_desired[odo_edges[i].first](1,2),
						m_poses_desired[odo_edges[i].first](2,0),
						m_poses_desired[odo_edges[i].first](2,1),
						m_poses_desired[odo_edges[i].first](2,2),
						m_poses_desired[odo_edges[i].second](0,3),
						m_poses_desired[odo_edges[i].second](1,3),
						m_poses_desired[odo_edges[i].second](2,3),
						m_poses_desired[odo_edges[i].second](0,0),
						m_poses_desired[odo_edges[i].second](0,1),
						m_poses_desired[odo_edges[i].second](0,2),
						m_poses_desired[odo_edges[i].second](1,0),
						m_poses_desired[odo_edges[i].second](1,1),
						m_poses_desired[odo_edges[i].second](1,2),
						m_poses_desired[odo_edges[i].second](2,0),
						m_poses_desired[odo_edges[i].second](2,1),
						m_poses_desired[odo_edges[i].second](2,2));

				Eigen::Matrix<double, 12, 1> delta;
				relative_pose_obs_eq_wc(
						delta,
						m_poses[odo_edges[i].first](0,3),
						m_poses[odo_edges[i].first](1,3),
						m_poses[odo_edges[i].first](2,3),
						m_poses[odo_edges[i].first](0,0),
						m_poses[odo_edges[i].first](0,1),
						m_poses[odo_edges[i].first](0,2),
						m_poses[odo_edges[i].first](1,0),
						m_poses[odo_edges[i].first](1,1),
						m_poses[odo_edges[i].first](1,2),
						m_poses[odo_edges[i].first](2,0),
						m_poses[odo_edges[i].first](2,1),
						m_poses[odo_edges[i].first](2,2),
						m_poses[odo_edges[i].second](0,3),
						m_poses[odo_edges[i].second](1,3),
						m_poses[odo_edges[i].second](2,3),
						m_poses[odo_edges[i].second](0,0),
						m_poses[odo_edges[i].second](0,1),
						m_poses[odo_edges[i].second](0,2),
						m_poses[odo_edges[i].second](1,0),
						m_poses[odo_edges[i].second](1,1),
						m_poses[odo_edges[i].second](1,2),
						m_poses[odo_edges[i].second](2,0),
						m_poses[odo_edges[i].second](2,1),
						m_poses[odo_edges[i].second](2,2),
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0),
						relative_pose_measurement_odo(6,0),
						relative_pose_measurement_odo(7,0),
						relative_pose_measurement_odo(8,0),
						relative_pose_measurement_odo(9,0),
						relative_pose_measurement_odo(10,0),
						relative_pose_measurement_odo(11,0));

				Eigen::Matrix<double, 12, 24, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_wc_jacobian(jacobian,
						m_poses[odo_edges[i].first](0,3),
						m_poses[odo_edges[i].first](1,3),
						m_poses[odo_edges[i].first](2,3),
						m_poses[odo_edges[i].first](0,0),
						m_poses[odo_edges[i].first](0,1),
						m_poses[odo_edges[i].first](0,2),
						m_poses[odo_edges[i].first](1,0),
						m_poses[odo_edges[i].first](1,1),
						m_poses[odo_edges[i].first](1,2),
						m_poses[odo_edges[i].first](2,0),
						m_poses[odo_edges[i].first](2,1),
						m_poses[odo_edges[i].first](2,2),
						m_poses[odo_edges[i].second](0,3),
						m_poses[odo_edges[i].second](1,3),
						m_poses[odo_edges[i].second](2,3),
						m_poses[odo_edges[i].second](0,0),
						m_poses[odo_edges[i].second](0,1),
						m_poses[odo_edges[i].second](0,2),
						m_poses[odo_edges[i].second](1,0),
						m_poses[odo_edges[i].second](1,1),
						m_poses[odo_edges[i].second](1,2),
						m_poses[odo_edges[i].second](2,0),
						m_poses[odo_edges[i].second](2,1),
						m_poses[odo_edges[i].second](2,2));

				int ir = tripletListB.size();

				int ic_1 = odo_edges[i].first * 12;
				int ic_2 = odo_edges[i].second * 12;

				for(size_t row = 0 ; row < 12; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));
					tripletListA.emplace_back(ir + row, ic_1 + 6, -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_1 + 7, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_1 + 8, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_1 + 9, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_1 + 10, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_1 + 11, -jacobian(row,11));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,13));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,14));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,15));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,16));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,17));
					tripletListA.emplace_back(ir + row, ic_2 + 6, -jacobian(row,18));
					tripletListA.emplace_back(ir + row, ic_2 + 7, -jacobian(row,19));
					tripletListA.emplace_back(ir + row, ic_2 + 8, -jacobian(row,20));
					tripletListA.emplace_back(ir + row, ic_2 + 9, -jacobian(row,21));
					tripletListA.emplace_back(ir + row, ic_2 + 10, -jacobian(row,22));
					tripletListA.emplace_back(ir + row, ic_2 + 11, -jacobian(row,23));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));
				tripletListB.emplace_back(ir + 6, 0, delta(6,0));
				tripletListB.emplace_back(ir + 7, 0, delta(7,0));
				tripletListB.emplace_back(ir + 8, 0, delta(8,0));
				tripletListB.emplace_back(ir + 9, 0, delta(9,0));
				tripletListB.emplace_back(ir + 10, 0, delta(10,0));
				tripletListB.emplace_back(ir + 11, 0, delta(11,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
				tripletListP.emplace_back(ir + 6, ir + 6, 1);
				tripletListP.emplace_back(ir + 7, ir + 7, 1);
				tripletListP.emplace_back(ir + 8, ir + 8, 1);
				tripletListP.emplace_back(ir + 9, ir + 9, 1);
				tripletListP.emplace_back(ir + 10, ir + 10, 1);
				tripletListP.emplace_back(ir + 11, ir + 11, 1);
			}

			for(size_t i = 0 ; i < loop_edges.size(); i++){
				Eigen::Matrix<double, 12, 1> relative_pose_measurement_loop;
				relative_pose_wc(relative_pose_measurement_loop,
						m_poses_desired[loop_edges[i].first](0,3),
						m_poses_desired[loop_edges[i].first](1,3),
						m_poses_desired[loop_edges[i].first](2,3),
						m_poses_desired[loop_edges[i].first](0,0),
						m_poses_desired[loop_edges[i].first](0,1),
						m_poses_desired[loop_edges[i].first](0,2),
						m_poses_desired[loop_edges[i].first](1,0),
						m_poses_desired[loop_edges[i].first](1,1),
						m_poses_desired[loop_edges[i].first](1,2),
						m_poses_desired[loop_edges[i].first](2,0),
						m_poses_desired[loop_edges[i].first](2,1),
						m_poses_desired[loop_edges[i].first](2,2),
						m_poses_desired[loop_edges[i].second](0,3),
						m_poses_desired[loop_edges[i].second](1,3),
						m_poses_desired[loop_edges[i].second](2,3),
						m_poses_desired[loop_edges[i].second](0,0),
						m_poses_desired[loop_edges[i].second](0,1),
						m_poses_desired[loop_edges[i].second](0,2),
						m_poses_desired[loop_edges[i].second](1,0),
						m_poses_desired[loop_edges[i].second](1,1),
						m_poses_desired[loop_edges[i].second](1,2),
						m_poses_desired[loop_edges[i].second](2,0),
						m_poses_desired[loop_edges[i].second](2,1),
						m_poses_desired[loop_edges[i].second](2,2));

				Eigen::Matrix<double, 12, 1> delta;
				relative_pose_obs_eq_wc(
						delta,
						m_poses[loop_edges[i].first](0,3),
						m_poses[loop_edges[i].first](1,3),
						m_poses[loop_edges[i].first](2,3),
						m_poses[loop_edges[i].first](0,0),
						m_poses[loop_edges[i].first](0,1),
						m_poses[loop_edges[i].first](0,2),
						m_poses[loop_edges[i].first](1,0),
						m_poses[loop_edges[i].first](1,1),
						m_poses[loop_edges[i].first](1,2),
						m_poses[loop_edges[i].first](2,0),
						m_poses[loop_edges[i].first](2,1),
						m_poses[loop_edges[i].first](2,2),
						m_poses[loop_edges[i].second](0,3),
						m_poses[loop_edges[i].second](1,3),
						m_poses[loop_edges[i].second](2,3),
						m_poses[loop_edges[i].second](0,0),
						m_poses[loop_edges[i].second](0,1),
						m_poses[loop_edges[i].second](0,2),
						m_poses[loop_edges[i].second](1,0),
						m_poses[loop_edges[i].second](1,1),
						m_poses[loop_edges[i].second](1,2),
						m_poses[loop_edges[i].second](2,0),
						m_poses[loop_edges[i].second](2,1),
						m_poses[loop_edges[i].second](2,2),
						relative_pose_measurement_loop(0,0),
						relative_pose_measurement_loop(1,0),
						relative_pose_measurement_loop(2,0),
						relative_pose_measurement_loop(3,0),
						relative_pose_measurement_loop(4,0),
						relative_pose_measurement_loop(5,0),
						relative_pose_measurement_loop(6,0),
						relative_pose_measurement_loop(7,0),
						relative_pose_measurement_loop(8,0),
						relative_pose_measurement_loop(9,0),
						relative_pose_measurement_loop(10,0),
						relative_pose_measurement_loop(11,0));

				Eigen::Matrix<double, 12, 24, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_wc_jacobian(jacobian,
						m_poses[loop_edges[i].first](0,3),
						m_poses[loop_edges[i].first](1,3),
						m_poses[loop_edges[i].first](2,3),
						m_poses[loop_edges[i].first](0,0),
						m_poses[loop_edges[i].first](0,1),
						m_poses[loop_edges[i].first](0,2),
						m_poses[loop_edges[i].first](1,0),
						m_poses[loop_edges[i].first](1,1),
						m_poses[loop_edges[i].first](1,2),
						m_poses[loop_edges[i].first](2,0),
						m_poses[loop_edges[i].first](2,1),
						m_poses[loop_edges[i].first](2,2),
						m_poses[loop_edges[i].second](0,3),
						m_poses[loop_edges[i].second](1,3),
						m_poses[loop_edges[i].second](2,3),
						m_poses[loop_edges[i].second](0,0),
						m_poses[loop_edges[i].second](0,1),
						m_poses[loop_edges[i].second](0,2),
						m_poses[loop_edges[i].second](1,0),
						m_poses[loop_edges[i].second](1,1),
						m_poses[loop_edges[i].second](1,2),
						m_poses[loop_edges[i].second](2,0),
						m_poses[loop_edges[i].second](2,1),
						m_poses[loop_edges[i].second](2,2));

				int ir = tripletListB.size();

				int ic_1 = loop_edges[i].first * 12;
				int ic_2 = loop_edges[i].second * 12;

				for(size_t row = 0 ; row < 12; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));
					tripletListA.emplace_back(ir + row, ic_1 + 6, -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_1 + 7, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_1 + 8, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_1 + 9, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_1 + 10, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_1 + 11, -jacobian(row,11));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,13));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,14));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,15));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,16));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,17));
					tripletListA.emplace_back(ir + row, ic_2 + 6, -jacobian(row,18));
					tripletListA.emplace_back(ir + row, ic_2 + 7, -jacobian(row,19));
					tripletListA.emplace_back(ir + row, ic_2 + 8, -jacobian(row,20));
					tripletListA.emplace_back(ir + row, ic_2 + 9, -jacobian(row,21));
					tripletListA.emplace_back(ir + row, ic_2 + 10, -jacobian(row,22));
					tripletListA.emplace_back(ir + row, ic_2 + 11, -jacobian(row,23));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));
				tripletListB.emplace_back(ir + 6, 0, delta(6,0));
				tripletListB.emplace_back(ir + 7, 0, delta(7,0));
				tripletListB.emplace_back(ir + 8, 0, delta(8,0));
				tripletListB.emplace_back(ir + 9, 0, delta(9,0));
				tripletListB.emplace_back(ir + 10, 0, delta(10,0));
				tripletListB.emplace_back(ir + 11, 0, delta(11,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
				tripletListP.emplace_back(ir + 6, ir + 6, 1);
				tripletListP.emplace_back(ir + 7, ir + 7, 1);
				tripletListP.emplace_back(ir + 8, ir + 8, 1);
				tripletListP.emplace_back(ir + 9, ir + 9, 1);
				tripletListP.emplace_back(ir + 10, ir + 10, 1);
				tripletListP.emplace_back(ir + 11, ir + 11, 1);
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);
			tripletListA.emplace_back(ir + 3 , 3, 1);
			tripletListA.emplace_back(ir + 4 , 4, 1);
			tripletListA.emplace_back(ir + 5 , 5, 1);
			tripletListA.emplace_back(ir + 6 , 6, 1);
			tripletListA.emplace_back(ir + 7 , 7, 1);
			tripletListA.emplace_back(ir + 8 , 8, 1);
			tripletListA.emplace_back(ir + 9 , 9, 1);
			tripletListA.emplace_back(ir + 10 , 10, 1);
			tripletListA.emplace_back(ir + 11 , 11, 1);

			tripletListP.emplace_back(ir     , ir,     10000000000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 10000000000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 10000000000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 10000000000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 10000000000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 10000000000000);
			tripletListP.emplace_back(ir + 6 , ir + 6, 10000000000000);
			tripletListP.emplace_back(ir + 7 , ir + 7, 10000000000000);
			tripletListP.emplace_back(ir + 8 , ir + 8, 10000000000000);
			tripletListP.emplace_back(ir + 9 , ir + 9, 10000000000000);
			tripletListP.emplace_back(ir + 10 , ir + 10, 10000000000000);
			tripletListP.emplace_back(ir + 11 , ir + 11, 10000000000000);


			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);
			tripletListB.emplace_back(ir + 6 , 0, 0);
			tripletListB.emplace_back(ir + 7 , 0, 0);
			tripletListB.emplace_back(ir + 8 , 0, 0);
			tripletListB.emplace_back(ir + 9 , 0, 0);
			tripletListB.emplace_back(ir + 10 , 0, 0);
			tripletListB.emplace_back(ir + 11 , 0, 0);


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 12);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 12 , m_poses.size() * 12);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 12 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;
			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 12 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					m_poses[i](0,3) += h_x[counter++];
					m_poses[i](1,3) += h_x[counter++];
					m_poses[i](2,3) += h_x[counter++];
					m_poses[i](0,0) += h_x[counter++];
					m_poses[i](0,1) += h_x[counter++];
					m_poses[i](0,2) += h_x[counter++];
					m_poses[i](1,0) += h_x[counter++];
					m_poses[i](1,1) += h_x[counter++];
					m_poses[i](1,2) += h_x[counter++];
					m_poses[i](2,0) += h_x[counter++];
					m_poses[i](2,1) += h_x[counter++];
					m_poses[i](2,2) += h_x[counter++];

					orthogonalize_rotation(m_poses[i]);
				}
				std::cout << "optimizing without rotation matrix parametrization finished" << std::endl;
			}else{
				std::cout << "optimizing without rotation matrix parametrization FAILED" << std::endl;
			}
			break;
		}
		case 'a':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<TaitBryanPose> poses;
			std::vector<TaitBryanPose> poses_desired;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_tait_bryan_from_affine_matrix(m_poses[i]));
			}
			for(size_t i = 0 ; i < m_poses_desired.size(); i++){
				poses_desired.push_back(pose_tait_bryan_from_affine_matrix(m_poses_desired[i]));
			}

			for(size_t i = 0 ; i < odo_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_odo;
				relative_pose_tait_bryan_wc_case1(relative_pose_measurement_odo,
						poses_desired[odo_edges[i].first].px,
						poses_desired[odo_edges[i].first].py,
						poses_desired[odo_edges[i].first].pz,
						poses_desired[odo_edges[i].first].om,
						poses_desired[odo_edges[i].first].fi,
						poses_desired[odo_edges[i].first].ka,
						poses_desired[odo_edges[i].second].px,
						poses_desired[odo_edges[i].second].py,
						poses_desired[odo_edges[i].second].pz,
						poses_desired[odo_edges[i].second].om,
						poses_desired[odo_edges[i].second].fi,
						poses_desired[odo_edges[i].second].ka);

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_2_obs_eq_tait_bryan_wc_case1(
						delta,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].om,
						poses[odo_edges[i].first].fi,
						poses[odo_edges[i].first].ka,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].om,
						poses[odo_edges[i].second].fi,
						poses[odo_edges[i].second].ka,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_2_obs_eq_tait_bryan_wc_case1_jacobian(jacobian,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].om,
						poses[odo_edges[i].first].fi,
						poses[odo_edges[i].first].ka,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].om,
						poses[odo_edges[i].second].fi,
						poses[odo_edges[i].second].ka,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0));

				int ir = tripletListB.size();

				int ic_1 = odo_edges[i].first * 6;
				int ic_2 = odo_edges[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
			}

			for(size_t i = 0 ; i < loop_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_loop;
				relative_pose_tait_bryan_wc_case1(relative_pose_measurement_loop,
						poses_desired[loop_edges[i].first].px,
						poses_desired[loop_edges[i].first].py,
						poses_desired[loop_edges[i].first].pz,
						poses_desired[loop_edges[i].first].om,
						poses_desired[loop_edges[i].first].fi,
						poses_desired[loop_edges[i].first].ka,
						poses_desired[loop_edges[i].second].px,
						poses_desired[loop_edges[i].second].py,
						poses_desired[loop_edges[i].second].pz,
						poses_desired[loop_edges[i].second].om,
						poses_desired[loop_edges[i].second].fi,
						poses_desired[loop_edges[i].second].ka);

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_tait_bryan_wc_case1(
						delta,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].om,
						poses[loop_edges[i].first].fi,
						poses[loop_edges[i].first].ka,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].om,
						poses[loop_edges[i].second].fi,
						poses[loop_edges[i].second].ka,
						relative_pose_measurement_loop(0,0),
						relative_pose_measurement_loop(1,0),
						relative_pose_measurement_loop(2,0),
						relative_pose_measurement_loop(3,0),
						relative_pose_measurement_loop(4,0),
						relative_pose_measurement_loop(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_tait_bryan_wc_case1_jacobian(jacobian,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].om,
						poses[loop_edges[i].first].fi,
						poses[loop_edges[i].first].ka,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].om,
						poses[loop_edges[i].second].fi,
						poses[loop_edges[i].second].ka);

				int ir = tripletListB.size();

				int ic_1 = loop_edges[i].first * 6;
				int ic_2 = loop_edges[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
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


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 6 , m_poses.size() * 6);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 6 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 6 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];
					m_poses[i] = affine_matrix_from_pose_tait_bryan(pose);
				}
				std::cout << "optimizing with tait bryan finished" << std::endl;
			}else{
				std::cout << "optimizing with tait bryan FAILED" << std::endl;
			}

			break;
		}
		case 's':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<RodriguesPose> poses;
			std::vector<RodriguesPose> poses_desired;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_rodrigues_from_affine_matrix(m_poses[i]));
			}
			for(size_t i = 0 ; i < m_poses_desired.size(); i++){
				poses_desired.push_back(pose_rodrigues_from_affine_matrix(m_poses_desired[i]));
			}

			for(size_t i = 0 ; i < odo_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_odo;
				relative_pose_rodrigues_wc(relative_pose_measurement_odo,
						poses_desired[odo_edges[i].first].px,
						poses_desired[odo_edges[i].first].py,
						poses_desired[odo_edges[i].first].pz,
						poses_desired[odo_edges[i].first].sx,
						poses_desired[odo_edges[i].first].sy,
						poses_desired[odo_edges[i].first].sz,
						poses_desired[odo_edges[i].second].px,
						poses_desired[odo_edges[i].second].py,
						poses_desired[odo_edges[i].second].pz,
						poses_desired[odo_edges[i].second].sx,
						poses_desired[odo_edges[i].second].sy,
						poses_desired[odo_edges[i].second].sz);

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_2_obs_eq_rodrigues_wc(
						delta,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].sx,
						poses[odo_edges[i].first].sy,
						poses[odo_edges[i].first].sz,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].sx,
						poses[odo_edges[i].second].sy,
						poses[odo_edges[i].second].sz,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_2_obs_eq_rodrigues_wc_jacobian(jacobian,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].sx,
						poses[odo_edges[i].first].sy,
						poses[odo_edges[i].first].sz,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].sx,
						poses[odo_edges[i].second].sy,
						poses[odo_edges[i].second].sz,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0));

				int ir = tripletListB.size();

				int ic_1 = odo_edges[i].first * 6;
				int ic_2 = odo_edges[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
			}

			for(size_t i = 0 ; i < loop_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_loop;
				relative_pose_rodrigues_wc(relative_pose_measurement_loop,
						poses_desired[loop_edges[i].first].px,
						poses_desired[loop_edges[i].first].py,
						poses_desired[loop_edges[i].first].pz,
						poses_desired[loop_edges[i].first].sx,
						poses_desired[loop_edges[i].first].sy,
						poses_desired[loop_edges[i].first].sz,
						poses_desired[loop_edges[i].second].px,
						poses_desired[loop_edges[i].second].py,
						poses_desired[loop_edges[i].second].pz,
						poses_desired[loop_edges[i].second].sx,
						poses_desired[loop_edges[i].second].sy,
						poses_desired[loop_edges[i].second].sz);

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_rodrigues_wc(
						delta,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].sx,
						poses[loop_edges[i].first].sy,
						poses[loop_edges[i].first].sz,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].sx,
						poses[loop_edges[i].second].sy,
						poses[loop_edges[i].second].sz,
						relative_pose_measurement_loop(0,0),
						relative_pose_measurement_loop(1,0),
						relative_pose_measurement_loop(2,0),
						relative_pose_measurement_loop(3,0),
						relative_pose_measurement_loop(4,0),
						relative_pose_measurement_loop(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_rodrigues_wc_jacobian(jacobian,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].sx,
						poses[loop_edges[i].first].sy,
						poses[loop_edges[i].first].sz,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].sx,
						poses[loop_edges[i].second].sy,
						poses[loop_edges[i].second].sz);

				int ir = tripletListB.size();

				int ic_1 = loop_edges[i].first * 6;
				int ic_2 = loop_edges[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
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


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 6 , m_poses.size() * 6);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 6 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 6 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					RodriguesPose pose = pose_rodrigues_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.sx += h_x[counter++];
					pose.sy += h_x[counter++];
					pose.sz += h_x[counter++];

					m_poses[i] = affine_matrix_from_pose_rodrigues(pose);

					Eigen::Vector3d vx(m_poses[i](0,0), m_poses[i](1,0), m_poses[i](2,0));
					Eigen::Vector3d vy(m_poses[i](0,1), m_poses[i](1,1), m_poses[i](2,1));
					Eigen::Vector3d vz(m_poses[i](0,2), m_poses[i](1,2), m_poses[i](2,2));

					std::cout << std::setprecision(15);
					std::cout << "norm: "<< vx.norm() << " " << vy.norm() << " " << vz.norm() << " " <<
							vx.dot(vy) << " " << vy.dot(vz) << " " << vx.dot(vz) << std::endl;
				}
				std::cout << "optimizing with rodrigues finished" << std::endl;
			}else{
				std::cout << "optimizing with rodrigues FAILED" << std::endl;
			}

			break;
		}
		case 'd':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<QuaternionPose> poses;
			std::vector<QuaternionPose> poses_desired;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_quaternion_from_affine_matrix(m_poses[i]));
			}
			for(size_t i = 0 ; i < m_poses_desired.size(); i++){
				poses_desired.push_back(pose_quaternion_from_affine_matrix(m_poses_desired[i]));
			}

			for(size_t i = 0 ; i < odo_edges.size(); i++){
				Eigen::Matrix<double, 7, 1> relative_pose_measurement_odo;
				relative_pose_quaternion_wc(relative_pose_measurement_odo,
						poses_desired[odo_edges[i].first].px,
						poses_desired[odo_edges[i].first].py,
						poses_desired[odo_edges[i].first].pz,
						poses_desired[odo_edges[i].first].q0,
						poses_desired[odo_edges[i].first].q1,
						poses_desired[odo_edges[i].first].q2,
						poses_desired[odo_edges[i].first].q3,
						poses_desired[odo_edges[i].second].px,
						poses_desired[odo_edges[i].second].py,
						poses_desired[odo_edges[i].second].pz,
						poses_desired[odo_edges[i].second].q0,
						poses_desired[odo_edges[i].second].q1,
						poses_desired[odo_edges[i].second].q2,
						poses_desired[odo_edges[i].second].q3);

				Eigen::Matrix<double, 7, 1> delta;
				relative_pose_2_obs_eq_quaternion_wc(
						delta,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].q0,
						poses[odo_edges[i].first].q1,
						poses[odo_edges[i].first].q2,
						poses[odo_edges[i].first].q3,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].q0,
						poses[odo_edges[i].second].q1,
						poses[odo_edges[i].second].q2,
						poses[odo_edges[i].second].q3,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0),
						relative_pose_measurement_odo(6,0));

				Eigen::Matrix<double, 7, 14, Eigen::RowMajor> jacobian;
				relative_pose_2_obs_eq_quaternion_wc_jacobian(jacobian,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].q0,
						poses[odo_edges[i].first].q1,
						poses[odo_edges[i].first].q2,
						poses[odo_edges[i].first].q3,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].q0,
						poses[odo_edges[i].second].q1,
						poses[odo_edges[i].second].q2,
						poses[odo_edges[i].second].q3,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0),
						relative_pose_measurement_odo(6,0));

				int ir = tripletListB.size();

				int ic_1 = odo_edges[i].first * 7;
				int ic_2 = odo_edges[i].second * 7;

				for(size_t row = 0 ; row < 7; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));
					tripletListA.emplace_back(ir + row, ic_1 + 6, -jacobian(row,6));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,11));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_2 + 6, -jacobian(row,13));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));
				tripletListB.emplace_back(ir + 6, 0, delta(6,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
				tripletListP.emplace_back(ir + 6, ir + 6, 1);
			}

			for(size_t i = 0 ; i < loop_edges.size(); i++){
				Eigen::Matrix<double, 7, 1> relative_pose_measurement_loop;
				relative_pose_quaternion_wc(relative_pose_measurement_loop,
						poses_desired[loop_edges[i].first].px,
						poses_desired[loop_edges[i].first].py,
						poses_desired[loop_edges[i].first].pz,
						poses_desired[loop_edges[i].first].q0,
						poses_desired[loop_edges[i].first].q1,
						poses_desired[loop_edges[i].first].q2,
						poses_desired[loop_edges[i].first].q3,
						poses_desired[loop_edges[i].second].px,
						poses_desired[loop_edges[i].second].py,
						poses_desired[loop_edges[i].second].pz,
						poses_desired[loop_edges[i].second].q0,
						poses_desired[loop_edges[i].second].q1,
						poses_desired[loop_edges[i].second].q2,
						poses_desired[loop_edges[i].second].q3);

				Eigen::Matrix<double, 7, 1> delta;
				relative_pose_obs_eq_quaternion_wc(
						delta,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].q0,
						poses[loop_edges[i].first].q1,
						poses[loop_edges[i].first].q2,
						poses[loop_edges[i].first].q3,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].q0,
						poses[loop_edges[i].second].q1,
						poses[loop_edges[i].second].q2,
						poses[loop_edges[i].second].q3,
						relative_pose_measurement_loop(0,0),
						relative_pose_measurement_loop(1,0),
						relative_pose_measurement_loop(2,0),
						relative_pose_measurement_loop(3,0),
						relative_pose_measurement_loop(4,0),
						relative_pose_measurement_loop(5,0),
						relative_pose_measurement_loop(6,0));

				Eigen::Matrix<double, 7, 14, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_quaternion_wc_jacobian(jacobian,
						poses[loop_edges[i].first].px,
						poses[loop_edges[i].first].py,
						poses[loop_edges[i].first].pz,
						poses[loop_edges[i].first].q0,
						poses[loop_edges[i].first].q1,
						poses[loop_edges[i].first].q2,
						poses[loop_edges[i].first].q3,
						poses[loop_edges[i].second].px,
						poses[loop_edges[i].second].py,
						poses[loop_edges[i].second].pz,
						poses[loop_edges[i].second].q0,
						poses[loop_edges[i].second].q1,
						poses[loop_edges[i].second].q2,
						poses[loop_edges[i].second].q3);

				int ir = tripletListB.size();

				int ic_1 = loop_edges[i].first * 7;
				int ic_2 = loop_edges[i].second * 7;

				for(size_t row = 0 ; row < 7; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));
					tripletListA.emplace_back(ir + row, ic_1 + 6, -jacobian(row,6));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,11));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_2 + 6, -jacobian(row,13));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));
				tripletListB.emplace_back(ir + 6, 0, delta(6,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
				tripletListP.emplace_back(ir + 6, ir + 6, 1);
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

			for(size_t i = 0 ; i < m_poses.size(); i++){
				int ic = i * 7;
				ir = tripletListB.size();
				QuaternionPose pose = pose_quaternion_from_affine_matrix(m_poses[i]);

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


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 7);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 7 , m_poses.size() * 7);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 7 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 7 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					QuaternionPose pose = pose_quaternion_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.q0 += h_x[counter++];
					pose.q1 += h_x[counter++];
					pose.q2 += h_x[counter++];
					pose.q3 += h_x[counter++];
					m_poses[i] = affine_matrix_from_pose_quaternion(pose);
				}
				std::cout << "optimizing with quaternions finished" << std::endl;
			}else{
				std::cout << "optimizing with quaternions FAILED" << std::endl;
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
	std::cout << "n: add noise to poses" << std::endl;
	std::cout << "t: optimize (Tait-Bryan)" << std::endl;
	std::cout << "r: optimize (Rodriguez)" << std::endl;
	std::cout << "q: optimize (Quaternion)" << std::endl;
	std::cout << "x: optimize (Without rotation matrix parametrization)" << std::endl;
	std::cout << "a: optimize (Tait-Bryan 2)" << std::endl;
	std::cout << "s: optimize (Rodriguez 2)" << std::endl;
	std::cout << "d: optimize (Quaternion 2)" << std::endl;

}







