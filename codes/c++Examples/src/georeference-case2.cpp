#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"
#include "relative_pose_tait_bryan_wc_jacobian.h"
#include "constraint_fixed_parameter_jacobian.h"

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

std::vector<std::pair<Eigen::Affine3d, int>> georeference_data;

int main(int argc, char *argv[]){

	for(size_t i = 0 ; i < 100; i++){
		TaitBryanPose p;
		p.px = i;
		p.py = -1;
		p.pz = 0.0;
		p.om = random(-0.01, 0.01);
		p.fi = random(-0.01, 0.01);
		p.ka = random(-0.01, 0.01);

		Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(p);
		m_poses.push_back(m);
	}
	for(size_t i = 0 ; i < 100; i++){
		TaitBryanPose p;
		p.px = i;
		p.py = 1;
		p.pz = 0.0;
		p.om = random(-0.01, 0.01);
		p.fi = random(-0.01, 0.01);
		p.ka = random(-0.01, 0.01);

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

	TaitBryanPose p;
	p.px = 0;
	p.py = 5;
	p.pz = 5;
	p.om = 0;
	p.fi = 0;
	p.ka = 0;
	Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(p);
	georeference_data.emplace_back(m, 4);

	p.px = 50;
	p.py = 3;
	p.pz = 4;
	m = affine_matrix_from_pose_tait_bryan(p);
	georeference_data.emplace_back(m, 54);

	p.px = 0;
	p.py = 8;
	p.pz = 5;
	m = affine_matrix_from_pose_tait_bryan(p);
	georeference_data.emplace_back(m, 102);

	p.px = 110;
	p.py = 8;
	p.pz = 7;
	m = affine_matrix_from_pose_tait_bryan(p);
	georeference_data.emplace_back(m, 197);

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
	glutCreateWindow("georeference-case2");
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

	for(size_t i = 0 ; i < georeference_data.size(); i++){
		glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(georeference_data[i].first(0,3), georeference_data[i].first(1,3), georeference_data[i].first(2,3));
		glVertex3f(georeference_data[i].first(0,3) + georeference_data[i].first(0,0), georeference_data[i].first(1,3) + georeference_data[i].first(1,0), georeference_data[i].first(2,3) + georeference_data[i].first(2,0));

		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(georeference_data[i].first(0,3), georeference_data[i].first(1,3), georeference_data[i].first(2,3));
		glVertex3f(georeference_data[i].first(0,3) + georeference_data[i].first(0,1), georeference_data[i].first(1,3) + georeference_data[i].first(1,1), georeference_data[i].first(2,3) + georeference_data[i].first(2,1));

		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(georeference_data[i].first(0,3), georeference_data[i].first(1,3), georeference_data[i].first(2,3));
		glVertex3f(georeference_data[i].first(0,3) + georeference_data[i].first(0,2), georeference_data[i].first(1,3) + georeference_data[i].first(1,2), georeference_data[i].first(2,3) + georeference_data[i].first(2,2));
		glEnd();
	}

	glColor3f(0.0f, 0.0f, 0.0f);
	glBegin(GL_LINES);
	for(size_t i = 0 ; i < georeference_data.size(); i++){
		glVertex3f(georeference_data[i].first(0,3), georeference_data[i].first(1,3), georeference_data[i].first(2,3));
		glVertex3f(m_poses[georeference_data[i].second](0,3), m_poses[georeference_data[i].second](1,3), m_poses[georeference_data[i].second](2,3) );
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
				pose.px += random(-0.1, 0.1);
				pose.py += random(-0.1, 0.1);
				pose.pz += random(-0.1, 0.1);
				pose.om += random(-0.01, 0.01);
				pose.fi += random(-0.01, 0.01);
				pose.ka += random(-0.01, 0.01);
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

				tripletListP.emplace_back(ir ,    ir,     1000);
				tripletListP.emplace_back(ir + 1, ir + 1, 1000);
				tripletListP.emplace_back(ir + 2, ir + 2, 1000);
				tripletListP.emplace_back(ir + 3, ir + 3, 1000);
				tripletListP.emplace_back(ir + 4, ir + 4, 1000);
				tripletListP.emplace_back(ir + 5, ir + 5, 1000);
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

			for(size_t i = 0; i < georeference_data.size(); i++){
				TaitBryanPose pose_gps = pose_tait_bryan_from_affine_matrix(georeference_data[i].first);
				TaitBryanPose pose_s = pose_tait_bryan_from_affine_matrix(m_poses[georeference_data[i].second]);

				Eigen::Matrix<double, 6, 1> relative_pose_measurement;
				relative_pose_measurement(0,0) = 0.0;
				relative_pose_measurement(1,0) = 0.0;
				relative_pose_measurement(2,0) = 0.0;
				relative_pose_measurement(3,0) = 0.0;
				relative_pose_measurement(4,0) = 0.0;
				relative_pose_measurement(5,0) = 0.0;


				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_tait_bryan_wc_case1(
						delta,
						pose_s.px,
						pose_s.py,
						pose_s.pz,
						pose_s.om,
						pose_s.fi,
						pose_s.ka,
						pose_gps.px,
						pose_gps.py,
						pose_gps.pz,
						pose_gps.om,
						pose_gps.fi,
						pose_gps.ka,
						relative_pose_measurement(0,0),
						relative_pose_measurement(1,0),
						relative_pose_measurement(2,0),
						relative_pose_measurement(3,0),
						relative_pose_measurement(4,0),
						relative_pose_measurement(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_tait_bryan_wc_case1_jacobian(jacobian,
						pose_s.px,
						pose_s.py,
						pose_s.pz,
						pose_s.om,
						pose_s.fi,
						pose_s.ka,
						pose_gps.px,
						pose_gps.py,
						pose_gps.pz,
						pose_gps.om,
						pose_gps.fi,
						pose_gps.ka);

				int ir = tripletListB.size();
				int ic_1 = georeference_data[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     cauchy(delta(0,0), 1));
				tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta(1,0), 1));
				tripletListP.emplace_back(ir + 2, ir + 2, cauchy(delta(2,0), 1));
				tripletListP.emplace_back(ir + 3, ir + 3, cauchy(delta(3,0), 1));
				tripletListP.emplace_back(ir + 4, ir + 4, cauchy(delta(4,0), 1));
				tripletListP.emplace_back(ir + 5, ir + 5, cauchy(delta(5,0), 1));
			}

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
			}
			std::cout << "optimizing with tait bryan finished" << std::endl;
			break;
		}
		case 'y':{
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

				tripletListP.emplace_back(ir ,    ir,     1000);
				tripletListP.emplace_back(ir + 1, ir + 1, 1000);
				tripletListP.emplace_back(ir + 2, ir + 2, 1000);
				tripletListP.emplace_back(ir + 3, ir + 3, 1000);
				tripletListP.emplace_back(ir + 4, ir + 4, 1000);
				tripletListP.emplace_back(ir + 5, ir + 5, 1000);
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

			for(size_t i = 0; i < georeference_data.size(); i++){
				TaitBryanPose pose_gps = pose_tait_bryan_from_affine_matrix(georeference_data[i].first);
				TaitBryanPose pose_s = pose_tait_bryan_from_affine_matrix(m_poses[georeference_data[i].second]);

				Eigen::Matrix<double, 6, 1> relative_pose_measurement;
				relative_pose_measurement(0,0) = 0.0;
				relative_pose_measurement(1,0) = 0.0;
				relative_pose_measurement(2,0) = 0.0;
				relative_pose_measurement(3,0) = 0.0;
				relative_pose_measurement(4,0) = 0.0;
				relative_pose_measurement(5,0) = 0.0;


				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_tait_bryan_wc_case1(
						delta,
						pose_s.px,
						pose_s.py,
						pose_s.pz,
						pose_s.om,
						pose_s.fi,
						pose_s.ka,
						pose_gps.px,
						pose_gps.py,
						pose_gps.pz,
						pose_gps.om,
						pose_gps.fi,
						pose_gps.ka,
						relative_pose_measurement(0,0),
						relative_pose_measurement(1,0),
						relative_pose_measurement(2,0),
						relative_pose_measurement(3,0),
						relative_pose_measurement(4,0),
						relative_pose_measurement(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_tait_bryan_wc_case1_jacobian(jacobian,
						pose_s.px,
						pose_s.py,
						pose_s.pz,
						pose_s.om,
						pose_s.fi,
						pose_s.ka,
						pose_gps.px,
						pose_gps.py,
						pose_gps.pz,
						pose_gps.om,
						pose_gps.fi,
						pose_gps.ka);

				int ir = tripletListB.size();
				int ic_1 = georeference_data[i].second * 6;
				int ic_2 = m_poses.size() * 6 + i * 6;

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

				tripletListP.emplace_back(ir ,    ir,     cauchy(delta(0,0), 1));
				tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta(1,0), 1));
				tripletListP.emplace_back(ir + 2, ir + 2, cauchy(delta(2,0), 1));
				tripletListP.emplace_back(ir + 3, ir + 3, cauchy(delta(3,0), 1));
				tripletListP.emplace_back(ir + 4, ir + 4, cauchy(delta(4,0), 1));
				tripletListP.emplace_back(ir + 5, ir + 5, cauchy(delta(5,0), 1));
			}

			for(size_t i = 0; i < georeference_data.size(); i++){
				int ir = tripletListB.size();
				int ic = m_poses.size() * 6 + i * 6;
				TaitBryanPose pose_gps = pose_tait_bryan_from_affine_matrix(georeference_data[i].first);
				double residual;
				residual_constraint_fixed_optimization_parameter(residual, pose_gps.px, pose_gps.px);
				Eigen::Matrix<double, 1, 1> jacobian;
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, pose_gps.px, pose_gps.px);

				tripletListA.emplace_back(ir, ic  , -jacobian(0,0));
				tripletListB.emplace_back(ir, 0, residual);
				tripletListP.emplace_back(ir, ir, 1000000000);

				residual_constraint_fixed_optimization_parameter(residual, pose_gps.py, pose_gps.py);
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, pose_gps.py, pose_gps.py);
				tripletListA.emplace_back(ir + 1, ic + 1, -jacobian(0,0));
				tripletListB.emplace_back(ir + 1, 0, residual);
				tripletListP.emplace_back(ir + 1, ir + 1, 1000000000);

				residual_constraint_fixed_optimization_parameter(residual, pose_gps.pz, pose_gps.pz);
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, pose_gps.pz, pose_gps.pz);
				tripletListA.emplace_back(ir + 2, ic + 2, -jacobian(0,0));
				tripletListB.emplace_back(ir + 2, 0, residual);
				tripletListP.emplace_back(ir + 2, ir + 2, 1000000000);

				residual_constraint_fixed_optimization_parameter(residual, pose_gps.om, pose_gps.om);
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, pose_gps.om, pose_gps.om);
				tripletListA.emplace_back(ir + 3, ic + 3, -jacobian(0,0));
				tripletListB.emplace_back(ir + 3, 0, residual);
				tripletListP.emplace_back(ir + 3, ir + 3, 1000000000);

				residual_constraint_fixed_optimization_parameter(residual, pose_gps.fi, pose_gps.fi);
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, pose_gps.fi, pose_gps.fi);
				tripletListA.emplace_back(ir + 4, ic + 4, -jacobian(0,0));
				tripletListB.emplace_back(ir + 4, 0, residual);
				tripletListP.emplace_back(ir + 4, ir + 4, 1000000000);

				residual_constraint_fixed_optimization_parameter(residual, pose_gps.ka, pose_gps.ka);
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, pose_gps.ka, pose_gps.ka);
				tripletListA.emplace_back(ir + 5, ic + 5, -jacobian(0,0));
				tripletListB.emplace_back(ir + 5, 0, residual);
				tripletListP.emplace_back(ir + 5, ir + 5, 1000000000);
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 6 + georeference_data.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 6 + georeference_data.size() * 6, m_poses.size() * 6 + georeference_data.size() * 6);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 6 + georeference_data.size() * 6 , 1);

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

			if(h_x.size() == 6 * m_poses.size() + georeference_data.size() * 6){
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
	std::cout << "t: optimize (Tait-Bryan and relative pose to GPS (half A))" << std::endl;
	std::cout << "y: optimize (Tait-Bryan and relative pose to GPS (full A))" << std::endl;
}







