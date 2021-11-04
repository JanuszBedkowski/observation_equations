#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "../../smoothness_tait_bryan_wc_jacobian.h"
#include "../../smoothness_rodrigues_wc_jacobian.h"
#include "../../smoothness_quaternion_wc_jacobian.h"
#include "../../quaternion_constraint_jacobian.h"

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
	glutCreateWindow("smoothness");
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
	glBegin(GL_LINE_STRIP);
	for(size_t i = 0; i < m_poses.size(); i++){
		glVertex3f(m_poses[i](0,3), m_poses[i](1,3), m_poses[i](2,3) );
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
			for(size_t i = 1 ; i < m_poses.size(); i++){
				TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m_poses[i]);
				pose.px += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.py += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.pz += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.om += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.001;
				pose.fi += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.001;
				pose.ka += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.001;
				m_poses[i] = affine_matrix_from_pose_tait_bryan(pose);
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<TaitBryanPose> poses;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_tait_bryan_from_affine_matrix(m_poses[i]));
			}

			for(size_t i = 1; i < m_poses.size() - 1; i++){
				Eigen::Matrix<double, 6, 1> delta;
				smoothness_obs_eq_tait_bryan_wc(delta,
						poses[i-1].px,
						poses[i-1].py,
						poses[i-1].pz,
						poses[i-1].om,
						poses[i-1].fi,
						poses[i-1].ka,
						poses[i].px,
						poses[i].py,
						poses[i].pz,
						poses[i].om,
						poses[i].fi,
						poses[i].ka,
						poses[i+1].px,
						poses[i+1].py,
						poses[i+1].pz,
						poses[i+1].om,
						poses[i+1].fi,
						poses[i+1].ka);

				Eigen::Matrix<double, 6, 18, Eigen::RowMajor> jacobian;
				smoothness_obs_eq_tait_bryan_wc_jacobian(jacobian,
						poses[i-1].px,
						poses[i-1].py,
						poses[i-1].pz,
						poses[i-1].om,
						poses[i-1].fi,
						poses[i-1].ka,
						poses[i].px,
						poses[i].py,
						poses[i].pz,
						poses[i].om,
						poses[i].fi,
						poses[i].ka,
						poses[i+1].px,
						poses[i+1].py,
						poses[i+1].pz,
						poses[i+1].om,
						poses[i+1].fi,
						poses[i+1].ka);

				int ir = tripletListB.size();

				int ic_1 = (i-1) * 6;
				int ic_2 = i * 6;
				int ic_3 = (i+1) * 6;

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

					tripletListA.emplace_back(ir + row, ic_3    , -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_3 + 1, -jacobian(row,13));
					tripletListA.emplace_back(ir + row, ic_3 + 2, -jacobian(row,14));
					tripletListA.emplace_back(ir + row, ic_3 + 3, -jacobian(row,15));
					tripletListA.emplace_back(ir + row, ic_3 + 4, -jacobian(row,16));
					tripletListA.emplace_back(ir + row, ic_3 + 5, -jacobian(row,17));
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

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_rodrigues_from_affine_matrix(m_poses[i]));
			}

			for(size_t i = 1; i < m_poses.size() - 1; i++){
				Eigen::Matrix<double, 6, 1> delta;
				smoothness_obs_eq_rodrigues_wc(delta,
						poses[i-1].px,
						poses[i-1].py,
						poses[i-1].pz,
						poses[i-1].sx,
						poses[i-1].sy,
						poses[i-1].sz,
						poses[i].px,
						poses[i].py,
						poses[i].pz,
						poses[i].sx,
						poses[i].sy,
						poses[i].sz,
						poses[i+1].px,
						poses[i+1].py,
						poses[i+1].pz,
						poses[i+1].sx,
						poses[i+1].sy,
						poses[i+1].sz);

				Eigen::Matrix<double, 6, 18, Eigen::RowMajor> jacobian;
				smoothness_obs_eq_rodrigues_wc_jacobian(jacobian,
						poses[i-1].px,
						poses[i-1].py,
						poses[i-1].pz,
						poses[i-1].sx,
						poses[i-1].sy,
						poses[i-1].sz,
						poses[i].px,
						poses[i].py,
						poses[i].pz,
						poses[i].sx,
						poses[i].sy,
						poses[i].sz,
						poses[i+1].px,
						poses[i+1].py,
						poses[i+1].pz,
						poses[i+1].sx,
						poses[i+1].sy,
						poses[i+1].sz);

				int ir = tripletListB.size();

				int ic_1 = (i-1) * 6;
				int ic_2 = i * 6;
				int ic_3 = (i+1) * 6;

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

					tripletListA.emplace_back(ir + row, ic_3    , -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_3 + 1, -jacobian(row,13));
					tripletListA.emplace_back(ir + row, ic_3 + 2, -jacobian(row,14));
					tripletListA.emplace_back(ir + row, ic_3 + 3, -jacobian(row,15));
					tripletListA.emplace_back(ir + row, ic_3 + 4, -jacobian(row,16));
					tripletListA.emplace_back(ir + row, ic_3 + 5, -jacobian(row,17));
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

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_quaternion_from_affine_matrix(m_poses[i]));
			}

			for(size_t i = 1; i < m_poses.size() - 1; i++){
				Eigen::Matrix<double, 7, 1> delta;
				smoothness_obs_eq_quaternion_wc(delta,
						poses[i-1].px,
						poses[i-1].py,
						poses[i-1].pz,
						poses[i-1].q0,
						poses[i-1].q1,
						poses[i-1].q2,
						poses[i-1].q3,
						poses[i].px,
						poses[i].py,
						poses[i].pz,
						poses[i].q0,
						poses[i].q1,
						poses[i].q2,
						poses[i].q3,
						poses[i+1].px,
						poses[i+1].py,
						poses[i+1].pz,
						poses[i+1].q0,
						poses[i+1].q1,
						poses[i+1].q2,
						poses[i+1].q3);

				Eigen::Matrix<double, 7, 21, Eigen::RowMajor> jacobian;
				smoothness_obs_eq_quaternion_wc_jacobian(jacobian,
						poses[i-1].px,
						poses[i-1].py,
						poses[i-1].pz,
						poses[i-1].q0,
						poses[i-1].q1,
						poses[i-1].q2,
						poses[i-1].q3,
						poses[i].px,
						poses[i].py,
						poses[i].pz,
						poses[i].q0,
						poses[i].q1,
						poses[i].q2,
						poses[i].q3,
						poses[i+1].px,
						poses[i+1].py,
						poses[i+1].pz,
						poses[i+1].q0,
						poses[i+1].q1,
						poses[i+1].q2,
						poses[i+1].q3);

				int ir = tripletListB.size();

				int ic_1 = (i-1) * 7;
				int ic_2 = i * 7;
				int ic_3 = (i+1) * 7;

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

					tripletListA.emplace_back(ir + row, ic_3    , -jacobian(row,14));
					tripletListA.emplace_back(ir + row, ic_3 + 1, -jacobian(row,15));
					tripletListA.emplace_back(ir + row, ic_3 + 2, -jacobian(row,16));
					tripletListA.emplace_back(ir + row, ic_3 + 3, -jacobian(row,17));
					tripletListA.emplace_back(ir + row, ic_3 + 4, -jacobian(row,18));
					tripletListA.emplace_back(ir + row, ic_3 + 5, -jacobian(row,19));
					tripletListA.emplace_back(ir + row, ic_3 + 6, -jacobian(row,20));
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
}







