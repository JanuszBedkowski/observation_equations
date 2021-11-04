#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"
#include "../../point_to_point_vo_scale_tait_bryan_wc_jacobian.h"
#include "../../point_to_point_vo_scale_rodrigues_wc_jacobian.h"
#include "../../point_to_point_vo_scale_quaternion_wc_jacobian.h"
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


std::vector<Eigen::Affine3d> trajectory;
std::vector<Eigen::Vector3d> points_global;
std::vector<std::vector<Eigen::Vector3d>> points_local;

double scale = 1.0;

int main(int argc, char *argv[]){

	for(size_t i = 0 ; i < 10000; i++){
		Eigen::Vector3d p;
		p.x() = ((float(rand()%1000000))/1000000.0f - 0.5) * 100;
		p.y() = ((float(rand()%1000000))/1000000.0f - 0.5) * 100;
		p.z() = ((float(rand()%1000000))/1000000.0f - 0.5) * 100;
		points_global.push_back(p);
	}

	TaitBryanPose pose;
	pose.px = -4;
	pose.py = 0;
	pose.pz = 0;
	pose.om = 0;
	pose.fi = 0;
	pose.ka = 0;
	trajectory.push_back(affine_matrix_from_pose_tait_bryan(pose));

	pose.px = -2;
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

	pose.px = 2;
	pose.py = 0;
	pose.pz = 0;
	pose.om = 0;
	pose.fi = 0;
	pose.ka = 0;
	trajectory.push_back(affine_matrix_from_pose_tait_bryan(pose));

	pose.px = 4;
	pose.py = 0;
	pose.pz = 0;
	pose.om = 0;
	pose.fi = 0;
	pose.ka = 0;
	trajectory.push_back(affine_matrix_from_pose_tait_bryan(pose));


	for(size_t i = 0; i < trajectory.size(); i++){
		Eigen::Affine3d m_inv = trajectory[i].inverse();
		std::vector<Eigen::Vector3d> pp;

		for(size_t j = 0 ; j < points_global.size(); j++){
			Eigen::Vector3d vt = m_inv * points_global[j];
			pp.push_back(vt);
		}
		points_local.push_back(pp);
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
	glutCreateWindow("optimizing_scale_VO_with_lidar_data");
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
		Eigen::Affine3d m = trajectory[i];

		m(0,3) *= scale;
		m(1,3) *= scale;
		m(2,3) *= scale;

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

	glColor3f(1,0,0);
	glPointSize(3);
	glBegin(GL_POINTS);

	for(size_t i = 0 ; i < trajectory.size(); i++){
		Eigen::Affine3d m = trajectory[i];

		m(0,3) *= scale;
		m(1,3) *= scale;
		m(2,3) *= scale;

		if(i == 0)glColor3f(1,0,0);
		if(i == 1)glColor3f(0,1,0);
		if(i == 2)glColor3f(0,0,1);
		if(i == 3)glColor3f(1,0,1);
		if(i == 4)glColor3f(1,1,0);

		for(size_t j = 0 ; j < points_local[i].size(); j++){
			Eigen::Vector3d vt;
			vt = m * points_local[i][j];
			glVertex3f(vt.x(), vt.y(), vt.z());
		}
	}
	glEnd();
	glPointSize(1);

	glColor3f(0.8,0.8,0.8);
	glBegin(GL_LINE_STRIP);
	for(size_t i = 0 ; i < trajectory.size(); i++){
		Eigen::Affine3d m = trajectory[i];
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
		case '-':{
			scale -= 0.01;
			if(scale < 0.01)scale = 0.01;
			std::cout << "scale: " << scale << std::endl;
			break;
		}
		case '=':{
			scale += 0.01;
			std::cout << "scale: " << scale << std::endl;
			break;
		}
		case 't':{
			scale += 0.0000000001;

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < points_local.size(); i++){
				for(size_t j = 0 ; j < points_local.size(); j++){
					if(i != j){

						TaitBryanPose pose_1 = pose_tait_bryan_from_affine_matrix(trajectory[i]);
						TaitBryanPose pose_2 = pose_tait_bryan_from_affine_matrix(trajectory[j]);

						for(size_t k = 0 ; k < points_local[i].size(); k++){
							Eigen::Vector3d &p_1 = points_local[i][k];
							Eigen::Vector3d &p_2 = points_local[j][k];
							double delta_x;
							double delta_y;
							double delta_z;
							point_to_point_vo_scale_tait_bryan_wc(delta_x, delta_y, delta_z, pose_1.px, pose_1.py, pose_1.pz, pose_1.om, pose_1.fi, pose_1.ka, pose_2.px, pose_2.py, pose_2.pz, pose_2.om, pose_2.fi, pose_2.ka, p_1.x(), p_1.y(), p_1.z(), p_2.x(), p_2.y(), p_2.z(), scale);

							Eigen::Matrix<double, 3, 1> jacobian;
							point_to_point_vo_scale_tait_bryan_wc_jacobian(jacobian, pose_1.px, pose_1.py, pose_1.pz, pose_2.px, pose_2.py, pose_2.pz);

							int ir = tripletListB.size();
							int ic = 0;

							tripletListA.emplace_back(ir     , ic, -jacobian(0,0));
							tripletListA.emplace_back(ir + 1 , ic, -jacobian(1,0));
							tripletListA.emplace_back(ir + 2 , ic, -jacobian(2,0));

							tripletListP.emplace_back(ir    , ir    ,  cauchy(delta_x, 1));
							tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta_y, 1));
							tripletListP.emplace_back(ir + 2, ir + 2,  cauchy(delta_z, 1));

							tripletListB.emplace_back(ir    , 0,  delta_x);
							tripletListB.emplace_back(ir + 1, 0,  delta_y);
							tripletListB.emplace_back(ir + 2, 0,  delta_z);
						}
					}
				}
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), 1);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(1, 1);
			Eigen::SparseMatrix<double> AtPB(1, 1);

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

			if(h_x.size() == 1){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				scale += h_x[0];
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}

			break;
		}
		case 'r':{
			scale += 0.0000000001;

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < points_local.size(); i++){
				for(size_t j = 0 ; j < points_local.size(); j++){
					if(i != j){
						RodriguesPose pose_1 = pose_rodrigues_from_affine_matrix(trajectory[i]);
						RodriguesPose pose_2 = pose_rodrigues_from_affine_matrix(trajectory[j]);

						pose_1.sx += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;
						pose_1.sy += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;
						pose_1.sz += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;

						pose_2.sx += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;
						pose_2.sy += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;
						pose_2.sz += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.000001;

						for(size_t k = 0 ; k < points_local[i].size(); k++){
							Eigen::Vector3d &p_1 = points_local[i][k];
							Eigen::Vector3d &p_2 = points_local[j][k];
							double delta_x;
							double delta_y;
							double delta_z;
							point_to_point_vo_scale_rodrigues_wc(delta_x, delta_y, delta_z, pose_1.px, pose_1.py, pose_1.pz, pose_1.sx, pose_1.sy, pose_1.sz, pose_2.px, pose_2.py, pose_2.pz, pose_2.sx, pose_2.sy, pose_2.sz, p_1.x(), p_1.y(), p_1.z(), p_2.x(), p_2.y(), p_2.z(), scale);

							Eigen::Matrix<double, 3, 1> jacobian;
							point_to_point_vo_scale_rodrigues_wc_jacobian(jacobian, pose_1.px, pose_1.py, pose_1.pz, pose_2.px, pose_2.py, pose_2.pz);

							int ir = tripletListB.size();
							int ic = 0;

							tripletListA.emplace_back(ir     , ic, -jacobian(0,0));
							tripletListA.emplace_back(ir + 1 , ic, -jacobian(1,0));
							tripletListA.emplace_back(ir + 2 , ic, -jacobian(2,0));

							tripletListP.emplace_back(ir    , ir    ,  cauchy(delta_x, 1));
							tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta_y, 1));
							tripletListP.emplace_back(ir + 2, ir + 2,  cauchy(delta_z, 1));

							tripletListB.emplace_back(ir    , 0,  delta_x);
							tripletListB.emplace_back(ir + 1, 0,  delta_y);
							tripletListB.emplace_back(ir + 2, 0,  delta_z);
						}
					}
				}
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), 1);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(1, 1);
			Eigen::SparseMatrix<double> AtPB(1, 1);

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

			if(h_x.size() == 1){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				scale += h_x[0];
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}

			break;
		}
		case 'q':{
			scale += 0.0000000001;

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < points_local.size(); i++){
				for(size_t j = 0 ; j < points_local.size(); j++){
					if(i != j){
						QuaternionPose pose_1 = pose_quaternion_from_affine_matrix(trajectory[i]);
						QuaternionPose pose_2 = pose_quaternion_from_affine_matrix(trajectory[j]);

						for(size_t k = 0 ; k < points_local[i].size(); k++){
							Eigen::Vector3d &p_1 = points_local[i][k];
							Eigen::Vector3d &p_2 = points_local[j][k];
							double delta_x;
							double delta_y;
							double delta_z;
							point_to_point_vo_scale_quaternion_wc(delta_x, delta_y, delta_z, pose_1.px, pose_1.py, pose_1.pz, pose_1.q0, pose_1.q1, pose_1.q2, pose_1.q3, pose_2.px, pose_2.py, pose_2.pz, pose_2.q0, pose_2.q1, pose_2.q2, pose_2.q3, p_1.x(), p_1.y(), p_1.z(), p_2.x(), p_2.y(), p_2.z(), scale);

							Eigen::Matrix<double, 3, 1> jacobian;
							point_to_point_vo_scale_rodrigues_wc_jacobian(jacobian, pose_1.px, pose_1.py, pose_1.pz, pose_2.px, pose_2.py, pose_2.pz);

							int ir = tripletListB.size();
							int ic = 0;

							tripletListA.emplace_back(ir     , ic, -jacobian(0,0));
							tripletListA.emplace_back(ir + 1 , ic, -jacobian(1,0));
							tripletListA.emplace_back(ir + 2 , ic, -jacobian(2,0));

							tripletListP.emplace_back(ir    , ir    ,  cauchy(delta_x, 1));
							tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta_y, 1));
							tripletListP.emplace_back(ir + 2, ir + 2,  cauchy(delta_z, 1));

							tripletListB.emplace_back(ir    , 0,  delta_x);
							tripletListB.emplace_back(ir + 1, 0,  delta_y);
							tripletListB.emplace_back(ir + 2, 0,  delta_z);
						}
					}
				}
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), 1);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(1, 1);
			Eigen::SparseMatrix<double> AtPB(1, 1);

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

			if(h_x.size() == 1){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				scale += h_x[0];
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
	std::cout << "r: optimize (Rodrigues)" << std::endl;
	std::cout << "q: optimize (Quaternion)" << std::endl;
	std::cout << "-: scale -= 0.01" << std::endl;
	std::cout << "=: scale += 0.01" << std::endl;
}







