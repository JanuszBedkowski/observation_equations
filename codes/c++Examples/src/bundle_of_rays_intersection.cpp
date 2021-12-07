#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "ray_intersection_observation_equation_jacobian.h"

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -10.0;
float translate_x, translate_y = 0.0;

std::vector<Eigen::Affine3d> bundle_of_rays;
Eigen::Vector3d intersection(0,0,0);

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

int main(int argc, char *argv[]){
	for(size_t i = 0; i < 10; i++){
		TaitBryanPose pose;

		pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + 5;
		pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + 5;
		pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + 5;

		pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;
		pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;
		pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;

		bundle_of_rays.push_back(affine_matrix_from_pose_tait_bryan(pose));
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
	glutCreateWindow("bundle_of_rays_intersection");
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

	glColor3f(0,1,0);
	glBegin(GL_LINES);
	for(size_t i = 0; i < bundle_of_rays.size(); i++){
		Eigen::Vector3d z_begin(0, 0,-100);
		Eigen::Vector3d z_end(0, 0, 100);

		Eigen::Vector3d z_begin_t = bundle_of_rays[i] * z_begin;
		Eigen::Vector3d z_end_t = bundle_of_rays[i] * z_end;

		glVertex3f(z_begin_t.x(), z_begin_t.y(), z_begin_t.z());
		glVertex3f(z_end_t.x(), z_end_t.y(), z_end_t.z());
	}
	glEnd();

	glPointSize(10);
	glColor3f(1,0,0);
	glBegin(GL_POINTS);
	glVertex3f(intersection.x(), intersection.y(), intersection.z());
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
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < bundle_of_rays.size(); i++){
				Eigen::Vector3d vx(bundle_of_rays[i](0,0), bundle_of_rays[i](1,0), bundle_of_rays[i](2,0));
				Eigen::Vector3d vy(bundle_of_rays[i](0,1), bundle_of_rays[i](1,1), bundle_of_rays[i](2,1));

				Eigen::Matrix<double, 2, 1> delta;
				ray_intersection_observation_equation(delta,
						intersection.x(), intersection.y(), intersection.z(),
						bundle_of_rays[i](0,3), bundle_of_rays[i](1,3), bundle_of_rays[i](2,3),
						vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z());
				Eigen::Matrix<double, 2, 3> delta_jacobian;
				ray_intersection_observation_equation_jacobian(delta_jacobian,
						intersection.x(), intersection.y(), intersection.z(),
						bundle_of_rays[i](0,3), bundle_of_rays[i](1,3), bundle_of_rays[i](2,3),
						vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z());

				int ir = tripletListB.size();
				tripletListA.emplace_back(ir, 0, -delta_jacobian(0,0));
				tripletListA.emplace_back(ir, 1, -delta_jacobian(0,1));
				tripletListA.emplace_back(ir, 2, -delta_jacobian(0,2));
				tripletListA.emplace_back(ir + 1, 0, -delta_jacobian(1,0));
				tripletListA.emplace_back(ir + 1, 1, -delta_jacobian(1,1));
				tripletListA.emplace_back(ir + 1, 2, -delta_jacobian(1,2));

				tripletListP.emplace_back(ir    , ir    ,  1);
				tripletListP.emplace_back(ir + 1, ir + 1,  1);

				tripletListB.emplace_back(ir    , 0,  delta(0,0));
				tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(3, 3);
			Eigen::SparseMatrix<double> AtPB(3, 1);

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

			if(h_x.size() == 3){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				intersection.x() += h_x[0];
				intersection.y() += h_x[1];
				intersection.z() += h_x[2];
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'n':{
			for(size_t i = 0; i < bundle_of_rays.size(); i++){
				TaitBryanPose pose;

				pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01;
				pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01;
				pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01;

				pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;
				pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;
				pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;

				bundle_of_rays[i] = bundle_of_rays[i] * affine_matrix_from_pose_tait_bryan(pose);
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
	std::cout << "n: modify rays" << std::endl;
	std::cout << "t: optimize" << std::endl;
}
