#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "distance_to_circle_jacobian.h"

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -5.0;
float translate_x, translate_y = 0.0;

struct Circle{
	double x;
	double y;
	double r;

	Circle(){
		x = 0.0;
		y = 0.0;
		r = 1.0;
	}
};

Circle circle;
std::vector<Eigen::Vector2d> points;

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

void draw_circle(const Circle &circle, float angle_step_rad, float r, float g, float b);

int main(int argc, char *argv[]){

	Eigen::Affine3d m_center = Eigen::Affine3d::Identity();
	m_center(0,3) = 0.2;
	m_center(1,3) = 0.4;
	for(double angle = 0.0; angle <= 2.0 * M_PI; angle += 0.1){
		TaitBryanPose rot;
		rot.ka = angle;

		Eigen::Affine3d mrot = affine_matrix_from_pose_tait_bryan(rot);

		Eigen::Vector3d p(2 + ((float(rand()%1000000))/1000000.0 - 0.05) * 2.0 * 0.1, 0, 0);
		Eigen::Affine3d m = m_center * mrot;
		Eigen::Vector3d pt = m * p;

		points.push_back(Eigen::Vector2d(pt.x(), pt.y()));
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
	glutCreateWindow("distance_to_circle");
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

	draw_circle(circle, 0.01, 0, 0, 0);

	glPointSize(5);
	glColor3f(1,0,0);
	glBegin(GL_POINTS);
		for(const auto& p:points){
			glVertex2f(p.x(), p.y());
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
		case 'o':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < points.size(); i++){
				double delta;
				Eigen::Matrix<double, 1, 3, Eigen::RowMajor> jacobian;
				observation_equation_distance_to_circle(delta, points[i].x(), points[i].y(), circle.x, circle.y, circle.r);
				observation_equation_distance_to_circle_jacobian(jacobian, points[i].x(), points[i].y(), circle.x, circle.y, circle.r);

				int ir = tripletListB.size();

				tripletListA.emplace_back(ir     , 0, -jacobian(0,0));
				tripletListA.emplace_back(ir     , 1, -jacobian(0,1));
				tripletListA.emplace_back(ir     , 2, -jacobian(0,2));

				tripletListP.emplace_back(ir    , ir    ,  1);

				tripletListB.emplace_back(ir    , 0,  delta);
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

			std::cout << "results" << std::endl;
			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 3){
				std::cout << "OPTIMIZATION SUCCESS" << std::endl;
				circle.x += h_x[0];
				circle.y += h_x[1];
				circle.r += h_x[2];
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
	std::cout << "o: optimize" << std::endl;
}

void draw_circle(const Circle &circle, float angle_step_rad, float r, float g, float b){
	glColor3f(r,g,b);
	Eigen::Affine3d m_center = Eigen::Affine3d::Identity();
	m_center(0,3) = circle.x;
	m_center(1,3) = circle.y;


	glBegin(GL_LINE_STRIP);
	for(double angle = 0.0; angle <= 2.0 * M_PI; angle += angle_step_rad){
		TaitBryanPose rot;
		rot.ka = angle;

		Eigen::Affine3d mrot = affine_matrix_from_pose_tait_bryan(rot);

		Eigen::Vector3d p(circle.r, 0,0);
		Eigen::Affine3d m = m_center * mrot;
		Eigen::Vector3d pt = m * p;

		glVertex3f(pt.x(), pt.y(), pt.z());

	}
	glEnd();
}





