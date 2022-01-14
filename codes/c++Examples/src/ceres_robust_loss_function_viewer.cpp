#include <GL/freeglut.h>
#include <iostream>
#include <vector>
#include <cmath>

#include "ceres_robust_loss_function.h"

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

std::vector<double> r;
std::vector<double> TrivialLoss;
std::vector<double> TrivialLoss_first_derivative;
std::vector<double> TrivialLoss_second_derivative;

std::vector<double> HuberLoss;
std::vector<double> HuberLoss_first_derivative;
std::vector<double> HuberLoss_second_derivative;

std::vector<double> SoftLOneLoss;
std::vector<double> SoftLOneLoss_first_derivative;
std::vector<double> SoftLOneLoss_second_derivative;

std::vector<double> CauchyLoss;
std::vector<double> CauchyLoss_first_derivative;
std::vector<double> CauchyLoss_second_derivative;

std::vector<double> ArctanLoss;
std::vector<double> ArctanLoss_first_derivative;
std::vector<double> ArctanLoss_second_derivative;

std::vector<double> TolerantLoss;
std::vector<double> TolerantLoss_first_derivative;
std::vector<double> TolerantLoss_second_derivative;

#define RENDER_RHO 1
#define RENDER_RHO_FIRST_DERIVATIVE 2
#define RENDER_RHO_SECOND_DERIVATIVE 3

int render_type = RENDER_RHO;

int main(int argc, char *argv[]){
	if(argc == 1){
		std::cout << "USAGE: " << argv[0] << " param(optional)" << std::endl;
		std::cout << "param=1 - print rho" << std::endl;
		std::cout << "param=2 - print rho first derivative" << std::endl;
		std::cout << "param=3 - print rho second derivative" << std::endl;
	}

	for(double residuum = -10.0; residuum <= 10.0; residuum += 0.01){
		r.push_back(residuum);
	}

	for(size_t i = 0; i < r.size(); i++){
		TrivialLoss.push_back(trivial_loss(r[i]*r[i]));
		TrivialLoss_first_derivative.push_back(trivial_loss_first_derivative(r[i]*r[i]));
		TrivialLoss_second_derivative.push_back(trivial_loss_second_derivative(r[i]*r[i]));

		if(r[i]*r[i] < 1){
			HuberLoss.push_back(huber_loss_less_1(r[i]*r[i]));
			HuberLoss_first_derivative.push_back(huber_loss_less_1_first_derivative(r[i]*r[i]));
			HuberLoss_second_derivative.push_back(huber_loss_less_1_second_derivative(r[i]*r[i]));
		}else{
			HuberLoss.push_back(huber_loss_more_1(r[i]*r[i]));
			HuberLoss_first_derivative.push_back(huber_loss_more_1_first_derivative(r[i]*r[i]));
			HuberLoss_second_derivative.push_back(huber_loss_more_1_second_derivative(r[i]*r[i]));
		}

		SoftLOneLoss.push_back(soft_lone_loss(r[i]*r[i]));
		SoftLOneLoss_first_derivative.push_back(soft_lone_loss_first_derivative(r[i]*r[i]));
		SoftLOneLoss_second_derivative.push_back(soft_lone_loss_second_derivative(r[i]*r[i]));

		CauchyLoss.push_back(cauchy_loss(r[i]*r[i]));
		CauchyLoss_first_derivative.push_back(cauchy_loss_first_derivative(r[i]*r[i]));
		CauchyLoss_second_derivative.push_back(cauchy_loss_second_derivative(r[i]*r[i]));

		ArctanLoss.push_back(arctan_loss(r[i]*r[i]));
		ArctanLoss_first_derivative.push_back(arctan_loss_first_derivative(r[i]*r[i]));
		ArctanLoss_second_derivative.push_back(arctan_loss_second_derivative(r[i]*r[i]));

		TolerantLoss.push_back(tolerant_loss(r[i]*r[i],1,1));
		TolerantLoss_first_derivative.push_back(tolerant_loss_first_derivative(r[i]*r[i],1,1));
		TolerantLoss_second_derivative.push_back(tolerant_loss_second_derivative(r[i]*r[i],1,1));
	}

	if(argc == 2){
		std::cout << "r,$L_1$,$L_2$,$L_1-L_2$,$L_p$,Fair,Huber,Cauchy,Geman-McClure,Welsch,Tukey" << std::endl;
		/*if(atoi(argv[1]) == 1){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << l1_rho[i] << "," << l2_rho[i] << "," << l1l2_rho[i] << "," << lp_rho[i] << "," <<
						fair_rho[i] << "," << huber_rho[i] << "," << cauchy_rho[i] << "," <<  geman_mcclure_rho[i] << "," <<
						welsch_rho[i] << "," << tukey_rho[i] << std::endl;
			}
		}
		if(atoi(argv[1]) == 2){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << l1_upsilon[i] << "," << l2_upsilon[i] << "," << l1l2_upsilon[i] << "," << lp_upsilon[i] << "," <<
						fair_upsilon[i] << "," << huber_upsilon[i] << "," << cauchy_upsilon[i] << "," <<  geman_mcclure_upsilon[i] << "," <<
						welsch_upsilon[i] << "," << tukey_upsilon[i] << std::endl;
			}
		}
		if(atoi(argv[1]) == 3){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << l1_w[i] << "," << l2_w[i] << "," << l1l2_w[i] << "," << lp_w[i] << "," <<
						fair_w[i] << "," << huber_w[i] << "," << cauchy_w[i] << "," <<  geman_mcclure_w[i] << "," <<
						welsch_w[i] << "," << tukey_w[i] << std::endl;
			}
		}*/
		return 0;
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
	glutCreateWindow("ceres_robust_loss_function_viewer");
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


	switch(render_type){
		case RENDER_RHO:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < TrivialLoss.size(); i++){
				glVertex3f(r[i], TrivialLoss[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < HuberLoss.size(); i++){
				glVertex3f(r[i], HuberLoss[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < SoftLOneLoss.size(); i++){
				glVertex3f(r[i], SoftLOneLoss[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < CauchyLoss.size(); i++){
				glVertex3f(r[i], CauchyLoss[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < ArctanLoss.size(); i++){
				glVertex3f(r[i], ArctanLoss[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < TolerantLoss.size(); i++){
				glVertex3f(r[i], TolerantLoss[i], 0);
			}
			glEnd();
			break;
		}
		case RENDER_RHO_FIRST_DERIVATIVE:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < TrivialLoss_first_derivative.size(); i++){
				glVertex3f(r[i], TrivialLoss_first_derivative[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < HuberLoss_first_derivative.size(); i++){
				glVertex3f(r[i], HuberLoss_first_derivative[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < SoftLOneLoss_first_derivative.size(); i++){
				glVertex3f(r[i], SoftLOneLoss_first_derivative[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < CauchyLoss_first_derivative.size(); i++){
				glVertex3f(r[i], CauchyLoss_first_derivative[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < ArctanLoss_first_derivative.size(); i++){
				glVertex3f(r[i], ArctanLoss_first_derivative[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < TolerantLoss_first_derivative.size(); i++){
				glVertex3f(r[i], TolerantLoss_first_derivative[i], 0);
			}
			glEnd();
			break;
		}
		case RENDER_RHO_SECOND_DERIVATIVE:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < TrivialLoss_second_derivative.size(); i++){
				glVertex3f(r[i], TrivialLoss_second_derivative[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < HuberLoss_second_derivative.size(); i++){
				glVertex3f(r[i], HuberLoss_second_derivative[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < SoftLOneLoss_second_derivative.size(); i++){
				glVertex3f(r[i], SoftLOneLoss_second_derivative[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < CauchyLoss_second_derivative.size(); i++){
				glVertex3f(r[i], CauchyLoss_second_derivative[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < ArctanLoss_second_derivative.size(); i++){
				glVertex3f(r[i], ArctanLoss_second_derivative[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < TolerantLoss_second_derivative.size(); i++){
				glVertex3f(r[i], TolerantLoss_second_derivative[i], 0);
			}
			glEnd();
			break;
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
		case '1':{
			render_type = RENDER_RHO;
			break;
		}
		case '2':{
			render_type = RENDER_RHO_FIRST_DERIVATIVE;
			break;
		}
		case '3':{
			render_type = RENDER_RHO_SECOND_DERIVATIVE;
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
	std::cout << "1: RENDER_RHO" << std::endl;
	std::cout << "2: RENDER_UPSILON" << std::endl;
	std::cout << "3: RENDER_W" << std::endl;
}

