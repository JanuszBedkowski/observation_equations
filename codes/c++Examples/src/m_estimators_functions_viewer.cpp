#include <GL/freeglut.h>
#include <iostream>
#include <vector>

#include "m_estimators.h"

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
std::vector<double> l1_rho;
std::vector<double> l1_upsilon;
std::vector<double> l1_w;

std::vector<double> l2_rho;
std::vector<double> l2_upsilon;
std::vector<double> l2_w;

std::vector<double> l1l2_rho;
std::vector<double> l1l2_upsilon;
std::vector<double> l1l2_w;

std::vector<double> lp_rho;
std::vector<double> lp_upsilon;
std::vector<double> lp_w;

std::vector<double> fair_rho;
std::vector<double> fair_upsilon;
std::vector<double> fair_w;

std::vector<double> huber_rho;
std::vector<double> huber_upsilon;
std::vector<double> huber_w;

std::vector<double> cauchy_rho;
std::vector<double> cauchy_upsilon;
std::vector<double> cauchy_w;

std::vector<double> geman_mcclure_rho;
std::vector<double> geman_mcclure_upsilon;
std::vector<double> geman_mcclure_w;

std::vector<double> welsch_rho;
std::vector<double> welsch_upsilon;
std::vector<double> welsch_w;

std::vector<double> tukey_rho;
std::vector<double> tukey_upsilon;
std::vector<double> tukey_w;

#define RENDER_RHO 1
#define RENDER_UPSILON 2
#define RENDER_W 3

int render_type = RENDER_RHO;

int main(int argc, char *argv[]){
	if(argc == 1){
		std::cout << "USAGE: " << argv[0] << " param(optional)" << std::endl;
		std::cout << "param=1 - print rho" << std::endl;
		std::cout << "param=2 - print upsilon" << std::endl;
		std::cout << "param=3 - print w" << std::endl;
	}

	for(double residuum = -10.0; residuum <= 10.0; residuum += 0.01){
		r.push_back(residuum);
	}

	for(size_t i = 0; i < r.size(); i++){
		l1_rho.push_back(get_l1_rho(r[i]));
		l1_upsilon.push_back(get_l1_upsilon(r[i]));
		l1_w.push_back(get_l1_w(r[i]));

		l2_rho.push_back(get_l2_rho(r[i]));
		l2_upsilon.push_back(get_l2_upsilon(r[i]));
		l2_w.push_back(get_l2_w(r[i]));

		l1l2_rho.push_back(get_l1l2_rho(r[i]));
		l1l2_upsilon.push_back(get_l1l2_upsilon(r[i]));
		l1l2_w.push_back(get_l1l2_w(r[i]));

		double nu = 1.2;
		lp_rho.push_back(get_lp_rho(r[i],nu));
		lp_upsilon.push_back(get_lp_upsilon(r[i],nu));
		lp_w.push_back(get_lp_w(r[i],nu));

		double c = 1.3998;
		fair_rho.push_back(get_fair_rho(r[i],c));
		fair_upsilon.push_back(get_fair_upsilon(r[i],c));
		fair_w.push_back(get_fair_w(r[i],c));

		//double k = 1.345;
		double k = 1.0;
		huber_rho.push_back(get_huber_rho(r[i],k));
		huber_upsilon.push_back(get_huber_upsilon(r[i],k));
		huber_w.push_back(get_huber_w(r[i],k));

		c = 1;
		cauchy_rho.push_back(get_cauchy_rho(r[i],c));
		cauchy_upsilon.push_back(get_cauchy_upsilon(r[i],c));
		cauchy_w.push_back(get_cauchy_w(r[i],c));

		geman_mcclure_rho.push_back(get_geman_mcclure_rho(r[i]));
		geman_mcclure_upsilon.push_back(get_geman_mcclure_upsilon(r[i]));
		geman_mcclure_w.push_back(get_geman_mcclure_w(r[i]));

		c = 1;
		welsch_rho.push_back(get_welsch_rho(r[i], c));
		welsch_upsilon.push_back(get_welsch_upsilon(r[i], c));
		welsch_w.push_back(get_welsch_w(r[i], c));

		c = 1;
		tukey_rho.push_back(get_tukey_rho(r[i], c));
		tukey_upsilon.push_back(get_tukey_upsilon(r[i], c));
		tukey_w.push_back(get_tukey_w(r[i], c));
	}

	if(argc == 2){
		std::cout << "r,L_1,L_2,L_1-L_2,L_p,Fair,Huber,Cauchy,Geman-McClure,Welsch,Tukey" << std::endl;
		if(atoi(argv[1]) == 1){
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
		}
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
	glutCreateWindow("m_estimators_functions_viewer");
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
			for(size_t i = 0 ; i < l1_rho.size(); i++){
				glVertex3f(r[i], l1_rho[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < l2_rho.size(); i++){
				glVertex3f(r[i], l2_rho[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < l1l2_rho.size(); i++){
				glVertex3f(r[i], l1l2_rho[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < lp_rho.size(); i++){
				glVertex3f(r[i], lp_rho[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < fair_rho.size(); i++){
				glVertex3f(r[i], fair_rho[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < huber_rho.size(); i++){
				glVertex3f(r[i], huber_rho[i], 0);
			}
			glEnd();

			glColor3f(0.6,0.2,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < cauchy_rho.size(); i++){
				glVertex3f(r[i], cauchy_rho[i], 0);
			}
			glEnd();

			glColor3f(0.1,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < geman_mcclure_rho.size(); i++){
				glVertex3f(r[i], geman_mcclure_rho[i], 0);
			}
			glEnd();

			glColor3f(0.1,0.1,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < welsch_rho.size(); i++){
				glVertex3f(r[i], welsch_rho[i], 0);
			}
			glEnd();

			glColor3f(0.9,0.1,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < tukey_rho.size(); i++){
				glVertex3f(r[i], tukey_rho[i], 0);
			}
			glEnd();
			break;
		}
		case RENDER_UPSILON:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < l1_upsilon.size(); i++){
				glVertex3f(r[i], l1_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < l2_upsilon.size(); i++){
				glVertex3f(r[i], l2_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < l1l2_upsilon.size(); i++){
				glVertex3f(r[i], l1l2_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < lp_upsilon.size(); i++){
				glVertex3f(r[i], lp_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < fair_upsilon.size(); i++){
				glVertex3f(r[i], fair_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < huber_upsilon.size(); i++){
				glVertex3f(r[i], huber_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0.6,0.2,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < cauchy_upsilon.size(); i++){
				glVertex3f(r[i], cauchy_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0.1,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < geman_mcclure_upsilon.size(); i++){
				glVertex3f(r[i], geman_mcclure_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0.1,0.1,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < welsch_upsilon.size(); i++){
				glVertex3f(r[i], welsch_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0.9,0.1,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < tukey_upsilon.size(); i++){
				glVertex3f(r[i], tukey_upsilon[i], 0);
			}
			glEnd();

			break;
		}
		case RENDER_W:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < l1_w.size(); i++){
				glVertex3f(r[i], l1_w[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < l2_w.size(); i++){
				glVertex3f(r[i], l2_w[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < l1l2_w.size(); i++){
				glVertex3f(r[i], l1l2_w[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < lp_w.size(); i++){
				glVertex3f(r[i], lp_w[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < fair_w.size(); i++){
				glVertex3f(r[i], fair_w[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < huber_w.size(); i++){
				glVertex3f(r[i], huber_w[i], 0);
			}
			glEnd();

			glColor3f(0.6,0.2,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < cauchy_w.size(); i++){
				glVertex3f(r[i], cauchy_w[i], 0);
			}
			glEnd();

			glColor3f(0.1,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < geman_mcclure_w.size(); i++){
				glVertex3f(r[i], geman_mcclure_w[i], 0);
			}
			glEnd();

			glColor3f(0.1,0.1,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < welsch_w.size(); i++){
				glVertex3f(r[i], welsch_w[i], 0);
			}
			glEnd();

			glColor3f(0.9,0.1,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < tukey_w.size(); i++){
				glVertex3f(r[i], tukey_w[i], 0);
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
			render_type = RENDER_UPSILON;
			break;
		}
		case '3':{
			render_type = RENDER_W;
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

