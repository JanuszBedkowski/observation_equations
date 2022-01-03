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

std::vector<double> barronminf_rho;
std::vector<double> barronminf_upsilon;
std::vector<double> barronminf_w;

std::vector<double> barronm2_rho;
std::vector<double> barronm2_upsilon;
std::vector<double> barronm2_w;

std::vector<double> barron0_rho;
std::vector<double> barron0_upsilon;
std::vector<double> barron0_w;

std::vector<double> barron15_rho;
std::vector<double> barron15_upsilon;
std::vector<double> barron15_w;

std::vector<double> barron1_rho;
std::vector<double> barron1_upsilon;
std::vector<double> barron1_w;

std::vector<double> barron32_rho;
std::vector<double> barron32_upsilon;
std::vector<double> barron32_w;

std::vector<double> barron2_rho;
std::vector<double> barron2_upsilon;
std::vector<double> barron2_w;

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

	float c = 1.0;
	for(size_t i = 0; i < r.size(); i++){
		barronminf_rho.push_back(get_barron_rho(r[i], -0.00000000001, c));
		barronminf_upsilon.push_back(get_barron_upsilon(r[i], -0.00000000001, c));
		barronminf_w.push_back(get_barron_w(r[i], -0.00000000001, c));

		barronm2_rho.push_back(get_barron_rho(r[i], -2, c));
		barronm2_upsilon.push_back(get_barron_upsilon(r[i], -2, c));
		barronm2_w.push_back(get_barron_w(r[i], -2, c));

		barron0_rho.push_back(get_barron_rho(r[i], 0, c));
		barron0_upsilon.push_back(get_barron_upsilon(r[i], 0, c));
		barron0_w.push_back(get_barron_w(r[i], 0, c));

		barron15_rho.push_back(get_barron_rho(r[i], 0.5, c));
		barron15_upsilon.push_back(get_barron_upsilon(r[i], 0.5, c));
		barron15_w.push_back(get_barron_w(r[i], 0.5, c));

		barron1_rho.push_back(get_barron_rho(r[i], 1, c));
		barron1_upsilon.push_back(get_barron_upsilon(r[i], 1, c));
		barron1_w.push_back(get_barron_w(r[i], 1, c));

		barron32_rho.push_back(get_barron_rho(r[i], 1.5, c));
		barron32_upsilon.push_back(get_barron_upsilon(r[i], 1.5, c));
		barron32_w.push_back(get_barron_w(r[i], 1.5, c));

		barron2_rho.push_back(get_barron_rho(r[i], 2, c));
		barron2_upsilon.push_back(get_barron_upsilon(r[i], 2, c));
		barron2_w.push_back(get_barron_w(r[i], 2, c));
	}

	if(argc == 2){
		std::cout << "r,/alpha=-/infty,/alpha=-2,/alpha=0,/alpha=0.5,/alpha=1,/alpha=1.5,/alpha=2" << std::endl;
		if(atoi(argv[1]) == 1){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << barronminf_rho[i] << "," << barronm2_rho[i] << "," << barron0_rho[i] << "," <<
						barron15_rho[i] << "," << barron1_rho[i] << "," << barron32_rho[i] << "," << barron2_rho[i] << std::endl;
			}
		}
		if(atoi(argv[1]) == 2){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << barronminf_upsilon[i] << "," << barronm2_upsilon[i] << "," << barron0_upsilon[i] << "," <<
						barron15_upsilon[i] << "," << barron1_upsilon[i] << "," << barron32_upsilon[i] << "," << barron2_upsilon[i] << std::endl;
			}
		}
		if(atoi(argv[1]) == 3){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << barronminf_w[i] << "," << barronm2_w[i] << "," << barron0_w[i] << "," <<
						barron15_w[i] << "," << barron1_w[i] << "," << barron32_w[i] << "," << barron2_w[i] << std::endl;
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
	glutCreateWindow("adaptive_robust_loss_function_viewer");
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
			for(size_t i = 0 ; i < barronminf_rho.size(); i++){
				glVertex3f(r[i], barronminf_rho[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barronm2_rho.size(); i++){
				glVertex3f(r[i], barronm2_rho[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron0_rho.size(); i++){
				glVertex3f(r[i], barron0_rho[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron15_rho.size(); i++){
				glVertex3f(r[i], barron15_rho[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron1_rho.size(); i++){
				glVertex3f(r[i], barron1_rho[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron32_rho.size(); i++){
				glVertex3f(r[i], barron32_rho[i], 0);
			}
			glEnd();

			glColor3f(0.6,0.2,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron2_rho.size(); i++){
				glVertex3f(r[i], barron2_rho[i], 0);
			}
			glEnd();
			break;
		}
		case RENDER_UPSILON:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barronminf_upsilon.size(); i++){
				glVertex3f(r[i], barronminf_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barronm2_upsilon.size(); i++){
				glVertex3f(r[i], barronm2_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron0_upsilon.size(); i++){
				glVertex3f(r[i], barron0_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron15_upsilon.size(); i++){
				glVertex3f(r[i], barron15_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron1_upsilon.size(); i++){
				glVertex3f(r[i], barron1_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron32_upsilon.size(); i++){
				glVertex3f(r[i], barron32_upsilon[i], 0);
			}
			glEnd();

			glColor3f(0.6,0.2,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron2_upsilon.size(); i++){
				glVertex3f(r[i], barron2_upsilon[i], 0);
			}
			glEnd();
			break;
		}
		case RENDER_W:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barronminf_w.size(); i++){
				glVertex3f(r[i], barronminf_w[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barronm2_w.size(); i++){
				glVertex3f(r[i], barronm2_w[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron0_w.size(); i++){
				glVertex3f(r[i], barron0_w[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron15_w.size(); i++){
				glVertex3f(r[i], barron15_w[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron1_w.size(); i++){
				glVertex3f(r[i], barron1_w[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron32_w.size(); i++){
				glVertex3f(r[i], barron32_w[i], 0);
			}
			glEnd();

			glColor3f(0.6,0.2,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron2_w.size(); i++){
				glVertex3f(r[i], barron2_w[i], 0);
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

