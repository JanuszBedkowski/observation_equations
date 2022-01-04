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
std::vector<double> barronminf_mpdf;
std::vector<double> barronminf_rho_a_tilde;

std::vector<double> barronm2_rho;
std::vector<double> barronm2_upsilon;
std::vector<double> barronm2_w;
std::vector<double> barronm2_mpdf;
std::vector<double> barronm2_rho_a_tilde;

std::vector<double> barron0_rho;
std::vector<double> barron0_upsilon;
std::vector<double> barron0_w;
std::vector<double> barron0_mpdf;
std::vector<double> barron0_rho_a_tilde;

std::vector<double> barron05_rho;
std::vector<double> barron05_upsilon;
std::vector<double> barron05_w;
std::vector<double> barron05_mpdf;
std::vector<double> barron05_rho_a_tilde;

std::vector<double> barron1_rho;
std::vector<double> barron1_upsilon;
std::vector<double> barron1_w;
std::vector<double> barron1_mpdf;
std::vector<double> barron1_rho_a_tilde;

std::vector<double> barron15_rho;
std::vector<double> barron15_upsilon;
std::vector<double> barron15_w;
std::vector<double> barron15_mpdf;
std::vector<double> barron15_rho_a_tilde;

std::vector<double> barron2_rho;
std::vector<double> barron2_upsilon;
std::vector<double> barron2_w;
std::vector<double> barron2_mpdf;
std::vector<double> barron2_rho_a_tilde;

#define RENDER_RHO 1
#define RENDER_UPSILON 2
#define RENDER_W 3
#define RENDER_P_TILDE 4
#define RENDER_RHO_A_TILDE 5

int render_type = RENDER_RHO;

int main(int argc, char *argv[]){
	if(argc == 1){
		std::cout << "USAGE: " << argv[0] << " param(optional)" << std::endl;
		std::cout << "optional params" << std::endl;
		std::cout << "'1' print rho" << std::endl;
		std::cout << "'2' print upsilon" << std::endl;
		std::cout << "'3' print w" << std::endl;
		std::cout << "'4' print mpdf" << std::endl;
		std::cout << "'5' print rho_a_tide" << std::endl;
	}

	for(double residuum = -10.0; residuum <= 10.0; residuum += 0.01){
		r.push_back(residuum);
	}

	double minf = -1000000000.0;

	float c = 1.0;
	float num_steps = 1000000;
	double Z_tilde_minf = get_approximate_partition_function(-10.0, 10.0, minf, c, num_steps);
	double Z_tilde_m2   = get_approximate_partition_function(-10.0, 10.0, -2  , c, num_steps);
	double Z_tilde_0    = get_approximate_partition_function(-10.0, 10.0, 0   , c, num_steps);
	double Z_tilde_05   = get_approximate_partition_function(-10.0, 10.0, 0.5 , c, num_steps);
	double Z_tilde_1    = get_approximate_partition_function(-10.0, 10.0, 1.0 , c, num_steps);
	double Z_tilde_15   = get_approximate_partition_function(-10.0, 10.0, 1.5 , c, num_steps);
	double Z_tilde_2    = get_approximate_partition_function(-10.0, 10.0, 2.0 , c, num_steps);


	for(size_t i = 0; i < r.size(); i++){
		barronminf_rho.push_back(get_barron_rho(r[i], minf, c));
		barronminf_upsilon.push_back(get_barron_upsilon(r[i], minf, c));
		barronminf_w.push_back(get_barron_w(r[i], minf, c));
		barronminf_mpdf.push_back(1/(c*Z_tilde_minf)*exp(-get_barron_rho(r[i], minf, c)));
		barronminf_rho_a_tilde.push_back(get_truncated_robust_kernel(r[i], minf, c, Z_tilde_minf));

		barronm2_rho.push_back(get_barron_rho(r[i], -2, c));
		barronm2_upsilon.push_back(get_barron_upsilon(r[i], -2, c));
		barronm2_w.push_back(get_barron_w(r[i], -2, c));
		barronm2_mpdf.push_back(1/(c*Z_tilde_m2)*exp(-get_barron_rho(r[i], -2, c)));
		barronm2_rho_a_tilde.push_back(get_truncated_robust_kernel(r[i], -2, c, Z_tilde_m2));

		barron0_rho.push_back(get_barron_rho(r[i], 0, c));
		barron0_upsilon.push_back(get_barron_upsilon(r[i], 0, c));
		barron0_w.push_back(get_barron_w(r[i], 0, c));
		barron0_mpdf.push_back(1/(c*Z_tilde_0)*exp(-get_barron_rho(r[i], 0, c)));
		barron0_rho_a_tilde.push_back(get_truncated_robust_kernel(r[i], 0, c, Z_tilde_0));

		barron05_rho.push_back(get_barron_rho(r[i], 0.5, c));
		barron05_upsilon.push_back(get_barron_upsilon(r[i], 0.5, c));
		barron05_w.push_back(get_barron_w(r[i], 0.5, c));
		barron05_mpdf.push_back(1/(c*Z_tilde_05)*exp(-get_barron_rho(r[i], 0.5, c)));
		barron05_rho_a_tilde.push_back(get_truncated_robust_kernel(r[i], 0.5, c, Z_tilde_05));

		barron1_rho.push_back(get_barron_rho(r[i], 1, c));
		barron1_upsilon.push_back(get_barron_upsilon(r[i], 1, c));
		barron1_w.push_back(get_barron_w(r[i], 1, c));
		barron1_mpdf.push_back(1/(c*Z_tilde_1)*exp(-get_barron_rho(r[i], 1.0, c)));
		barron1_rho_a_tilde.push_back(get_truncated_robust_kernel(r[i], 1.0, c, Z_tilde_1));

		barron15_rho.push_back(get_barron_rho(r[i], 1.5, c));
		barron15_upsilon.push_back(get_barron_upsilon(r[i], 1.5, c));
		barron15_w.push_back(get_barron_w(r[i], 1.5, c));
		barron15_mpdf.push_back(1/(c*Z_tilde_15)*exp(-get_barron_rho(r[i], 1.5, c)));
		barron15_rho_a_tilde.push_back(get_truncated_robust_kernel(r[i], 1.5, c, Z_tilde_15));

		barron2_rho.push_back(get_barron_rho(r[i], 2, c));
		barron2_upsilon.push_back(get_barron_upsilon(r[i], 2, c));
		barron2_w.push_back(get_barron_w(r[i], 2, c));
		barron2_mpdf.push_back(1/(c*Z_tilde_2)*exp(-get_barron_rho(r[i], 2.0, c)));
		barron2_rho_a_tilde.push_back(get_truncated_robust_kernel(r[i], 2.0, c, Z_tilde_2));

	}

	if(argc == 2){
		std::cout << "r,$\\alpha=-\\infty$,$\\alpha=-2$,$\\alpha=0$,$\\alpha=0.5$,$\\alpha=1$,$\\alpha=1.5$,$\\alpha=2$" << std::endl;
		if(atoi(argv[1]) == 1){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << barronminf_rho[i] << "," << barronm2_rho[i] << "," << barron0_rho[i] << "," <<
						barron05_rho[i] << "," << barron1_rho[i] << "," << barron15_rho[i] << "," << barron2_rho[i] << std::endl;
			}
		}
		if(atoi(argv[1]) == 2){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << barronminf_upsilon[i] << "," << barronm2_upsilon[i] << "," << barron0_upsilon[i] << "," <<
						barron05_upsilon[i] << "," << barron1_upsilon[i] << "," << barron15_upsilon[i] << "," << barron2_upsilon[i] << std::endl;
			}
		}
		if(atoi(argv[1]) == 3){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << barronminf_w[i] << "," << barronm2_w[i] << "," << barron0_w[i] << "," <<
						barron05_w[i] << "," << barron1_w[i] << "," << barron15_w[i] << "," << barron2_w[i] << std::endl;
			}
		}
		if(atoi(argv[1]) == 4){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << barronminf_mpdf[i] << "," << barronm2_mpdf[i] << "," << barron0_mpdf[i] << "," <<
						barron05_mpdf[i] << "," << barron1_mpdf[i] << "," << barron15_mpdf[i] << "," << barron2_mpdf[i] << std::endl;
			}
		}

		if(atoi(argv[1]) == 5){
			for(size_t i = 0; i < r.size(); i++){
				std::cout << r[i] << "," << barronminf_rho_a_tilde[i] << "," << barronm2_rho_a_tilde[i] << "," << barron0_rho_a_tilde[i] << "," <<
						barron05_rho_a_tilde[i] << "," << barron1_rho_a_tilde[i] << "," << barron15_rho_a_tilde[i] << "," << barron2_rho_a_tilde[i] << std::endl;
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
			for(size_t i = 0 ; i < barron05_rho.size(); i++){
				glVertex3f(r[i], barron05_rho[i], 0);
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
			for(size_t i = 0 ; i < barron15_rho.size(); i++){
				glVertex3f(r[i], barron15_rho[i], 0);
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
			for(size_t i = 0 ; i < barron05_upsilon.size(); i++){
				glVertex3f(r[i], barron05_upsilon[i], 0);
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
			for(size_t i = 0 ; i < barron15_upsilon.size(); i++){
				glVertex3f(r[i], barron15_upsilon[i], 0);
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
			for(size_t i = 0 ; i < barron05_w.size(); i++){
				glVertex3f(r[i], barron05_w[i], 0);
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
			for(size_t i = 0 ; i < barron15_w.size(); i++){
				glVertex3f(r[i], barron15_w[i], 0);
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
		case RENDER_P_TILDE:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barronminf_mpdf.size(); i++){
				glVertex3f(r[i], barronminf_mpdf[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barronm2_mpdf.size(); i++){
				glVertex3f(r[i], barronm2_mpdf[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron0_mpdf.size(); i++){
				glVertex3f(r[i], barron0_mpdf[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron05_mpdf.size(); i++){
				glVertex3f(r[i], barron05_mpdf[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron1_mpdf.size(); i++){
				glVertex3f(r[i], barron1_mpdf[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron15_mpdf.size(); i++){
				glVertex3f(r[i], barron15_mpdf[i], 0);
			}
			glEnd();

			glColor3f(0.6,0.2,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron2_mpdf.size(); i++){
				glVertex3f(r[i], barron2_mpdf[i], 0);
			}
			glEnd();
			break;
		}
		case RENDER_RHO_A_TILDE:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barronminf_rho_a_tilde.size(); i++){
				glVertex3f(r[i], barronminf_rho_a_tilde[i], 0);
			}
			glEnd();

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barronm2_rho_a_tilde.size(); i++){
				glVertex3f(r[i], barronm2_rho_a_tilde[i], 0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron0_rho_a_tilde.size(); i++){
				glVertex3f(r[i], barron0_rho_a_tilde[i], 0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron05_rho_a_tilde.size(); i++){
				glVertex3f(r[i], barron05_rho_a_tilde[i], 0);
			}
			glEnd();

			glColor3f(0.5,0.9,0);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron1_rho_a_tilde.size(); i++){
				glVertex3f(r[i], barron1_rho_a_tilde[i], 0);
			}
			glEnd();

			glColor3f(0.0,0.6,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron15_rho_a_tilde.size(); i++){
				glVertex3f(r[i], barron15_rho_a_tilde[i], 0);
			}
			glEnd();

			glColor3f(0.6,0.2,0.4);
			glBegin(GL_LINE_STRIP);
			for(size_t i = 0 ; i < barron2_rho_a_tilde.size(); i++){
				glVertex3f(r[i], barron2_rho_a_tilde[i], 0);
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
		case '4':{
			render_type = RENDER_P_TILDE;
			break;
		}
		case '5':{
			render_type = RENDER_RHO_A_TILDE;
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
	std::cout << "4: RENDER_P_TILDE" << std::endl;
	std::cout << "5: RENDER_RHO_A_TILDE" << std::endl;
}




