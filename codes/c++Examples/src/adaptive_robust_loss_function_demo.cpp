#include <GL/freeglut.h>
#include <iostream>
#include <vector>
#include <Eigen/Eigen>

#include "m_estimators.h"
#include "example_func_ax_plus_b_eq_y_jacobian.h"

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

double a = 0.3;
double b = 1.0;

double a_barron = 0.3;
double b_barron = 1.0;

double a_cauchy = 0.3;
double b_cauchy = 1.0;

double a_huber = 0.3;
double b_huber = 1.0;

//double barron_c = 0.25;
double barron_c = 0.05;
double barron_alpha = 1;

std::vector<std::pair<double, double>> input_data;

int main(int argc, char *argv[]){
	double y;
	for(double x = -10; x <= 10; x+= 0.01){
		example_func_ax_plus_b(y,  x,  a,  b);
		y += ((double(rand()%1000000)/1000000.0) - 0.5) * 2.0 * 0.1;
		input_data.emplace_back(x,y);
	}

	//outliers
	for(double x = -10; x <= 0; x+= 0.02){
		example_func_ax_plus_b(y,  x,  a,  b);
		y += ((double(rand()%1000000)/1000000.0) - 0.5) * 2.0 + 5;
		input_data.emplace_back(x,y);
	}
	for(double x = 0; x <= 10; x+= 0.02){
		example_func_ax_plus_b(y,  x,  a,  b);
		y += ((double(rand()%1000000)/1000000.0) - 0.5) * 2.0 - 5;
		input_data.emplace_back(x,y);
	}

	a = 0.5;
	a_barron = 0.5;
	a_cauchy = 0.5;
	a_huber = 0.5;

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
	glutCreateWindow("adaptive_robust_loss_function_demo");
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

	glPointSize(4);
	glLineWidth(5);
	glColor3f(1,0,0);
	glBegin(GL_LINE_STRIP);
	double y;
	for(double x = -10; x < 10; x+=0.01){
		example_func_ax_plus_b(y,  x,  a_barron,  b_barron);
		glVertex3f(x, y, 0);
	}
	glEnd();

	glColor3f(0,1,0);
	glBegin(GL_LINE_STRIP);
	for(double x = -10; x < 10; x+=0.01){
		example_func_ax_plus_b(y,  x,  a_cauchy,  b_cauchy);
		glVertex3f(x, y, 0);
	}
	glEnd();

	glColor3f(0,0,1);
	glBegin(GL_LINE_STRIP);
	for(double x = -10; x < 10; x+=0.01){
		example_func_ax_plus_b(y,  x,  a_huber,  b_huber);
		glVertex3f(x, y, 0);
	}
	glEnd();

	glColor3f(0,0,0);
	glBegin(GL_LINE_STRIP);
	for(double x = -10; x < 10; x+=0.01){
		example_func_ax_plus_b(y,  x,  a,  b);
		glVertex3f(x, y, 0);
	}
	glEnd();


	glColor3f(0.5,0.5,0.9);
	glBegin(GL_POINTS);
	for(const auto &d:input_data){
		glVertex3f(d.first, d.second, 0);
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
		case 'o':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < input_data.size() ; i++){
				double delta;
				observation_equation_example_func_ax_plus_b_eq_y(delta, input_data[i].first, input_data[i].second, a, b);

				Eigen::Matrix<double, 1, 2> jacobian;
				observation_equation_example_func_ax_plus_b_eq_y_jacobian(jacobian, input_data[i].first, input_data[i].second, a, b);

				int ir = tripletListB.size();

				tripletListA.emplace_back(ir, 0, -jacobian(0,0));
				tripletListA.emplace_back(ir, 1, -jacobian(0,1));
				tripletListP.emplace_back(ir, ir,  1);
				tripletListB.emplace_back(ir, 0,  delta);

			}

			int number_of_columns = 2;
			Eigen::SparseMatrix<double> matA(tripletListB.size(), number_of_columns);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(number_of_columns, number_of_columns);
			Eigen::SparseMatrix<double> AtPB(number_of_columns, 1);

			{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP) * matA;
			AtPB = (AtP) * matB;
			}

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

			if(h_x.size() == number_of_columns){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				a += h_x[0];
				b += h_x[1];

			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'b':{
			//double c = 0.01;
			double min_sum = 1000000000.0;
			for(double alpha = -10; alpha <=2; alpha += 0.1){
				double Z_tilde = get_approximate_partition_function(-10, 10, alpha, barron_c, 100);
				double sum = 0;
				for(size_t i = 0; i < input_data.size() ; i++){
					double delta;
					observation_equation_example_func_ax_plus_b_eq_y(delta, input_data[i].first, input_data[i].second, a, b);
					sum += get_truncated_robust_kernel(delta, alpha, barron_c, Z_tilde);
				}
				if(sum < min_sum){
					min_sum = sum;
					barron_alpha = alpha;
				}
			}
			std::cout << "barron_alpha: " << barron_alpha << std::endl;

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < input_data.size() ; i++){
				double delta;
				observation_equation_example_func_ax_plus_b_eq_y(delta, input_data[i].first, input_data[i].second, a_barron, b_barron);

				Eigen::Matrix<double, 1, 2> jacobian;
				observation_equation_example_func_ax_plus_b_eq_y_jacobian(jacobian, input_data[i].first, input_data[i].second, a_barron, b_barron);

				int ir = tripletListB.size();

				tripletListA.emplace_back(ir, 0, -jacobian(0,0));
				tripletListA.emplace_back(ir, 1, -jacobian(0,1));
				tripletListP.emplace_back(ir, ir, get_barron_w(delta, barron_alpha, barron_c));
				tripletListB.emplace_back(ir, 0,  delta);
			}

			int number_of_columns = 2;
			Eigen::SparseMatrix<double> matA(tripletListB.size(), number_of_columns);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(number_of_columns, number_of_columns);
			Eigen::SparseMatrix<double> AtPB(number_of_columns, 1);

			{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP) * matA;
			AtPB = (AtP) * matB;
			}

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

			if(h_x.size() == number_of_columns){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				a_barron += h_x[0];
				b_barron += h_x[1];

			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'c':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < input_data.size() ; i++){
				double delta;
				observation_equation_example_func_ax_plus_b_eq_y(delta, input_data[i].first, input_data[i].second, a_cauchy, b_cauchy);

				Eigen::Matrix<double, 1, 2> jacobian;
				observation_equation_example_func_ax_plus_b_eq_y_jacobian(jacobian, input_data[i].first, input_data[i].second, a_cauchy, b_cauchy);

				int ir = tripletListB.size();

				tripletListA.emplace_back(ir, 0, -jacobian(0,0));
				tripletListA.emplace_back(ir, 1, -jacobian(0,1));
				tripletListP.emplace_back(ir, ir, get_barron_w(delta, 0, barron_c));
				tripletListB.emplace_back(ir, 0,  delta);
			}

			int number_of_columns = 2;
			Eigen::SparseMatrix<double> matA(tripletListB.size(), number_of_columns);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(number_of_columns, number_of_columns);
			Eigen::SparseMatrix<double> AtPB(number_of_columns, 1);

			{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP) * matA;
			AtPB = (AtP) * matB;
			}

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

			if(h_x.size() == number_of_columns){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				a_cauchy += h_x[0];
				b_cauchy += h_x[1];

			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'l':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < input_data.size() ; i++){
				double delta;
				observation_equation_example_func_ax_plus_b_eq_y(delta, input_data[i].first, input_data[i].second, a_huber, b_huber);

				Eigen::Matrix<double, 1, 2> jacobian;
				observation_equation_example_func_ax_plus_b_eq_y_jacobian(jacobian, input_data[i].first, input_data[i].second, a_huber, b_huber);

				int ir = tripletListB.size();

				tripletListA.emplace_back(ir, 0, -jacobian(0,0));
				tripletListA.emplace_back(ir, 1, -jacobian(0,1));
				tripletListP.emplace_back(ir, ir, get_barron_w(delta, 1, barron_c));
				tripletListB.emplace_back(ir, 0,  delta);
			}

			int number_of_columns = 2;
			Eigen::SparseMatrix<double> matA(tripletListB.size(), number_of_columns);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(number_of_columns, number_of_columns);
			Eigen::SparseMatrix<double> AtPB(number_of_columns, 1);

			{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP) * matA;
			AtPB = (AtP) * matB;
			}

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

			if(h_x.size() == number_of_columns){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				a_huber += h_x[0];
				b_huber += h_x[1];

			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case '-':{
			barron_alpha -= 0.1;
			std::cout << "barron_alpha " << barron_alpha << std::endl;
			break;
		}
		case '=':{
			barron_alpha += 0.1;
			std::cout << "barron_alpha " << barron_alpha << std::endl;
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
	std::cout << "'o': optimize" << std::endl;
	std::cout << "'b': optimize robust Barron" << std::endl;
	std::cout << "'c': optimize robust Cauchy" << std::endl;
	std::cout << "'l': optimize robust L1L2" << std::endl;
	std::cout << "'-': barron_alpha -= 0.1" << std::endl;
	std::cout << "'=': barron_alpha += 0.1" << std::endl;
}




