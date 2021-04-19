#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "../../example_func_x_jacobian.h"
#include "../../constraints_jacobian.h"

#define RENDER_PSI 0
#define RENDER_OBJECTIVE_FUNC 1
#define RENDER_OBJECTIVE_FUNC_WITH_CONTRAINT 2

int render_type = RENDER_PSI;


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

double x_result = 0;
double x_trg = 1;

std::vector<double> path_result;
double contraint_a = 1.0;
double contraint_x_trg = 3;
double weight_constraint = 1;
int constraint_type = 0;

void obs_eq_constraint(double &delta, double a, double x, double x_trg, int type)
{
	if(type == 0){
		observation_equation_constraint(delta, a, x, x_trg);
	}else{
		observation_equation_sq_constraint(delta, a, x, x_trg);
	}
}
void obs_eq_constraint_jacobian(Eigen::Matrix<double, 1, 1> &j, double a, double x, double x_trg, int type)
{
	if(type == 0){
		observation_equation_constraint_jacobian(j, a, x, x_trg);
	}else{
		observation_equation_sq_constraint_jacobian(j, a, x, x_trg);
	}
}

int main(int argc, char *argv[]){
	path_result.push_back(x_result);
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
	glutCreateWindow("simple_optimization_problem_func_x_with_constraints");
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

	glLineWidth(2);

	switch(render_type){
		case RENDER_PSI:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(double x = -1000; x <=1000; x+=0.01){
				double y;
				example_func_x(y,x);
				glVertex3f(x,y,0);
			}
			glEnd();

			glPointSize(20);
			glColor3f(0,0,1);
			glBegin(GL_POINTS);
				double y;
				example_func_x(y,x_result);
				glVertex3f(x_result, y, 0);
			glEnd();
			glPointSize(1);

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(const auto &x:path_result){
				double y;
				example_func_x(y,x);
				glVertex3f(x,y,0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(double x = -1000; x <=1000; x+=0.01){
				double c;
				obs_eq_constraint(c,contraint_a,x,contraint_x_trg, constraint_type);
				glVertex3f(x,c,0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINES);
				glVertex3f(-1000, x_trg, 0);
				glVertex3f(1000, x_trg, 0);
			glEnd();
		break;
		}
		case RENDER_OBJECTIVE_FUNC:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(double x = -1000; x <=1000; x+=0.01){
				double y;
				observation_equation_example_func_x(y,x,x_trg);
				glVertex3f(x,y*y,0);
			}
			glEnd();

			glPointSize(20);
			glColor3f(0,0,1);
			glBegin(GL_POINTS);
				double y;
				observation_equation_example_func_x(y,x_result,x_trg);
				glVertex3f(x_result, y*y, 0);
			glEnd();
			glPointSize(1);

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(const auto &x:path_result){
				double y;
				observation_equation_example_func_x(y,x,x_trg);
				glVertex3f(x,y*y,0);
			}
			glEnd();

			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(double x = -1000; x <=1000; x+=0.01){
				double c;
				obs_eq_constraint(c,contraint_a,x,contraint_x_trg, constraint_type);
				glVertex3f(x,c*c,0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINES);
				glVertex3f(-1000, 0, 0);
				glVertex3f(1000, 0, 0);
			glEnd();
			break;
		}
		case RENDER_OBJECTIVE_FUNC_WITH_CONTRAINT:{
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for(double x = -1000; x <=1000; x+=0.01){
				double y;
				observation_equation_example_func_x(y,x,x_trg);
				double c;
				obs_eq_constraint(c,contraint_a,x,contraint_x_trg, constraint_type);
				glVertex3f(x,y*y + c*c,0);
			}
			glEnd();

			glPointSize(20);
			glColor3f(0,0,1);
			glBegin(GL_POINTS);
				double y;
				observation_equation_example_func_x(y,x_result,x_trg);
				double c;
				obs_eq_constraint(c,contraint_a,x_result,contraint_x_trg, constraint_type);
				glVertex3f(x_result, y*y + c*c,0);
			glEnd();
			glPointSize(1);

			glColor3f(0,1,0);
			glBegin(GL_LINE_STRIP);
			for(const auto &x:path_result){
				double y;
				observation_equation_example_func_x(y,x,x_trg);
				double c;
				obs_eq_constraint(c,contraint_a,x,contraint_x_trg, constraint_type);
				glVertex3f(x,y*y + c*c,0);
			}
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINES);
				glVertex3f(-1000, 0, 0);
				glVertex3f(1000, 0, 0);
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
		case 'o':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			double delta;
			observation_equation_example_func_x(delta, x_result, x_trg);

			Eigen::Matrix<double, 1, 1> jacobian;
			observation_equation_example_func_x_jacobian(jacobian, x_result);

			int ir = 0;
			int ic = 0;
			tripletListA.emplace_back(ir, ic, -jacobian(0,0));
			tripletListP.emplace_back(ir, ir,  1);
			tripletListB.emplace_back(ir, 0,  delta);

			int number_of_columns = 1;
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

				x_result += h_x[0];
				path_result.push_back(x_result);
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'c':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			double delta;
			observation_equation_example_func_x(delta, x_result, x_trg);

			Eigen::Matrix<double, 1, 1> jacobian;
			observation_equation_example_func_x_jacobian(jacobian, x_result);

			int ir = 0;
			int ic = 0;
			tripletListA.emplace_back(ir, ic, -jacobian(0,0));
			tripletListP.emplace_back(ir, ir,  1.0);
			tripletListB.emplace_back(ir, 0,  delta);

			ir = tripletListB.size();
			obs_eq_constraint(delta, contraint_a, x_result, contraint_x_trg, constraint_type);
		    obs_eq_constraint_jacobian(jacobian, contraint_a, x_result, contraint_x_trg, constraint_type);
		    tripletListA.emplace_back(ir, ic,  -jacobian(0,0));
			tripletListP.emplace_back(ir, ir,  weight_constraint);
			tripletListB.emplace_back(ir, 0,  delta);

			int number_of_columns = 1;
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

				x_result += h_x[0];
				path_result.push_back(x_result);
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'r':{
			x_result = ((rand()%1000000)/1000000.0f - 0.5) * 2.0 * 20;
			path_result.clear();
			path_result.push_back(x_result);
			break;
		}
		case '-':{
			contraint_a -= 0.01;
			if(contraint_a < 0)contraint_a= 0.01;
			break;
		}
		case '=':{
			contraint_a += 0.01;
			break;
		}
		case 'd':{
			contraint_x_trg -= 0.1;
			break;
		}
		case 'i':{
			contraint_x_trg += 0.1;
			break;
		}
		case '1':{
			render_type = RENDER_PSI;
			break;
		}
		case '2':{
			render_type = RENDER_OBJECTIVE_FUNC;
			break;
		}
		case '3':{
			render_type = RENDER_OBJECTIVE_FUNC_WITH_CONTRAINT;
			break;
		}
		case 'y':{
			x_trg += 0.1;
			std::cout << "x_trg: " << x_trg << std::endl;
			break;
		}
		case 't':{
			x_trg -= 0.1;
			std::cout << "x_trg: " << x_trg << std::endl;
			break;
		}
		case 'q':{
			x_result-=0.1;
			path_result.clear();
			path_result.push_back(x_result);
			break;
		}
		case 'w':{
			x_result+=0.1;
			path_result.clear();
			path_result.push_back(x_result);
			break;
		}
		case '4':{
			constraint_type = 0;
			break;
		}
		case '5':{
			constraint_type = 1;
			break;
		}
		case '6':{
			weight_constraint /=10;
			std::cout << "weight_constraint: " << weight_constraint << std::endl;
			break;
		}
		case '7':{
			weight_constraint *=10;
			std::cout << "weight_constraint: " << weight_constraint << std::endl;
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
	std::cout << "c: optimize with constaint" << std::endl;
	std::cout << "-: contraint_a -= 0.01" << std::endl;
	std::cout << "=: contraint_a += 0.01" << std::endl;
	std::cout << "d: contraint_x_trg -= 0.1" << std::endl;
	std::cout << "i: contraint_x_trg += 0.1" << std::endl;
	std::cout << "1: RENDER_PSI" << std::endl;
	std::cout << "2: RENDER_OBJECTIVE_FUNC" << std::endl;
	std::cout << "3: RENDER_OBJECTIVE_FUNC_WITH_CONTRAINT" << std::endl;
	std::cout << "y: x_trg += 0.1" << std::endl;
	std::cout << "t: x_trg -= 0.1" << std::endl;
	std::cout << "q: x_result-=0.1" << std::endl;
	std::cout << "w: x_result+=0.1" << std::endl;
	std::cout << "4: constraint_type linear" << std::endl;
	std::cout << "5: constraint_type squared" << std::endl;
	std::cout << "6: weight_constraint /=10" << std::endl;
	std::cout << "7: weight_constraint *=10" << std::endl;
}







