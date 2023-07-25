#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "example_func_xy_jacobian.h"
#include "constraints_jacobian.h"
#include <transformations.h>

#define RENDER_PSI 0
#define RENDER_OBJECTIVE_FUNC 1
#define RENDER_OBJECTIVE_FUNC_WITH_CONTRAINT 2

int render_type = RENDER_PSI;

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -20.0;
float translate_x, translate_y = 0.0;

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

double x_result = 5;
double y_result = 5;

double psi_trg = 0.99;

std::vector<std::pair<double, double>> path_result;
double contraint_a = 1;
double contraint_x_trg = 0;
double contraint_y_trg = 0;
double weight_constraint = 1;
int constraint_type = 0;
double angle_heading = 0;
double update_scale = 1.0;
double w_longitudal_motion = 0.01;

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
	path_result.emplace_back(x_result, y_result);

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
	glutCreateWindow("simple_optimization_problem_func_xy_with_constraints");
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

	glLineWidth(1);

	/*glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(100.0f, 0.0f, 0.0f);

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 100.0f, 0.0f);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 100.0f);
	glEnd();*/

	float scale = 10;
	switch(render_type){
		case RENDER_PSI:{
			glColor3f(1,0,0);

			float step = 0.1;
			for(double x = -10; x <=10; x+=step){
				for(double y = -10; y <=10; y+=step){
					double psi1;
					example_func_xy(psi1, x, y);

					double psi2;
					example_func_xy(psi2, x+step, y);

					double psi3;
					example_func_xy(psi3, x+step, y+step);

					double psi4;
					example_func_xy(psi4, x, y+step);

					glColor3f(psi1,1-psi1,0);
					glBegin(GL_LINE_STRIP);
						glVertex3f(x,y,psi1*scale);
						glVertex3f(x+step,y,psi2*scale);
						glVertex3f(x+step,y+step,psi3*scale);
						glVertex3f(x,y+step,psi4*scale);
						glVertex3f(x,y,psi1*scale);
					glEnd();
				}
			}

			glPointSize(20);
			glColor3f(0,0,1);
			glBegin(GL_POINTS);
				double psi;
				example_func_xy(psi,x_result,y_result);
				glVertex3f(x_result, y_result, psi*scale);
			glEnd();
			glPointSize(1);

			//heading
			glLineWidth(3);
			Eigen::Affine2d m_rot = Eigen::Affine2d::Identity();

			double angle_rad = angle_heading * M_PI/180.0;
			double sh = sin(angle_rad);
			double ch = cos(angle_rad);

			m_rot(0,0) = ch;
			m_rot(1,0) = sh;

			m_rot(0,1) = -sh;
			m_rot(1,1) =  ch;

			glColor3f(0.1,0.5,0.6);
			glBegin(GL_LINES);
				Eigen::Vector2d v(1,0);
				Eigen::Vector2d vt = m_rot * v;
				example_func_xy(psi,x_result,y_result);
				glVertex3f(x_result, y_result, psi*scale);
				glVertex3f(x_result + vt.x(), y_result + vt.y(), psi*scale);
			glEnd();
			glLineWidth(1);



			glLineWidth(2);
			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(const auto &p:path_result){
				double psi;
				example_func_xy(psi,p.first,p.second);
				glVertex3f(p.first,p.second, psi*scale);
			}
			glEnd();
			glLineWidth(1);

			glColor3f(0,0,0);
			glBegin(GL_LINES);
			for(float x = -10; x <= 10; x+=1){
				glVertex3f(-10, x, psi_trg*scale);
				glVertex3f(10, x, psi_trg*scale);

				glVertex3f(x, -10, psi_trg*scale);
				glVertex3f(x, 10, psi_trg*scale);
			}
			glEnd();
		break;
		}
		case RENDER_OBJECTIVE_FUNC:{

			glColor3f(1,0,0);
			float step = 0.1;
			for(double x = -10; x <=10; x+=step){
				for(double y = -10; y <=10; y+=step){
					double psi;
					example_func_xy(psi, x, y);
					glColor3f(psi,1-psi,0);

					glBegin(GL_LINE_STRIP);
					    double psi1;
						observation_equation_example_func_xy(psi1, x, y, psi_trg);

						double psi2;
						observation_equation_example_func_xy(psi2, x+step, y, psi_trg);

						double psi3;
						observation_equation_example_func_xy(psi3, x+step, y+step, psi_trg);

						double psi4;
						observation_equation_example_func_xy(psi4, x, y+step, psi_trg);

						glVertex3f(x,y,           psi1*psi1*scale);
						glVertex3f(x+step,y,      psi2*psi2*scale);
						glVertex3f(x+step,y+step, psi3*psi3*scale);
						glVertex3f(x,y+step,      psi4*psi4*scale);
						glVertex3f(x,y,           psi1*psi1*scale);
					glEnd();
				}
			}

			for(double x = -10; x <=10; x+=step){
				for(double y = -10; y <=10; y+=step){
					double cx_delta1;
					obs_eq_constraint(cx_delta1, contraint_a, x, contraint_x_trg, constraint_type);

					double cy_delta1;
					obs_eq_constraint(cy_delta1, contraint_a, y, contraint_y_trg, constraint_type);

					double cx_delta2;
					obs_eq_constraint(cx_delta2, contraint_a, x+step, contraint_x_trg, constraint_type);

					double cy_delta2;
					obs_eq_constraint(cy_delta2, contraint_a, y, contraint_y_trg, constraint_type);

					double cx_delta3;
					obs_eq_constraint(cx_delta3, contraint_a, x+step, contraint_x_trg, constraint_type);

					double cy_delta3;
					obs_eq_constraint(cy_delta3, contraint_a, y+step, contraint_y_trg, constraint_type);

					double cx_delta4;
					obs_eq_constraint(cx_delta4, contraint_a, x, contraint_x_trg, constraint_type);

					double cy_delta4;
					obs_eq_constraint(cy_delta4, contraint_a, y+step, contraint_y_trg, constraint_type);

					double psi;
					example_func_xy(psi, x, y);
					glColor3f(psi,1-psi,0);
					glBegin(GL_LINE_STRIP);
						glVertex3f(x,y,           cx_delta1*cx_delta1 + cy_delta1 * cy_delta1);
						glVertex3f(x+step,y,      cx_delta2*cx_delta2 + cy_delta2 * cy_delta2);
						glVertex3f(x+step,y+step, cx_delta3*cx_delta3 + cy_delta3 * cy_delta3);
						glVertex3f(x,y+step,      cx_delta4*cx_delta4 + cy_delta4 * cy_delta4);
						glVertex3f(x,y,           cx_delta1*cx_delta1 + cy_delta1 * cy_delta1);
					glEnd();
				}
			}


			glPointSize(20);
			glColor3f(0,0,1);
			glBegin(GL_POINTS);
				double psi;
				observation_equation_example_func_xy(psi, x_result, y_result, psi_trg);
				glVertex3f(x_result, y_result, psi*psi*scale);
			glEnd();
			glPointSize(1);

			//heading
			glLineWidth(3);
			Eigen::Affine2d m_rot = Eigen::Affine2d::Identity();

			double angle_rad = angle_heading * M_PI/180.0;
			double sh = sin(angle_rad);
			double ch = cos(angle_rad);

			m_rot(0,0) = ch;
			m_rot(1,0) = sh;

			m_rot(0,1) = -sh;
			m_rot(1,1) =  ch;

			glColor3f(0.1,0.5,0.6);
			glBegin(GL_LINES);
				Eigen::Vector2d v(1,0);
				Eigen::Vector2d vt = m_rot * v;
				observation_equation_example_func_xy(psi, x_result, y_result, psi_trg);
				glVertex3f(x_result, y_result, psi*psi*scale);
				glVertex3f(x_result + vt.x(), y_result + vt.y(), psi*psi*scale);
			glEnd();
			glLineWidth(1);


			glLineWidth(2);
			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(const auto &p:path_result){
				double psi;
				observation_equation_example_func_xy(psi, p.first, p.second, psi_trg);
				glVertex3f(p.first, p.second, psi*psi*scale);
			}
			glEnd();
			glLineWidth(1);
			break;
		}
		case RENDER_OBJECTIVE_FUNC_WITH_CONTRAINT:{
			glColor3f(1,0,0);
			float step = 0.1;
			for(double x = -10; x <=10; x+=step){
				for(double y = -10; y <=10; y+=step){
					double delta1;
					observation_equation_example_func_xy(delta1, x, y, psi_trg);

					double cx_delta1;
					obs_eq_constraint(cx_delta1, contraint_a, x, contraint_x_trg, constraint_type);

					double cy_delta1;
					obs_eq_constraint(cy_delta1, contraint_a, y, contraint_y_trg, constraint_type);

					double delta2;
					observation_equation_example_func_xy(delta2, x+step, y, psi_trg);

					double cx_delta2;
					obs_eq_constraint(cx_delta2, contraint_a, x+step, contraint_x_trg, constraint_type);

					double cy_delta2;
					obs_eq_constraint(cy_delta2, contraint_a, y, contraint_y_trg, constraint_type);

					double delta3;
					observation_equation_example_func_xy(delta3, x+step, y+step, psi_trg);

					double cx_delta3;
					obs_eq_constraint(cx_delta3, contraint_a, x+step, contraint_x_trg, constraint_type);

					double cy_delta3;
					obs_eq_constraint(cy_delta3, contraint_a, y+step, contraint_y_trg, constraint_type);

					double delta4;
					observation_equation_example_func_xy(delta4, x, y+step, psi_trg);

					double cx_delta4;
					obs_eq_constraint(cx_delta4, contraint_a, x, contraint_x_trg, constraint_type);

					double cy_delta4;
					obs_eq_constraint(cy_delta4, contraint_a, y+step, contraint_y_trg, constraint_type);

					double psi;
					example_func_xy(psi, x, y);
					glColor3f(psi,1-psi,0);
					glBegin(GL_LINE_STRIP);
						glVertex3f(x,y,           delta1*delta1 + cx_delta1*cx_delta1 + cy_delta1 * cy_delta1);
						glVertex3f(x+step,y,      delta2*delta2 + cx_delta2*cx_delta2 + cy_delta2 * cy_delta2);
						glVertex3f(x+step,y+step, delta3*delta3 + cx_delta3*cx_delta3 + cy_delta3 * cy_delta3);
						glVertex3f(x,y+step,      delta4*delta4 + cx_delta4*cx_delta4 + cy_delta4 * cy_delta4);
						glVertex3f(x,y,           delta1*delta1 + cx_delta1*cx_delta1 + cy_delta1 * cy_delta1);
					glEnd();
				}
			}

			glPointSize(20);
			glColor3f(0,0,1);
			glBegin(GL_POINTS);
				double delta;
				observation_equation_example_func_xy(delta, x_result, y_result, psi_trg);

				double cx_delta;
				obs_eq_constraint(cx_delta, contraint_a, x_result, contraint_x_trg, constraint_type);

				double cy_delta;
				obs_eq_constraint(cy_delta, contraint_a, y_result, contraint_y_trg, constraint_type);


				glVertex3f(x_result, y_result, delta*delta + cx_delta * cx_delta + cy_delta * cy_delta);
			glEnd();
			glPointSize(1);

			//heading
			glLineWidth(3);
			Eigen::Affine2d m_rot = Eigen::Affine2d::Identity();

			double angle_rad = angle_heading * M_PI/180.0;
			double sh = sin(angle_rad);
			double ch = cos(angle_rad);

			m_rot(0,0) = ch;
			m_rot(1,0) = sh;

			m_rot(0,1) = -sh;
			m_rot(1,1) =  ch;

			glColor3f(0.1,0.5,0.6);
			glBegin(GL_LINES);
				observation_equation_example_func_xy(delta, x_result, y_result, psi_trg);
				obs_eq_constraint(cx_delta, contraint_a, x_result, contraint_x_trg, constraint_type);
				obs_eq_constraint(cy_delta, contraint_a, y_result, contraint_y_trg, constraint_type);

				Eigen::Vector2d v(1,0);
				Eigen::Vector2d vt = m_rot * v;

				glVertex3f(x_result, y_result, delta*delta + cx_delta * cx_delta + cy_delta * cy_delta);
				glVertex3f(x_result+vt.x(), y_result+vt.y(), delta*delta + cx_delta * cx_delta + cy_delta * cy_delta);
			glEnd();
			glLineWidth(1);

			glLineWidth(2);
			glColor3f(0,0,1);
			glBegin(GL_LINE_STRIP);
			for(const auto &p:path_result){
				observation_equation_example_func_xy(delta, p.first, p.second, psi_trg);
				obs_eq_constraint(cx_delta, contraint_a, p.first, contraint_x_trg, constraint_type);
				obs_eq_constraint(cy_delta, contraint_a, p.second, contraint_y_trg, constraint_type);

				glVertex3f(p.first, p.second, delta*delta + cx_delta * cx_delta + cy_delta * cy_delta);
			}
			glEnd();
			glLineWidth(1);
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
			x_result += random(-0.0001, 0.0001);
			y_result += random(-0.0001, 0.0001);

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			double delta;
			observation_equation_example_func_xy(delta, x_result, y_result, psi_trg);

			Eigen::Matrix<double, 1, 2> jacobian;
			observation_equation_example_func_xy_jacobian(jacobian, x_result, y_result);

			if(jacobian(0,0) == 0.0 || jacobian(0,1) == 0.0 || delta == 0){
				return;
			}

			int ir = 0;
			int ic = 0;

			tripletListA.emplace_back(ir + 0, ic + 0, -jacobian(0,0));
			tripletListA.emplace_back(ir + 0, ic + 1, -jacobian(0,1));
			tripletListP.emplace_back(ir + 0, ir + 0,  1.0);
			tripletListB.emplace_back(ir + 0, 0,  delta);

			ic = 0;
			ir = tripletListB.size();
			tripletListA.emplace_back(ir, ic,  1);
			tripletListP.emplace_back(ir, ir,  0.000000001);
			tripletListB.emplace_back(ir, 0,  0);

			ic = 1;
			ir = tripletListB.size();
			tripletListA.emplace_back(ir, ic,  1);
			tripletListP.emplace_back(ir, ir,  0.000000001);
			tripletListB.emplace_back(ir, 0,  0);

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

				x_result += h_x[0] * update_scale;
				y_result += h_x[1] * update_scale;
				path_result.emplace_back(x_result, y_result);
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'c':{
			x_result += random(-0.0001, 0.0001);
			y_result += random(-0.0001, 0.0001);

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			double delta;
			observation_equation_example_func_xy(delta, x_result, y_result, psi_trg);

			Eigen::Matrix<double, 1, 2> jacobian;
			observation_equation_example_func_xy_jacobian(jacobian, x_result, y_result);

			if(jacobian(0,0) == 0.0 || jacobian(0,1) == 0.0 || delta == 0){
				return;
			}

			int ir = 0;
			int ic = 0;

			tripletListA.emplace_back(ir + 0, ic + 0, -jacobian(0,0));
			tripletListA.emplace_back(ir + 0, ic + 1, -jacobian(0,1));
			tripletListP.emplace_back(ir + 0, ir + 0,  1.0);
			tripletListB.emplace_back(ir + 0, 0,  delta);

			ic = 0;
			ir = tripletListB.size();
			double c_delta;
			Eigen::Matrix<double, 1, 1> c_jacobian;
			obs_eq_constraint(c_delta, contraint_a, x_result, contraint_x_trg, constraint_type);
			obs_eq_constraint_jacobian(c_jacobian, contraint_a, x_result, contraint_x_trg, constraint_type);
			tripletListA.emplace_back(ir, ic,  -c_jacobian(0,0));
			tripletListP.emplace_back(ir, ir,  weight_constraint);
			tripletListB.emplace_back(ir, 0,  c_delta);

			ic = 1;
			ir = tripletListB.size();
			obs_eq_constraint(c_delta, contraint_a, y_result, contraint_y_trg, constraint_type);
			obs_eq_constraint_jacobian(c_jacobian, contraint_a, y_result, contraint_y_trg, constraint_type);
			tripletListA.emplace_back(ir, ic,  -c_jacobian(0,0));
			tripletListP.emplace_back(ir, ir,  weight_constraint);
			tripletListB.emplace_back(ir, 0,  c_delta);

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

				x_result += h_x[0] * update_scale;
				y_result += h_x[1] * update_scale;
				path_result.emplace_back(x_result, y_result);
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'z':{
			x_result += random(-0.0001, 0.0001);
			y_result += random(-0.0001, 0.0001);

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			double delta;
			observation_equation_example_func_xy(delta, x_result, y_result, psi_trg);

			Eigen::Matrix<double, 1, 2> jacobian;
			observation_equation_example_func_xy_jacobian(jacobian, x_result, y_result);

			if(jacobian(0,0) == 0.0 || jacobian(0,1) == 0.0 || delta == 0){
				return;
			}

			int ir = 0;
			int ic = 0;

			tripletListA.emplace_back(ir + 0, ic + 0, -jacobian(0,0));
			tripletListA.emplace_back(ir + 0, ic + 1, -jacobian(0,1));
			tripletListP.emplace_back(ir + 0, ir + 0,  1.0);
			tripletListB.emplace_back(ir + 0, 0,  delta);

			//anizotropic motion weights
			Eigen::Matrix2d m_rot = Eigen::Matrix2d::Identity();

			double angle_rad = angle_heading * M_PI/180.0;
			double sh = sin(angle_rad);
			double ch = cos(angle_rad);

			m_rot(0,0) = ch;
			m_rot(1,0) = sh;

			m_rot(0,1) = -sh;
			m_rot(1,1) =  ch;

			Eigen::Matrix2d m_w = Eigen::Matrix2d::Identity();
			m_w(0,0) = w_longitudal_motion;
			m_w(1,1) = 1;

			Eigen::Matrix2d m_c = (m_rot*m_w)*m_rot.transpose();

			ic = 0;
			ir = tripletListB.size();
			tripletListA.emplace_back(ir, ic,  1);
			tripletListP.emplace_back(ir, ir,  m_c(0,0));
			tripletListP.emplace_back(ir, ir+1,  m_c(0,1));
			tripletListB.emplace_back(ir, 0,  0);

			ic = 1;
			ir = tripletListB.size();
			tripletListA.emplace_back(ir, ic,  1);
			tripletListP.emplace_back(ir, ir-1,  m_c(1,0));
			tripletListP.emplace_back(ir, ir,  m_c(1,1));
			tripletListB.emplace_back(ir, 0,  0);

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

				x_result += h_x[0] * update_scale;
				y_result += h_x[1] * update_scale;
				path_result.emplace_back(x_result, y_result);
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'x':{
			x_result += random(-0.0001, 0.0001);
			y_result += random(-0.0001, 0.0001);

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			double delta;
			observation_equation_example_func_xy(delta, x_result, y_result, psi_trg);

			Eigen::Matrix<double, 1, 2> jacobian;
			observation_equation_example_func_xy_jacobian(jacobian, x_result, y_result);

			if(jacobian(0,0) == 0.0 || jacobian(0,1) == 0.0 || delta == 0){
				return;
			}

			int ir = 0;
			int ic = 0;

			tripletListA.emplace_back(ir + 0, ic + 0, -jacobian(0,0));
			tripletListA.emplace_back(ir + 0, ic + 1, -jacobian(0,1));
			tripletListP.emplace_back(ir + 0, ir + 0,  1.0);
			tripletListB.emplace_back(ir + 0, 0,  delta);

			//anizotropic motion weights
			Eigen::Matrix2d m_rot = Eigen::Matrix2d::Identity();

			double angle_rad = angle_heading * M_PI/180.0;
			double sh = sin(angle_rad);
			double ch = cos(angle_rad);

			m_rot(0,0) = ch;
			m_rot(1,0) = sh;

			m_rot(0,1) = -sh;
			m_rot(1,1) =  ch;

			Eigen::Matrix2d m_w = Eigen::Matrix2d::Identity();
			m_w(0,0) = w_longitudal_motion;
			m_w(1,1) = 1;

			Eigen::Matrix2d m_c = (m_rot*m_w)*m_rot.transpose();

			ic = 0;
			ir = tripletListB.size();
			tripletListA.emplace_back(ir, ic,  1);
			tripletListP.emplace_back(ir, ir,  m_c(0,0));
			tripletListP.emplace_back(ir, ir+1,  m_c(0,1));
			tripletListB.emplace_back(ir, 0,  0);

			ic = 1;
			ir = tripletListB.size();
			tripletListA.emplace_back(ir, ic,  1);
			tripletListP.emplace_back(ir, ir-1,  m_c(1,0));
			tripletListP.emplace_back(ir, ir,  m_c(1,1));
			tripletListB.emplace_back(ir, 0,  0);

			ic = 0;
			ir = tripletListB.size();
			double c_delta;
			Eigen::Matrix<double, 1, 1> c_jacobian;
			obs_eq_constraint(c_delta, contraint_a, x_result, contraint_x_trg, constraint_type);
			obs_eq_constraint_jacobian(c_jacobian, contraint_a, x_result, contraint_x_trg, constraint_type);
			tripletListA.emplace_back(ir, ic,  -c_jacobian(0,0));
			tripletListP.emplace_back(ir, ir,  weight_constraint);
			tripletListB.emplace_back(ir, 0,  c_delta);

			ic = 1;
			ir = tripletListB.size();
			obs_eq_constraint(c_delta, contraint_a, y_result, contraint_y_trg, constraint_type);
			obs_eq_constraint_jacobian(c_jacobian, contraint_a, y_result, contraint_y_trg, constraint_type);
			tripletListA.emplace_back(ir, ic,  -c_jacobian(0,0));
			tripletListP.emplace_back(ir, ir,  weight_constraint);
			tripletListB.emplace_back(ir, 0,  c_delta);

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

				x_result += h_x[0] * update_scale;
				y_result += h_x[1] * update_scale;
				path_result.emplace_back(x_result, y_result);
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'r':{
			x_result = random(-10.0, 10.0);
			y_result = random(-10.0, 10.0);
			path_result.clear();
			path_result.emplace_back(x_result, y_result);
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
		case 'f':{
			contraint_y_trg -= 0.1;
			break;
		}
		case 'g':{
			contraint_y_trg += 0.1;
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
			psi_trg += 0.01;
			std::cout << "psi_trg: " << psi_trg << std::endl;
			break;
		}
		case 't':{
			psi_trg -= 0.01;
			std::cout << "psi_trg: " << psi_trg << std::endl;
			break;
		}
		case 'q':{
			x_result-=0.1;
			path_result.clear();
			path_result.emplace_back(x_result, y_result);
			break;
		}
		case 'w':{
			x_result+=0.1;
			path_result.clear();
			path_result.emplace_back(x_result, y_result);
			break;
		}
		case 'a':{
			y_result-=0.1;
			path_result.clear();
			path_result.emplace_back(x_result, y_result);
			break;
		}
		case 's':{
			y_result+=0.1;
			path_result.clear();
			path_result.emplace_back(x_result, y_result);
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
		case '8':{
			angle_heading += 1;
			break;
		}
		case '9':{
			angle_heading -= 1;
			break;
		}
		case '[':{
			update_scale -= 0.01;
			if(update_scale <= 0) update_scale = 0.01;
			std::cout << "update_scale: " << update_scale << std::endl;
			break;
		}
		case ']':{
			update_scale += 0.01;
			std::cout << "update_scale: " << update_scale << std::endl;
			break;
		}
		case 'k':{
			w_longitudal_motion *= 10;
			std::cout << "w_longitudal_motion: " << w_longitudal_motion << std::endl;
			break;
		}
		case 'l':{
			w_longitudal_motion /= 10;
			std::cout << "w_longitudal_motion: " << w_longitudal_motion << std::endl;
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
	std::cout << "z: optimize with anizotropic motion constaint" << std::endl;
	std::cout << "x: optimize with constaint and anizotropic motion constaint" << std::endl;
	std::cout << "-: contraint_a -= 0.01" << std::endl;
	std::cout << "=: contraint_a += 0.01" << std::endl;
	std::cout << "d: contraint_x_trg -= 0.1" << std::endl;
	std::cout << "i: contraint_x_trg += 0.1" << std::endl;
	std::cout << "f: contraint_y_trg -= 0.1" << std::endl;
	std::cout << "g: contraint_y_trg += 0.1" << std::endl;
	std::cout << "1: RENDER_PSI" << std::endl;
	std::cout << "2: RENDER_OBJECTIVE_FUNC" << std::endl;
	std::cout << "3: RENDER_OBJECTIVE_FUNC_WITH_CONTRAINT" << std::endl;
	std::cout << "y: psi_trg += 0.01" << std::endl;
	std::cout << "t: psi_trg -= 0.01" << std::endl;
	std::cout << "q: x_result-=0.1" << std::endl;
	std::cout << "w: x_result+=0.1" << std::endl;
	std::cout << "a: y_result-=0.1" << std::endl;
	std::cout << "s: y_result+=0.1" << std::endl;
	std::cout << "4: constraint_type linear" << std::endl;
	std::cout << "5: constraint_type squared" << std::endl;
	std::cout << "6: weight_constraint /=10" << std::endl;
	std::cout << "7: weight_constraint *=10" << std::endl;
	std::cout << "8: angle_heading += 1" << std::endl;
	std::cout << "9: angle_heading -= 1" << std::endl;
	std::cout << "[: pdate_scale -= 0.01" << std::endl;
	std::cout << "]: pdate_scale += 0.01" << std::endl;
	std::cout << "k: w_longitudal_motion *= 10" << std::endl;
	std::cout << "l: w_longitudal_motion /= 10" << std::endl;
	std::cout << "r: random initial guess" << std::endl;
}







