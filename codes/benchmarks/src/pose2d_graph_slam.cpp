#include <GL/freeglut.h>
#include <Eigen/Eigen>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"
#include "../../relative_pose_tait_bryan_wc_jacobian.h"
#include "../../relative_pose_rodrigues_wc_jacobian.h"
#include "../../relative_pose_quaternion_wc_jacobian.h"
#include "../../quaternion_constraint_jacobian.h"
#include "../../relative_pose_wc_jacobian.h"

#include "manif/SE2.h"

using manif::SE2d;
using manif::SE2Tangentd;

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

std::vector<Eigen::Affine3d> m_poses;
std::vector<std::tuple<int, int, Eigen::Affine3d>> edges_g2o;
std::vector<std::vector<double>> edges_g2o_w;

int main(int argc, char *argv[]){
	if(argc != 2){
		std::cout << "USAGE: " << argv[0] << " file.g2o" << std::endl;
		return 1;
	}

	std::ifstream g2o_file(argv[1]);

	std::string line;
    while (std::getline(g2o_file, line)) {
        std::stringstream line_stream(line);
        std::string class_element;
        line_stream >> class_element;
        if (class_element == "VERTEX_SE2"){
            int id=0;
            TaitBryanPose p;
            line_stream >> id;
            line_stream >> p.px;
            line_stream >> p.py;
            line_stream >> p.ka;

            Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(p);
            m_poses.push_back(m);
        }
        else if(class_element == "EDGE_SE2"){
            int pose_id1,pose_id2;
            line_stream>>pose_id1;
            line_stream>>pose_id2;
            TaitBryanPose p;
            line_stream >> p.px;
            line_stream >> p.py;
            line_stream >> p.ka;

            Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(p);

            edges_g2o.emplace_back(std::make_tuple(pose_id1,pose_id2, m));

            std::vector<double> w(6);
			for(size_t i = 0 ; i < 6; i++){
				line_stream >> w[i];
			}
			std::vector<double> ww(3);
			ww[0] = w[0];
			ww[1] = w[3];
			ww[2] = w[5];
			edges_g2o_w.push_back(ww);

        }
    }
    g2o_file.close();

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
	glutCreateWindow("pose2d_graph_slam");
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

    glBegin(GL_LINES);
    for(size_t i = 0; i < edges_g2o.size(); i++){
        const int pose_1 = std::get<0>(edges_g2o[i]);
        const int pose_2 = std::get<1>(edges_g2o[i]);
        glVertex3f(m_poses[pose_1](0,3), m_poses[pose_1](1,3), m_poses[pose_1](2,3) );
        glVertex3f(m_poses[pose_2](0,3), m_poses[pose_2](1,3), m_poses[pose_2](2,3) );
    }
    glEnd();

    glColor3f(1,0,1);
    glPointSize(5);
    glBegin(GL_POINTS);
    for(size_t i = 0; i < m_poses.size(); i++){
        glVertex3f(m_poses[i](0,3), m_poses[i](1,3), m_poses[i](2,3) );
    }
    glEnd();

    /*glBegin(GL_LINES);
    for(size_t i = 0; i < m_poses.size(); i++){

		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(m_poses[i](0,3), m_poses[i](1,3), m_poses[i](2,3));
		glVertex3f(m_poses[i](0,3) + m_poses[i](0,0), m_poses[i](1,3) + m_poses[i](1,0), m_poses[i](2,3) + + m_poses[i](2,0));

		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(m_poses[i](0,3), m_poses[i](1,3), m_poses[i](2,3));
		glVertex3f(m_poses[i](0,3) + m_poses[i](0,1), m_poses[i](1,3) + m_poses[i](1,1), m_poses[i](2,3) + + m_poses[i](2,1));

		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(m_poses[i](0,3), m_poses[i](1,3), m_poses[i](2,3));
		glVertex3f(m_poses[i](0,3) + m_poses[i](0,2), m_poses[i](1,3) + m_poses[i](1,2), m_poses[i](2,3) + + m_poses[i](2,2));

	}
    glEnd();*/


	glutSwapBuffers();
}

void keyboard(unsigned char key, int /*x*/, int /*y*/) {
	switch (key) {
		case (27): {
			glutDestroyWindow(glutGetWindow());
			return;
		}
		case 'n':{
			for(size_t i = 0 ; i < m_poses.size(); i++){
				TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m_poses[i]);
				pose.px += ((float(rand()%1000000))/1000000.0f - 0.5) * 1.1;
				pose.py += ((float(rand()%1000000))/1000000.0f - 0.5) * 1.1;
				pose.pz += ((float(rand()%1000000))/1000000.0f - 0.5) * 1.1;
				pose.om += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.11;
				pose.fi += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.11;
				pose.ka += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.11;
				m_poses[i] = affine_matrix_from_pose_tait_bryan(pose);
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<TaitBryanPose> poses;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_tait_bryan_from_affine_matrix(m_poses[i]));
			}

			for(size_t i = 0 ; i < edges_g2o.size(); i++){
				const int first = std::get<0>(edges_g2o[i]);
				const int second = std::get<1>(edges_g2o[i]);
				const Eigen::Affine3d& rel = std::get<2>(edges_g2o[i]);
				TaitBryanPose pose_rel = pose_tait_bryan_from_affine_matrix(rel);

				TaitBryanPose from = poses[first];
				TaitBryanPose to = poses[second];

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_tait_bryan_wc_case1(
					delta,
					poses[first].px,
					poses[first].py,
					poses[first].pz,
					poses[first].om,
					poses[first].fi,
					poses[first].ka,
					poses[second].px,
					poses[second].py,
					poses[second].pz,
					poses[second].om,
					poses[second].fi,
					poses[second].ka,
					pose_rel.px,
					pose_rel.py,
					pose_rel.pz,
					pose_rel.om,
					pose_rel.fi,
					pose_rel.ka);

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_tait_bryan_wc_case1_jacobian(jacobian,
					poses[first].px,
					poses[first].py,
					poses[first].pz,
					poses[first].om,
					poses[first].fi,
					poses[first].ka,
					poses[second].px,
					poses[second].py,
					poses[second].pz,
					poses[second].om,
					poses[second].fi,
					poses[second].ka);
				int ir = tripletListB.size();

				int ic_1 = first * 3;
				int ic_2 = second * 3;

				tripletListA.emplace_back(ir + 0, ic_1    , -jacobian(0,0));
				tripletListA.emplace_back(ir + 0, ic_1 + 1, -jacobian(0,1));
				tripletListA.emplace_back(ir + 0, ic_1 + 2, -jacobian(0,5));

				tripletListA.emplace_back(ir + 0, ic_2    , -jacobian(0,6));
				tripletListA.emplace_back(ir + 0, ic_2 + 1, -jacobian(0,7));
				tripletListA.emplace_back(ir + 0, ic_2 + 2, -jacobian(0,11));

				tripletListA.emplace_back(ir + 1, ic_1    , -jacobian(1,0));
				tripletListA.emplace_back(ir + 1, ic_1 + 1, -jacobian(1,1));
				tripletListA.emplace_back(ir + 1, ic_1 + 2, -jacobian(1,5));

				tripletListA.emplace_back(ir + 1, ic_2    , -jacobian(1,6));
				tripletListA.emplace_back(ir + 1, ic_2 + 1, -jacobian(1,7));
				tripletListA.emplace_back(ir + 1, ic_2 + 2, -jacobian(1,11));

				tripletListA.emplace_back(ir + 2, ic_1    , -jacobian(5,0));
				tripletListA.emplace_back(ir + 2, ic_1 + 1, -jacobian(5,1));
				tripletListA.emplace_back(ir + 2, ic_1 + 2, -jacobian(5,5));

				tripletListA.emplace_back(ir + 2, ic_2    , -jacobian(5,6));
				tripletListA.emplace_back(ir + 2, ic_2 + 1, -jacobian(5,7));
				tripletListA.emplace_back(ir + 2, ic_2 + 2, -jacobian(5,11));

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(5,0));

				/*if(abs(first-second)==1){
					tripletListP.emplace_back(ir ,    ir,     1000);
					tripletListP.emplace_back(ir + 1, ir + 1, 1000);
					tripletListP.emplace_back(ir + 2, ir + 2, 1000);
				}else{
					tripletListP.emplace_back(ir ,    ir,     1);
					tripletListP.emplace_back(ir + 1, ir + 1, 1);
					tripletListP.emplace_back(ir + 2, ir + 2, 0.000001);
				}*/
				tripletListP.emplace_back(ir ,    ir,     cauchy(delta(0,0),1) * edges_g2o_w[i][0]);
				tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta(1,0),1) * edges_g2o_w[i][1]);
				tripletListP.emplace_back(ir + 2, ir + 2, cauchy(delta(5,0),1) * edges_g2o_w[i][2]);

				//tripletListP.emplace_back(ir ,    ir,     edges_g2o_w[i][0]);
				//tripletListP.emplace_back(ir + 1, ir + 1, edges_g2o_w[i][1]);
				//tripletListP.emplace_back(ir + 2, ir + 2, edges_g2o_w[i][2]);

				//std::cout << delta(5,0) << " ";
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);

			tripletListP.emplace_back(ir     , ir,     1000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 1000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 1000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);

			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 3 , m_poses.size() * 3);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 3 , 1);

			{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP) * matA;
			AtPB = (AtP) * matB;
			}

			tripletListA.clear();
			tripletListP.clear();
			tripletListB.clear();
			//LM
			/*{
				for(size_t i = 0 ; i < m_poses.size() * 3; i++){
					tripletListA.emplace_back(i, i, 10);
				}
				Eigen::SparseMatrix<double> matLM(m_poses.size() * 3, m_poses.size() * 3);
				matLM.setFromTriplets(tripletListA.begin(), tripletListA.end());

				AtPA = AtPA + matLM;
			}*/

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;
			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			if(h_x.size() == 3 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.ka += h_x[counter++];
					m_poses[i] = affine_matrix_from_pose_tait_bryan(pose);
				}
				std::cout << "optimizing with tait bryan finished" << std::endl;
			}else{
				std::cout << "optimizing with tait bryan FAILED" << std::endl;
			}
			break;
		}
		case 'y':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<SE2d> X;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				TaitBryanPose p = pose_tait_bryan_from_affine_matrix(m_poses[i]);
				X.push_back(SE2d(p.px, p.py, p.ka));
			}


			for(size_t i = 0 ; i < edges_g2o.size(); i++){
				const int first = std::get<0>(edges_g2o[i]);
				const int second = std::get<1>(edges_g2o[i]);
				const Eigen::Affine3d& rel = std::get<2>(edges_g2o[i]);

				TaitBryanPose p = pose_tait_bryan_from_affine_matrix(rel);
				SE2d U = SE2d(p.px, p.py, p.ka);

				SE2Tangentd     d;
				SE2Tangentd     u;
				Eigen::Matrix<double, 3, 3>         J_d_xi, J_d_xj;

				SE2d         Xi,
							 Xj;

				Xi = X[first];
				Xj = X[second];

				d  = Xj.rminus(Xi, J_d_xj, J_d_xi);
				u = U.log();

				int ir = tripletListB.size();
				int ic_1 = first * 3;
				int ic_2 = second * 3;

				for(size_t row = 0 ; row < 3; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -J_d_xi(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -J_d_xi(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -J_d_xi(row,2));

					tripletListA.emplace_back(ir + row, ic_2    , -J_d_xj(row,0));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -J_d_xj(row,1));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -J_d_xj(row,2));
				}

				SE2Tangentd delta = d - u;

				tripletListB.emplace_back(ir,     0, delta.coeffs()(0));
				tripletListB.emplace_back(ir + 1, 0, delta.coeffs()(1));
				tripletListB.emplace_back(ir + 2, 0, delta.coeffs()(2));

				/*if(abs(first - second) == 1){
					tripletListP.emplace_back(ir ,    ir,     1000);
					tripletListP.emplace_back(ir + 1, ir + 1, 1000);
					tripletListP.emplace_back(ir + 2, ir + 2, 1000);
				}else{
					tripletListP.emplace_back(ir ,    ir,     1);
					tripletListP.emplace_back(ir + 1, ir + 1, 1);
					tripletListP.emplace_back(ir + 2, ir + 2, 0.000001);
				}*/

				tripletListP.emplace_back(ir ,    ir,     cauchy(delta.coeffs()(0),1) * edges_g2o_w[i][0]);
				tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta.coeffs()(1),1) * edges_g2o_w[i][1]);
				tripletListP.emplace_back(ir + 2, ir + 2, cauchy(delta.coeffs()(2),1) * edges_g2o_w[i][2]);
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);

			tripletListP.emplace_back(ir     , ir,     1000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 1000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 1000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);

			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 3 , m_poses.size() * 3);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 3 , 1);

			{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP) * matA;
			AtPB = (AtP) * matB;
			}

			tripletListA.clear();
			tripletListP.clear();
			tripletListB.clear();

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

			if(X.size() * 3 == h_x.size()){

				int counter = 0;
				for(size_t i = 0 ; i < X.size(); i++){
					SE2Tangentd     dx;
					dx.coeffs()(0) = h_x[counter++];
					dx.coeffs()(1) = h_x[counter++];
					dx.coeffs()(2) = h_x[counter++];
					X[i] = X[i] +  dx;
				}

				for (int i = 0 ; i < m_poses.size(); i++){
					TaitBryanPose p;
					p.px = X[i].translation()(0);
					p.py = X[i].translation()(1);
					p.ka = X[i].angle();

					m_poses[i] = affine_matrix_from_pose_tait_bryan(p);
				}
				std::cout << "h_x.size(): " << h_x.size() << std::endl;
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
			}else{
				std::cout << "h_x.size(): " << h_x.size() << std::endl;
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
	std::cout << "n: add noise to poses" << std::endl;
	std::cout << "t: optimize (Tait-Bryan)" << std::endl;
	std::cout << "y: optimize (Lie Algebra: manif lib)" << std::endl;
}
