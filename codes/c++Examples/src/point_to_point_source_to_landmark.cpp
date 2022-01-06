#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>
#include <Eigen/SparseQR>

#include "structures.h"
#include "transformations.h"
#include "point_to_point_source_to_target_tait_bryan_wc_jacobian.h"
#include "point_to_point_source_to_target_rodrigues_wc_jacobian.h"
#include "point_to_point_source_to_target_quaternion_wc_jacobian.h"
#include "quaternion_constraint_jacobian.h"
#include "cauchy.h"
#include "point_to_point_source_to_landmark_tait_bryan_wc_jacobian.h"
#include "point_to_point_source_to_landmark_rodrigues_wc_jacobian.h"
#include "point_to_point_source_to_landmark_quaternion_wc_jacobian.h"

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -100.0;
float translate_x, translate_y = 0.0;

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

struct Measurement{
	Eigen::Vector3d value;
	int index_landmark;
};

struct Node{
	Eigen::Affine3d pose;
	std::vector<Measurement> measurements;
};

std::vector<Eigen::Vector3d> landmarks;
std::vector<Node> nodes;
std::vector<Node> nodesInitial;

int main(int argc, char *argv[]){

	for(size_t i = 0; i < 10; i++){
		landmarks.emplace_back(Eigen::Vector3d((rand()%1000 - 500)*0.1, (rand()%1000 - 500)*0.1, (rand()%1000 - 500)*0.0001));
	}

	for(size_t i = 0; i < 10; i++){
		Eigen::Affine3d m = Eigen::Affine3d::Identity();
		m(0,3) = (rand()%1000 - 500)*0.001;
		m(1,3) = i*5;
		m(2,3) = 1;
		Node n;
		n.pose = m;

		TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m);
		pose.ka = M_PI /4.0;

		n.pose = affine_matrix_from_pose_tait_bryan(pose);

		nodes.emplace_back(n);
	}


	nodesInitial = nodes;

	for(size_t i = 0 ; i < nodes.size(); i++){
		for(size_t j = 0; j < landmarks.size(); j++){
			Measurement meas;
			meas.index_landmark = j;
			meas.value = nodes[i].pose.inverse() * landmarks[j];
			nodes[i].measurements.emplace_back(meas);
		}
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
	glutCreateWindow("point to point source to landmark");
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


	glColor3f(0, 1, 0);
	glBegin(GL_LINES);
	for(const auto &l:landmarks){
		glVertex3f(l.x() - 1, l.y(), l.z());
		glVertex3f(l.x() + 1, l.y(), l.z());

		glVertex3f(l.x(), l.y() - 1, l.z());
		glVertex3f(l.x(), l.y() + 1, l.z());
	}
	glEnd();

	glColor3f(1, 0, 0);
	glBegin(GL_LINE_STRIP);
		for(int i = 0 ; i < nodes.size(); i++){
			glVertex3f(nodes[i].pose(0,3), nodes[i].pose(1,3), nodes[i].pose(2,3));
		}
	glEnd();

	glColor3f(0,0,0);
	glBegin(GL_LINES);
		for(int i = 0 ; i < nodesInitial.size(); i++){
			glVertex3f(nodesInitial[i].pose(0,3)-1, nodesInitial[i].pose(1,3), nodesInitial[i].pose(2,3));
			glVertex3f(nodesInitial[i].pose(0,3)+1, nodesInitial[i].pose(1,3), nodesInitial[i].pose(2,3));

			glVertex3f(nodesInitial[i].pose(0,3), nodesInitial[i].pose(1,3) - 1, nodesInitial[i].pose(2,3));
			glVertex3f(nodesInitial[i].pose(0,3), nodesInitial[i].pose(1,3) + 1, nodesInitial[i].pose(2,3));
		}
	glEnd();

	glColor3f(0.6,0.6,0.6);
	glBegin(GL_LINES);
		for(int i = 0 ; i < nodesInitial.size(); i++){
			glVertex3f(nodesInitial[i].pose(0,3), nodesInitial[i].pose(1,3), nodesInitial[i].pose(2,3));
			glVertex3f(nodes[i].pose(0,3), nodes[i].pose(1,3), nodes[i].pose(2,3));
		}
	glEnd();




	glColor3f(0, 0, 1);
	glBegin(GL_LINES);

	for(size_t i = 0 ; i < nodes.size(); i++){
		for(size_t j = 0 ; j < nodes[i].measurements.size(); j++){
			Eigen::Vector3d meas = nodes[i].pose * nodes[i].measurements[j].value;

			glVertex3f(meas.x() - 0.1, meas.y(), meas.z());
			glVertex3f(meas.x() + 0.1, meas.y(), meas.z());

			glVertex3f(meas.x(), meas.y() - 0.1, meas.z());
			glVertex3f(meas.x(), meas.y() + 0.1, meas.z());
		}
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
		case 'n':{
			for(size_t i = 0 ; i < nodes.size(); i++){
				TaitBryanPose pose;
						pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.5;
						pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.5;
						pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.0000005;
						pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.0000005;
						pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.0000005;
						pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.005;
				nodes[i].pose = nodes[i].pose * affine_matrix_from_pose_tait_bryan(pose);
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < nodes.size(); i++){
				Eigen::Affine3d m = nodes[i].pose;
				TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m);

				for(size_t j = 0 ; j < nodes[i].measurements.size(); j++){

					Eigen::Vector3d &p_t = landmarks[nodes[i].measurements[j].index_landmark];
					Eigen::Vector3d &p_s = nodes[i].measurements[j].value;

					double delta_x;
					double delta_y;
					double delta_z;
					point_to_point_source_to_landmark_tait_bryan_wc(delta_x, delta_y, delta_z,
							pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
							p_s.x(), p_s.y(), p_s.z(), p_t.x(), p_t.y(), p_t.z());

					Eigen::Matrix<double, 3, 9, Eigen::RowMajor> jacobian;
					point_to_point_source_to_landmark_tait_bryan_wc_jacobian(jacobian,
							pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
							p_s.x(), p_s.y(), p_s.z());


					int ir = tripletListB.size();
					int ic = 6 * i;

					if(jacobian(0,0) != 0.0)tripletListA.emplace_back(ir + 0, ic + 0, -jacobian(0,0));
					if(jacobian(0,1) != 0.0)tripletListA.emplace_back(ir + 0, ic + 1, -jacobian(0,1));
					if(jacobian(0,2) != 0.0)tripletListA.emplace_back(ir + 0, ic + 2, -jacobian(0,2));
					if(jacobian(0,3) != 0.0)tripletListA.emplace_back(ir + 0, ic + 3, -jacobian(0,3));
					if(jacobian(0,4) != 0.0)tripletListA.emplace_back(ir + 0, ic + 4, -jacobian(0,4));
					if(jacobian(0,5) != 0.0)tripletListA.emplace_back(ir + 0, ic + 5, -jacobian(0,5));

					ic = nodes.size() * 6 + nodes[i].measurements[j].index_landmark * 3;
					if(jacobian(0,6) != 0.0)tripletListA.emplace_back(ir + 0, ic + 0, -jacobian(0,6));
					if(jacobian(0,7) != 0.0)tripletListA.emplace_back(ir + 0, ic + 1, -jacobian(0,7));
					if(jacobian(0,8) != 0.0)tripletListA.emplace_back(ir + 0, ic + 2, -jacobian(0,8));

					ic = 6 * i;
					if(jacobian(1,0) != 0.0)tripletListA.emplace_back(ir + 1, ic + 0, -jacobian(1,0));
					if(jacobian(1,1) != 0.0)tripletListA.emplace_back(ir + 1, ic + 1, -jacobian(1,1));
					if(jacobian(1,2) != 0.0)tripletListA.emplace_back(ir + 1, ic + 2, -jacobian(1,2));
					if(jacobian(1,3) != 0.0)tripletListA.emplace_back(ir + 1, ic + 3, -jacobian(1,3));
					if(jacobian(1,4) != 0.0)tripletListA.emplace_back(ir + 1, ic + 4, -jacobian(1,4));
					if(jacobian(1,5) != 0.0)tripletListA.emplace_back(ir + 1, ic + 5, -jacobian(1,5));

					ic = nodes.size() * 6 + nodes[i].measurements[j].index_landmark * 3;
					if(jacobian(1,6) != 0.0)tripletListA.emplace_back(ir + 1, ic + 0, -jacobian(1,6));
					if(jacobian(1,7) != 0.0)tripletListA.emplace_back(ir + 1, ic + 1, -jacobian(1,7));
					if(jacobian(1,8) != 0.0)tripletListA.emplace_back(ir + 1, ic + 2, -jacobian(1,8));

					ic = 6 * i;
					if(jacobian(2,0) != 0.0)tripletListA.emplace_back(ir + 2, ic + 0, -jacobian(2,0));
					if(jacobian(2,1) != 0.0)tripletListA.emplace_back(ir + 2, ic + 1, -jacobian(2,1));
					if(jacobian(2,2) != 0.0)tripletListA.emplace_back(ir + 2, ic + 2, -jacobian(2,2));
					if(jacobian(2,3) != 0.0)tripletListA.emplace_back(ir + 2, ic + 3, -jacobian(2,3));
					if(jacobian(2,4) != 0.0)tripletListA.emplace_back(ir + 2, ic + 4, -jacobian(2,4));
					if(jacobian(2,5) != 0.0)tripletListA.emplace_back(ir + 2, ic + 5, -jacobian(2,5));

					ic = nodes.size() * 6 + nodes[i].measurements[j].index_landmark * 3;
					if(jacobian(2,6) != 0.0)tripletListA.emplace_back(ir + 2, ic + 0, -jacobian(2,6));
					if(jacobian(2,7) != 0.0)tripletListA.emplace_back(ir + 2, ic + 1, -jacobian(2,7));
					if(jacobian(2,8) != 0.0)tripletListA.emplace_back(ir + 2, ic + 2, -jacobian(2,8));

					tripletListP.emplace_back(ir    , ir    , 1);
					tripletListP.emplace_back(ir + 1, ir + 1, 1);
					tripletListP.emplace_back(ir + 2, ir + 2, 1);

					tripletListB.emplace_back(ir    , 0,  delta_x);
					tripletListB.emplace_back(ir + 1, 0,  delta_y);
					tripletListB.emplace_back(ir + 2, 0,  delta_z);
				}
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);
			tripletListA.emplace_back(ir + 3 , 3, 1);
			tripletListA.emplace_back(ir + 4 , 4, 1);
			tripletListA.emplace_back(ir + 5 , 5, 1);

			tripletListP.emplace_back(ir     , ir,     10000000000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 10000000000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 10000000000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 10000000000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 10000000000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 10000000000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);

			Eigen::SparseMatrix<double> matA(tripletListB.size(), nodes.size() * 6 + landmarks.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(nodes.size() * 6 + landmarks.size() * 3, nodes.size() * 6 + landmarks.size() * 3);
			Eigen::SparseMatrix<double> AtPB(nodes.size() * 6 + landmarks.size() * 3, 1);

			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = AtP * matA;
			AtPB = AtP * matB;

			tripletListA.clear();
			tripletListP.clear();
			tripletListB.clear();

			Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(AtPA);
			Eigen::SparseMatrix<double> x = solver.solve(AtPB);

			std::vector<double> h_x;

			for (int k=0; k<x.outerSize(); ++k){
				for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it){
					std::cout << it.row() << " " << it.col() << " " << it.value() << std::endl;
					h_x.push_back(it.value());
				}
			}

			if(h_x.size() == nodes.size() * 6 + landmarks.size() * 3){
				std::cout << "optimization success" << std::endl;
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}

				int counter = 0;

				for(size_t i = 0 ; i < nodes.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(nodes[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];
					nodes[i].pose = affine_matrix_from_pose_tait_bryan(pose);
				}

				for(size_t i = 0 ; i < landmarks.size(); i++){
					landmarks[i].x() += h_x[counter++];
					landmarks[i].y() += h_x[counter++];
					landmarks[i].z() += h_x[counter++];
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'r':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < nodes.size(); i++){
				Eigen::Affine3d m = nodes[i].pose;
				RodriguesPose pose = pose_rodrigues_from_affine_matrix(m);

				for(size_t j = 0 ; j < nodes[i].measurements.size(); j++){

					Eigen::Vector3d &p_t = landmarks[nodes[i].measurements[j].index_landmark];
					Eigen::Vector3d &p_s = nodes[i].measurements[j].value;

					double delta_x;
					double delta_y;
					double delta_z;
					point_to_point_source_to_landmark_rodrigues_wc(delta_x, delta_y, delta_z,
							pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
							p_s.x(), p_s.y(), p_s.z(), p_t.x(), p_t.y(), p_t.z());

					Eigen::Matrix<double, 3, 9, Eigen::RowMajor> jacobian;
					point_to_point_source_to_landmark_rodrigues_wc_jacobian(jacobian,
							pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
							p_s.x(), p_s.y(), p_s.z());


					int ir = tripletListB.size();
					int ic = 6 * i;

					if(jacobian(0,0) != 0.0)tripletListA.emplace_back(ir + 0, ic + 0, -jacobian(0,0));
					if(jacobian(0,1) != 0.0)tripletListA.emplace_back(ir + 0, ic + 1, -jacobian(0,1));
					if(jacobian(0,2) != 0.0)tripletListA.emplace_back(ir + 0, ic + 2, -jacobian(0,2));
					if(jacobian(0,3) != 0.0)tripletListA.emplace_back(ir + 0, ic + 3, -jacobian(0,3));
					if(jacobian(0,4) != 0.0)tripletListA.emplace_back(ir + 0, ic + 4, -jacobian(0,4));
					if(jacobian(0,5) != 0.0)tripletListA.emplace_back(ir + 0, ic + 5, -jacobian(0,5));

					ic = nodes.size() * 6 + nodes[i].measurements[j].index_landmark * 3;
					if(jacobian(0,6) != 0.0)tripletListA.emplace_back(ir + 0, ic + 0, -jacobian(0,6));
					if(jacobian(0,7) != 0.0)tripletListA.emplace_back(ir + 0, ic + 1, -jacobian(0,7));
					if(jacobian(0,8) != 0.0)tripletListA.emplace_back(ir + 0, ic + 2, -jacobian(0,8));

					ic = 6 * i;
					if(jacobian(1,0) != 0.0)tripletListA.emplace_back(ir + 1, ic + 0, -jacobian(1,0));
					if(jacobian(1,1) != 0.0)tripletListA.emplace_back(ir + 1, ic + 1, -jacobian(1,1));
					if(jacobian(1,2) != 0.0)tripletListA.emplace_back(ir + 1, ic + 2, -jacobian(1,2));
					if(jacobian(1,3) != 0.0)tripletListA.emplace_back(ir + 1, ic + 3, -jacobian(1,3));
					if(jacobian(1,4) != 0.0)tripletListA.emplace_back(ir + 1, ic + 4, -jacobian(1,4));
					if(jacobian(1,5) != 0.0)tripletListA.emplace_back(ir + 1, ic + 5, -jacobian(1,5));

					ic = nodes.size() * 6 + nodes[i].measurements[j].index_landmark * 3;
					if(jacobian(1,6) != 0.0)tripletListA.emplace_back(ir + 1, ic + 0, -jacobian(1,6));
					if(jacobian(1,7) != 0.0)tripletListA.emplace_back(ir + 1, ic + 1, -jacobian(1,7));
					if(jacobian(1,8) != 0.0)tripletListA.emplace_back(ir + 1, ic + 2, -jacobian(1,8));

					ic = 6 * i;
					if(jacobian(2,0) != 0.0)tripletListA.emplace_back(ir + 2, ic + 0, -jacobian(2,0));
					if(jacobian(2,1) != 0.0)tripletListA.emplace_back(ir + 2, ic + 1, -jacobian(2,1));
					if(jacobian(2,2) != 0.0)tripletListA.emplace_back(ir + 2, ic + 2, -jacobian(2,2));
					if(jacobian(2,3) != 0.0)tripletListA.emplace_back(ir + 2, ic + 3, -jacobian(2,3));
					if(jacobian(2,4) != 0.0)tripletListA.emplace_back(ir + 2, ic + 4, -jacobian(2,4));
					if(jacobian(2,5) != 0.0)tripletListA.emplace_back(ir + 2, ic + 5, -jacobian(2,5));

					ic = nodes.size() * 6 + nodes[i].measurements[j].index_landmark * 3;
					if(jacobian(2,6) != 0.0)tripletListA.emplace_back(ir + 2, ic + 0, -jacobian(2,6));
					if(jacobian(2,7) != 0.0)tripletListA.emplace_back(ir + 2, ic + 1, -jacobian(2,7));
					if(jacobian(2,8) != 0.0)tripletListA.emplace_back(ir + 2, ic + 2, -jacobian(2,8));

					tripletListP.emplace_back(ir    , ir    , 1);
					tripletListP.emplace_back(ir + 1, ir + 1, 1);
					tripletListP.emplace_back(ir + 2, ir + 2, 1);

					tripletListB.emplace_back(ir    , 0,  delta_x);
					tripletListB.emplace_back(ir + 1, 0,  delta_y);
					tripletListB.emplace_back(ir + 2, 0,  delta_z);
				}
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);
			tripletListA.emplace_back(ir + 3 , 3, 1);
			tripletListA.emplace_back(ir + 4 , 4, 1);
			tripletListA.emplace_back(ir + 5 , 5, 1);

			tripletListP.emplace_back(ir     , ir,     10000000000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 10000000000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 10000000000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 10000000000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 10000000000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 10000000000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);

			Eigen::SparseMatrix<double> matA(tripletListB.size(), nodes.size() * 6 + landmarks.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(nodes.size() * 6 + landmarks.size() * 3, nodes.size() * 6 + landmarks.size() * 3);
			Eigen::SparseMatrix<double> AtPB(nodes.size() * 6 + landmarks.size() * 3, 1);

			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = AtP * matA;
			AtPB = AtP * matB;

			tripletListA.clear();
			tripletListP.clear();
			tripletListB.clear();

			Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(AtPA);
			Eigen::SparseMatrix<double> x = solver.solve(AtPB);

			std::vector<double> h_x;

			for (int k=0; k<x.outerSize(); ++k){
				for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it){
					std::cout << it.row() << " " << it.col() << " " << it.value() << std::endl;
					h_x.push_back(it.value());
				}
			}

			if(h_x.size() == nodes.size() * 6 + landmarks.size() * 3){
				std::cout << "optimization success" << std::endl;
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}

				int counter = 0;

				for(size_t i = 0 ; i < nodes.size(); i++){
					RodriguesPose pose = pose_rodrigues_from_affine_matrix(nodes[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.sx += h_x[counter++];
					pose.sy += h_x[counter++];
					pose.sz += h_x[counter++];
					nodes[i].pose = affine_matrix_from_pose_rodrigues(pose);
				}

				for(size_t i = 0 ; i < landmarks.size(); i++){
					landmarks[i].x() += h_x[counter++];
					landmarks[i].y() += h_x[counter++];
					landmarks[i].z() += h_x[counter++];
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'q':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < nodes.size(); i++){
				Eigen::Affine3d m = nodes[i].pose;
				QuaternionPose pose = pose_quaternion_from_affine_matrix(m);

				for(size_t j = 0 ; j < nodes[i].measurements.size(); j++){

					Eigen::Vector3d &p_t = landmarks[nodes[i].measurements[j].index_landmark];
					Eigen::Vector3d &p_s = nodes[i].measurements[j].value;

					double delta_x;
					double delta_y;
					double delta_z;
					point_to_point_source_to_landmark_quaternion_wc(delta_x, delta_y, delta_z,
							pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
							p_s.x(), p_s.y(), p_s.z(), p_t.x(), p_t.y(), p_t.z());

					Eigen::Matrix<double, 3, 10, Eigen::RowMajor> jacobian;
					point_to_point_source_to_landmark_quaternion_wc_jacobian(jacobian,
							pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
							p_s.x(), p_s.y(), p_s.z());


					int ir = tripletListB.size();
					int ic = 7 * i;

					if(jacobian(0,0) != 0.0)tripletListA.emplace_back(ir + 0, ic + 0, -jacobian(0,0));
					if(jacobian(0,1) != 0.0)tripletListA.emplace_back(ir + 0, ic + 1, -jacobian(0,1));
					if(jacobian(0,2) != 0.0)tripletListA.emplace_back(ir + 0, ic + 2, -jacobian(0,2));
					if(jacobian(0,3) != 0.0)tripletListA.emplace_back(ir + 0, ic + 3, -jacobian(0,3));
					if(jacobian(0,4) != 0.0)tripletListA.emplace_back(ir + 0, ic + 4, -jacobian(0,4));
					if(jacobian(0,5) != 0.0)tripletListA.emplace_back(ir + 0, ic + 5, -jacobian(0,5));
					if(jacobian(0,6) != 0.0)tripletListA.emplace_back(ir + 0, ic + 6, -jacobian(0,6));

					ic = nodes.size() * 7 + nodes[i].measurements[j].index_landmark * 3;
					if(jacobian(0,7) != 0.0)tripletListA.emplace_back(ir + 0, ic + 0, -jacobian(0,7));
					if(jacobian(0,8) != 0.0)tripletListA.emplace_back(ir + 0, ic + 1, -jacobian(0,8));
					if(jacobian(0,9) != 0.0)tripletListA.emplace_back(ir + 0, ic + 2, -jacobian(0,9));

					ic = 7 * i;
					if(jacobian(1,0) != 0.0)tripletListA.emplace_back(ir + 1, ic + 0, -jacobian(1,0));
					if(jacobian(1,1) != 0.0)tripletListA.emplace_back(ir + 1, ic + 1, -jacobian(1,1));
					if(jacobian(1,2) != 0.0)tripletListA.emplace_back(ir + 1, ic + 2, -jacobian(1,2));
					if(jacobian(1,3) != 0.0)tripletListA.emplace_back(ir + 1, ic + 3, -jacobian(1,3));
					if(jacobian(1,4) != 0.0)tripletListA.emplace_back(ir + 1, ic + 4, -jacobian(1,4));
					if(jacobian(1,5) != 0.0)tripletListA.emplace_back(ir + 1, ic + 5, -jacobian(1,5));
					if(jacobian(1,6) != 0.0)tripletListA.emplace_back(ir + 1, ic + 6, -jacobian(1,6));

					ic = nodes.size() * 7 + nodes[i].measurements[j].index_landmark * 3;
					if(jacobian(1,7) != 0.0)tripletListA.emplace_back(ir + 1, ic + 0, -jacobian(1,7));
					if(jacobian(1,8) != 0.0)tripletListA.emplace_back(ir + 1, ic + 1, -jacobian(1,8));
					if(jacobian(1,9) != 0.0)tripletListA.emplace_back(ir + 1, ic + 2, -jacobian(1,9));

					ic = 7 * i;
					if(jacobian(2,0) != 0.0)tripletListA.emplace_back(ir + 2, ic + 0, -jacobian(2,0));
					if(jacobian(2,1) != 0.0)tripletListA.emplace_back(ir + 2, ic + 1, -jacobian(2,1));
					if(jacobian(2,2) != 0.0)tripletListA.emplace_back(ir + 2, ic + 2, -jacobian(2,2));
					if(jacobian(2,3) != 0.0)tripletListA.emplace_back(ir + 2, ic + 3, -jacobian(2,3));
					if(jacobian(2,4) != 0.0)tripletListA.emplace_back(ir + 2, ic + 4, -jacobian(2,4));
					if(jacobian(2,5) != 0.0)tripletListA.emplace_back(ir + 2, ic + 5, -jacobian(2,5));
					if(jacobian(2,6) != 0.0)tripletListA.emplace_back(ir + 2, ic + 6, -jacobian(2,6));

					ic = nodes.size() * 7 + nodes[i].measurements[j].index_landmark * 3;
					if(jacobian(2,7) != 0.0)tripletListA.emplace_back(ir + 2, ic + 0, -jacobian(2,7));
					if(jacobian(2,8) != 0.0)tripletListA.emplace_back(ir + 2, ic + 1, -jacobian(2,8));
					if(jacobian(2,9) != 0.0)tripletListA.emplace_back(ir + 2, ic + 2, -jacobian(2,9));

					tripletListP.emplace_back(ir    , ir    , 1);
					tripletListP.emplace_back(ir + 1, ir + 1, 1);
					tripletListP.emplace_back(ir + 2, ir + 2, 1);

					tripletListB.emplace_back(ir    , 0,  delta_x);
					tripletListB.emplace_back(ir + 1, 0,  delta_y);
					tripletListB.emplace_back(ir + 2, 0,  delta_z);
				}
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);
			tripletListA.emplace_back(ir + 3 , 3, 1);
			tripletListA.emplace_back(ir + 4 , 4, 1);
			tripletListA.emplace_back(ir + 5 , 5, 1);
			tripletListA.emplace_back(ir + 6 , 6, 1);

			tripletListP.emplace_back(ir     , ir,     10000000000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 10000000000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 10000000000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 10000000000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 10000000000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 10000000000000);
			tripletListP.emplace_back(ir + 6 , ir + 6, 10000000000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);
			tripletListB.emplace_back(ir + 6 , 0, 0);

			for(size_t i = 0 ; i < nodes.size(); i++){
				int ic = i * 7;
				ir = tripletListB.size();
				QuaternionPose pose = pose_quaternion_from_affine_matrix(nodes[i].pose);

				double delta;
				quaternion_constraint(delta, pose.q0, pose.q1, pose.q2, pose.q3);

				Eigen::Matrix<double, 1, 4> jacobian;
				quaternion_constraint_jacobian(jacobian, pose.q0, pose.q1, pose.q2, pose.q3);

				tripletListA.emplace_back(ir, ic + 3 , -jacobian(0,0));
				tripletListA.emplace_back(ir, ic + 4 , -jacobian(0,1));
				tripletListA.emplace_back(ir, ic + 5 , -jacobian(0,2));
				tripletListA.emplace_back(ir, ic + 6 , -jacobian(0,3));

				tripletListP.emplace_back(ir, ir, 1000000.0);

				tripletListB.emplace_back(ir, 0, delta);
			}


			Eigen::SparseMatrix<double> matA(tripletListB.size(), nodes.size() * 7 + landmarks.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(nodes.size() * 7 + landmarks.size() * 3, nodes.size() * 7 + landmarks.size() * 3);
			Eigen::SparseMatrix<double> AtPB(nodes.size() * 7 + landmarks.size() * 3, 1);

			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = AtP * matA;
			AtPB = AtP * matB;

			tripletListA.clear();
			tripletListP.clear();
			tripletListB.clear();

			Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(AtPA);
			Eigen::SparseMatrix<double> x = solver.solve(AtPB);

			std::vector<double> h_x;

			for (int k=0; k<x.outerSize(); ++k){
				for (Eigen::SparseMatrix<double>::InnerIterator it(x,k); it; ++it){
					std::cout << it.row() << " " << it.col() << " " << it.value() << std::endl;
					h_x.push_back(it.value());
				}
			}

			if(h_x.size() == nodes.size() * 7 + landmarks.size() * 3){
				std::cout << "optimization success" << std::endl;
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}

				int counter = 0;

				for(size_t i = 0 ; i < nodes.size(); i++){
					QuaternionPose pose = pose_quaternion_from_affine_matrix(nodes[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.q0 += h_x[counter++];
					pose.q1 += h_x[counter++];
					pose.q2 += h_x[counter++];
					pose.q3 += h_x[counter++];
					nodes[i].pose = affine_matrix_from_pose_quaternion(pose);
				}
				for(size_t i = 0 ; i < landmarks.size(); i++){
					landmarks[i].x() += h_x[counter++];
					landmarks[i].y() += h_x[counter++];
					landmarks[i].z() += h_x[counter++];
				}
			}else{
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
	std::cout << "r: optimize (Rodrigues)" << std::endl;
	std::cout << "q: optimize (Quaternion)" << std::endl;
}



