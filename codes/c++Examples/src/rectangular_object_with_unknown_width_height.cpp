#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "../../rectangular_object_with_unknown_width_height_tait_bryan_wc_jacobian.h"
#include "../../rectangular_object_with_unknown_width_height_rodrigues_wc_jacobian.h"
#include "../../rectangular_object_with_unknown_width_height_quaternion_wc_jacobian.h"
#include "../../quaternion_constraint_jacobian.h"

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -10.0;
float translate_x, translate_y = 0.0;

std::vector<std::vector<Eigen::Affine3d>> bundle_of_rays;
struct Rectangle{
	std::vector<Eigen::Vector3d> corners_local;
	double scale_x;
	double scale_y;
	Eigen::Affine3d pose;
};
Rectangle rectangle;

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

int main(int argc, char *argv[]){

	std::vector<Eigen::Affine3d> br;
	for(size_t i = 0; i < 25; i++){
		TaitBryanPose pose;

		pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + 10;
		pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + -5;
		pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + -2;

		pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;
		pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;
		pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;

		br.push_back(affine_matrix_from_pose_tait_bryan(pose));
	}
	bundle_of_rays.push_back(br);
	br.clear();
	for(size_t i = 0; i < 25; i++){
		TaitBryanPose pose;

		pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + 10;
		pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 +  5;
		pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + -2;

		pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;
		pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;
		pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;

		br.push_back(affine_matrix_from_pose_tait_bryan(pose));
	}
	bundle_of_rays.push_back(br);
	br.clear();
	for(size_t i = 0; i < 25; i++){
		TaitBryanPose pose;

		pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + 10;
		pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 +  5;
		pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 +  5;

		pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;
		pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;
		pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;

		br.push_back(affine_matrix_from_pose_tait_bryan(pose));
	}
	bundle_of_rays.push_back(br);
	br.clear();
	for(size_t i = 0; i < 25; i++){
		TaitBryanPose pose;

		pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + 10;
		pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 + -5;
		pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.01 +  5;

		pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;
		pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;
		pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 1.5;

		br.push_back(affine_matrix_from_pose_tait_bryan(pose));
	}
	bundle_of_rays.push_back(br);
	br.clear();


	rectangle.scale_x = 1.0;
	rectangle.scale_y = 1.0;
	rectangle.pose = Eigen::Affine3d::Identity();
	rectangle.corners_local.emplace_back(-1,-1, 0);
	rectangle.corners_local.emplace_back( 1,-1, 0);
	rectangle.corners_local.emplace_back( 1, 1, 0);
	rectangle.corners_local.emplace_back(-1, 1, 0);



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
	glutCreateWindow("rectangular_object_with_unknown_width_height");
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

	glColor3f(0,1,0);
	glBegin(GL_LINES);
	for(size_t i = 0; i < bundle_of_rays.size(); i++){
		for(size_t j = 0; j < bundle_of_rays[i].size(); j++){
			Eigen::Vector3d z_begin(0, 0,-100);
			Eigen::Vector3d z_end(0, 0, 100);

			Eigen::Vector3d z_begin_t = bundle_of_rays[i][j] * z_begin;
			Eigen::Vector3d z_end_t = bundle_of_rays[i][j] * z_end;

			glVertex3f(z_begin_t.x(), z_begin_t.y(), z_begin_t.z());
			glVertex3f(z_end_t.x(), z_end_t.y(), z_end_t.z());
		}
	}
	glEnd();

	std::vector<Eigen::Vector3d> corners_global;
	for(size_t i = 0; i < rectangle.corners_local.size(); i++){
		Eigen::Vector3d v(rectangle.corners_local[i].x() * rectangle.scale_x, rectangle.corners_local[i].y() * rectangle.scale_y, rectangle.corners_local[i].z());
		corners_global.push_back(rectangle.pose * v);
	}

	glLineWidth(5);
	glColor3f(1,0,0);
	glBegin(GL_LINE_STRIP);
		for(size_t i = 0 ; i < corners_global.size(); i++){
			glVertex3f(corners_global[i].x(), corners_global[i].y(), corners_global[i].z());
		}
		glVertex3f(corners_global[0].x(), corners_global[0].y(), corners_global[0].z());
	glEnd();
	glLineWidth(1);
	glutSwapBuffers();
}


void keyboard(unsigned char key, int /*x*/, int /*y*/) {
	switch (key) {
		case (27): {
			glutDestroyWindow(glutGetWindow());
			return;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(rectangle.pose);

			for(size_t i = 0; i < bundle_of_rays.size(); i++){
				for(size_t j = 0; j < bundle_of_rays[i].size(); j++){
					Eigen::Vector3d vx(bundle_of_rays[i][j](0,0), bundle_of_rays[i][j](1,0), bundle_of_rays[i][j](2,0));
					Eigen::Vector3d vy(bundle_of_rays[i][j](0,1), bundle_of_rays[i][j](1,1), bundle_of_rays[i][j](2,1));

					Eigen::Matrix<double, 2, 1> delta;
					rectangular_object_with_unknown_width_height_tait_bryan_wc(delta,
							pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
							rectangle.corners_local[i].x(), rectangle.corners_local[i].y(), rectangle.corners_local[i].z(),
							bundle_of_rays[i][j](0,3), bundle_of_rays[i][j](1,3), bundle_of_rays[i][j](2,3),
							vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z(),
							rectangle.scale_x, rectangle.scale_y);

					Eigen::Matrix<double, 2, 8> delta_jacobian;
					rectangular_object_with_unknown_width_height_tait_bryan_wc_jacobian(delta_jacobian,
							pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
							rectangle.corners_local[i].x(), rectangle.corners_local[i].y(), rectangle.corners_local[i].z(),
							bundle_of_rays[i][j](0,3), bundle_of_rays[i][j](1,3), bundle_of_rays[i][j](2,3),
							vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z(),
							rectangle.scale_x, rectangle.scale_y);

					int ir = tripletListB.size();
					for(size_t k = 0 ; k < 8; k++){
						tripletListA.emplace_back(ir, k, -delta_jacobian(0,k));
					}
					for(size_t k = 0 ; k < 8; k++){
						tripletListA.emplace_back(ir + 1, k, -delta_jacobian(1,k));
					}

					tripletListP.emplace_back(ir    , ir    ,  1);
					tripletListP.emplace_back(ir + 1, ir + 1,  1);

					tripletListB.emplace_back(ir    , 0,  delta(0,0));
					tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
				}
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), 8);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(8, 8);
			Eigen::SparseMatrix<double> AtPB(8, 1);

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
					std::cout << it.value() << std::endl;
				}
			}

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			if(h_x.size() == 8){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;
				TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(rectangle.pose);
				pose.px += h_x[counter++];
				pose.py += h_x[counter++];
				pose.pz += h_x[counter++];
				pose.om += h_x[counter++];
				pose.fi += h_x[counter++];
				pose.ka += h_x[counter++];
				rectangle.pose = affine_matrix_from_pose_tait_bryan(pose);
				rectangle.scale_x += h_x[counter++];
				rectangle.scale_y += h_x[counter++];
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'r':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			TaitBryanPose pose_rand;
			pose_rand.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			pose_rand.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			pose_rand.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			pose_rand.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			pose_rand.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			pose_rand.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			rectangle.pose = rectangle.pose * affine_matrix_from_pose_tait_bryan(pose_rand);


			RodriguesPose pose = pose_rodrigues_from_affine_matrix(rectangle.pose);

			for(size_t i = 0; i < bundle_of_rays.size(); i++){
				for(size_t j = 0; j < bundle_of_rays[i].size(); j++){
					Eigen::Vector3d vx(bundle_of_rays[i][j](0,0), bundle_of_rays[i][j](1,0), bundle_of_rays[i][j](2,0));
					Eigen::Vector3d vy(bundle_of_rays[i][j](0,1), bundle_of_rays[i][j](1,1), bundle_of_rays[i][j](2,1));

					Eigen::Matrix<double, 2, 1> delta;
					rectangular_object_with_unknown_width_height_rodrigues_wc(delta,
							pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
							rectangle.corners_local[i].x(), rectangle.corners_local[i].y(), rectangle.corners_local[i].z(),
							bundle_of_rays[i][j](0,3), bundle_of_rays[i][j](1,3), bundle_of_rays[i][j](2,3),
							vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z(),
							rectangle.scale_x, rectangle.scale_y);

					Eigen::Matrix<double, 2, 8> delta_jacobian;
					rectangular_object_with_unknown_width_height_rodrigues_wc_jacobian(delta_jacobian,
							pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
							rectangle.corners_local[i].x(), rectangle.corners_local[i].y(), rectangle.corners_local[i].z(),
							bundle_of_rays[i][j](0,3), bundle_of_rays[i][j](1,3), bundle_of_rays[i][j](2,3),
							vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z(),
							rectangle.scale_x, rectangle.scale_y);

					int ir = tripletListB.size();
					for(size_t k = 0 ; k < 8; k++){
						tripletListA.emplace_back(ir, k, -delta_jacobian(0,k));
					}
					for(size_t k = 0 ; k < 8; k++){
						tripletListA.emplace_back(ir + 1, k, -delta_jacobian(1,k));
					}

					tripletListP.emplace_back(ir    , ir    ,  1);
					tripletListP.emplace_back(ir + 1, ir + 1,  1);

					tripletListB.emplace_back(ir    , 0,  delta(0,0));
					tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
				}
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), 8);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(8, 8);
			Eigen::SparseMatrix<double> AtPB(8, 1);

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
					std::cout << it.value() << std::endl;
				}
			}

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			if(h_x.size() == 8){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;
				RodriguesPose pose = pose_rodrigues_from_affine_matrix(rectangle.pose);
				pose.px += h_x[counter++];
				pose.py += h_x[counter++];
				pose.pz += h_x[counter++];
				pose.sx += h_x[counter++];
				pose.sy += h_x[counter++];
				pose.sz += h_x[counter++];
				rectangle.pose = affine_matrix_from_pose_rodrigues(pose);
				rectangle.scale_x += h_x[counter++];
				rectangle.scale_y += h_x[counter++];
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'q':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			TaitBryanPose pose_rand;
			pose_rand.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			pose_rand.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			pose_rand.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			pose_rand.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			pose_rand.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			pose_rand.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			rectangle.pose = rectangle.pose * affine_matrix_from_pose_tait_bryan(pose_rand);


			QuaternionPose pose = pose_quaternion_from_affine_matrix(rectangle.pose);

			for(size_t i = 0; i < bundle_of_rays.size(); i++){
				for(size_t j = 0; j < bundle_of_rays[i].size(); j++){
					Eigen::Vector3d vx(bundle_of_rays[i][j](0,0), bundle_of_rays[i][j](1,0), bundle_of_rays[i][j](2,0));
					Eigen::Vector3d vy(bundle_of_rays[i][j](0,1), bundle_of_rays[i][j](1,1), bundle_of_rays[i][j](2,1));

					Eigen::Matrix<double, 2, 1> delta;
					rectangular_object_with_unknown_width_height_quaternion_wc(delta,
							pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
							rectangle.corners_local[i].x(), rectangle.corners_local[i].y(), rectangle.corners_local[i].z(),
							bundle_of_rays[i][j](0,3), bundle_of_rays[i][j](1,3), bundle_of_rays[i][j](2,3),
							vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z(),
							rectangle.scale_x, rectangle.scale_y);

					Eigen::Matrix<double, 2, 9> delta_jacobian;
					rectangular_object_with_unknown_width_height_quaternion_wc_jacobian(delta_jacobian,
							pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
							rectangle.corners_local[i].x(), rectangle.corners_local[i].y(), rectangle.corners_local[i].z(),
							bundle_of_rays[i][j](0,3), bundle_of_rays[i][j](1,3), bundle_of_rays[i][j](2,3),
							vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z(),
							rectangle.scale_x, rectangle.scale_y);

					int ir = tripletListB.size();
					for(size_t k = 0 ; k < 9; k++){
						tripletListA.emplace_back(ir, k, -delta_jacobian(0,k));
					}
					for(size_t k = 0 ; k < 9; k++){
						tripletListA.emplace_back(ir + 1, k, -delta_jacobian(1,k));
					}

					tripletListP.emplace_back(ir    , ir    ,  1);
					tripletListP.emplace_back(ir + 1, ir + 1,  1);

					tripletListB.emplace_back(ir    , 0,  delta(0,0));
					tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
				}
			}


			int ir = tripletListB.size();

			double delta;
			quaternion_constraint(delta, pose.q0, pose.q1, pose.q2, pose.q3);

			Eigen::Matrix<double, 1, 4> jacobian;
			quaternion_constraint_jacobian(jacobian, pose.q0, pose.q1, pose.q2, pose.q3);

			tripletListA.emplace_back(ir, 3 , -jacobian(0,0));
			tripletListA.emplace_back(ir, 4 , -jacobian(0,1));
			tripletListA.emplace_back(ir, 5 , -jacobian(0,2));
			tripletListA.emplace_back(ir, 6 , -jacobian(0,3));

			tripletListP.emplace_back(ir, ir, 1000000.0);

			tripletListB.emplace_back(ir, 0, delta);

			Eigen::SparseMatrix<double> matA(tripletListB.size(), 9);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(9, 9);
			Eigen::SparseMatrix<double> AtPB(9, 1);

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
					std::cout << it.value() << std::endl;
				}
			}

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			if(h_x.size() == 9){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;
				QuaternionPose pose = pose_quaternion_from_affine_matrix(rectangle.pose);
				pose.px += h_x[counter++];
				pose.py += h_x[counter++];
				pose.pz += h_x[counter++];
				pose.q0 += h_x[counter++];
				pose.q1 += h_x[counter++];
				pose.q2 += h_x[counter++];
				pose.q3 += h_x[counter++];
				rectangle.pose = affine_matrix_from_pose_quaternion(pose);
				rectangle.scale_x += h_x[counter++];
				rectangle.scale_y += h_x[counter++];
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'n':{
			for(size_t i = 0; i < bundle_of_rays.size(); i++){
				for(size_t j = 0; j < bundle_of_rays[i].size(); j++){
					TaitBryanPose pose;

					pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.05;
					pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.05;
					pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 0.05;

					pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;
					pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;
					pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;

					bundle_of_rays[i][j] = bundle_of_rays[i][j] * affine_matrix_from_pose_tait_bryan(pose);
				}
			}
			TaitBryanPose pose;

			pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 5;
			pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 5;
			pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 5;

			pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;
			pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;
			pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 2;

			rectangle.pose = rectangle.pose * affine_matrix_from_pose_tait_bryan(pose);

			rectangle.scale_x += ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 5;
			rectangle.scale_y += ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 5;
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
	std::cout << "n: modify rays" << std::endl;
	std::cout << "t: optimize (Tait-Bryan)" << std::endl;
	std::cout << "r: optimize (Rodrigues)" << std::endl;
	std::cout << "q: optimize (Quaternion)" << std::endl;
}
