#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"
#include "../../perspective_camera_plucker_line_tait_bryan_wc_jacobian.h"
#include "../../perspective_camera_tait_bryan_wc_jacobian.h"
#include "../../perspective_camera_plucker_line_rodrigues_wc_jacobian.h"
#include "../../perspective_camera_plucker_line_quaternion_wc_jacobian.h"
#include "../../quaternion_constraint_jacobian.h"


struct LinePoint{
	double u;
	double v;
	int index_to_line;
};

struct Camera{
	Eigen::Affine3d pose;
	std::vector<std::pair<LinePoint, LinePoint>> uv;
};

typedef Eigen::Matrix<double, 6, 1> PLine;

std::vector<Camera> cameras;
PerspectiveCameraParams cam_params;
std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> lines;

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = -215, rotate_y = 273;
float translate_z = -68.0;
float translate_x = 4, translate_y = 30.0;

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();
PLine get_plucker_line(const Eigen::Vector3d &from, const Eigen::Vector3d &to);

int main(int argc, char *argv[]){

	if (false == initGL(&argc, argv)) {
		return 4;
	}

	cam_params.fx = 2000;
	cam_params.fy = 1000;
	cam_params.cx = 1000;
	cam_params.cy = 500;

	for(int i = -50 ; i < 50; i+=20){
		Camera c;
		c.pose = Eigen::Affine3d::Identity();
		c.pose(0,3) = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1.0;
		c.pose(1,3) = i;
		c.pose(2,3) = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1.0;
		cameras.push_back(c);
	}

	for(float x = 0; x <= 100; x += 5){
		Eigen::Vector3d a(x, -50.0, 400);
		Eigen::Vector3d b(x,  50.0, 400);
		lines.emplace_back(a,b);
	}
	for(float y = -50; y <= 50; y += 5){
		Eigen::Vector3d a(0, y, 400);
		Eigen::Vector3d b(100.0, y, 400);
		lines.emplace_back(a,b);
	}

	for(float x = -100; x <= 0; x += 5){
		Eigen::Vector3d a(x, -50.0, 440);
		Eigen::Vector3d b(x,  50.0, 440);
		lines.emplace_back(a,b);
	}
	for(float y = -50; y <= 50; y += 5){
		Eigen::Vector3d a(-100, y, 440);
		Eigen::Vector3d b(0.0, y, 440);
		lines.emplace_back(a,b);
	}

	for(size_t i = 0; i < cameras.size(); i++){
		for(size_t j = 0 ; j < lines.size(); j++){
			std::pair<LinePoint, LinePoint> uv;
			LinePoint kp;
			TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose);

			projection_perspective_camera_tait_bryan_wc(kp.u, kp.v, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka, lines[j].first.x(), lines[j].first.y(), lines[j].first.z());

			uv.first.u = kp.u;
			uv.first.v = kp.v;
			uv.first.index_to_line = j;

			projection_perspective_camera_tait_bryan_wc(kp.u, kp.v, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka, lines[j].second.x(), lines[j].second.y(), lines[j].second.z());

			uv.second.u = kp.u;
			uv.second.v = kp.v;
			uv.second.index_to_line = j;

			cameras[i].uv.push_back(uv);
		}
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
	glutCreateWindow("perspective_camera_external_orientation_plucker_line");
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

	for(size_t i = 0 ; i < cameras.size(); i++){
		Eigen::Affine3d m = cameras[i].pose;

		glBegin(GL_LINES);
			glColor3f(1.0f, 0.0f, 0.0f);
			glVertex3f(m(0,3), m(1,3), m(2,3));
			glVertex3f(m(0,3) + m(0,0), m(1,3) + m(1,0), m(2,3) + m(2,0));

			glColor3f(0.0f, 1.0f, 0.0f);
			glVertex3f(m(0,3), m(1,3), m(2,3));
			glVertex3f(m(0,3) + m(0,1), m(1,3) + m(1,1), m(2,3) + m(2,1));

			glColor3f(0.0f, 0.0f, 1.0f);
			glVertex3f(m(0,3), m(1,3), m(2,3));
			glVertex3f(m(0,3) + m(0,2), m(1,3) + m(1,2), m(2,3) + m(2,2));
		glEnd();
	}

	glColor3f(0,1,0);
	glBegin(GL_LINES);
		for(size_t i = 0; i < lines.size(); i++){
			glVertex3f(lines[i].first.x(), lines[i].first.y(), lines[i].first.z());
			glVertex3f(lines[i].second.x(), lines[i].second.y(), lines[i].second.z());
		}
	glEnd();

	glColor3f(0.8,0.8,0.8);
	glBegin(GL_LINE_STRIP);
	for(size_t i = 0 ; i < cameras.size(); i++){
		Eigen::Affine3d m = cameras[i].pose;
		glVertex3f(m(0,3), m(1,3), m(2,3));
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
		case 'c':{
			for(size_t i = 0; i < cameras.size(); i++){
				TaitBryanPose pose;
				pose.px = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1;
				pose.py = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1;
				pose.pz = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 1;
				pose.om = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.1;
				pose.fi = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.1;
				pose.ka = (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.1;

				Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(pose);
				cameras[i].pose = cameras[i].pose * m;
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras[i].uv.size(); j++){

					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose);

					PLine pl = get_plucker_line(lines[cameras[i].uv[j].first.index_to_line].first, lines[cameras[i].uv[j].first.index_to_line].second);
					Eigen::Matrix<double, 2, 1> delta;
					observation_equation_perspective_camera_plucker_line_tait_bryan_wc(delta, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
							pl(0,0), pl(1,0), pl(2,0), pl(3,0), pl(4,0), pl(5,0),
							cameras[i].uv[j].first.u, cameras[i].uv[j].first.v,
							cameras[i].uv[j].second.u, cameras[i].uv[j].second.v);

					Eigen::Matrix<double, 2, 6, Eigen::RowMajor> jacobian;
					observation_equation_perspective_camera_plucker_line_tait_bryan_wc_jacobian(jacobian, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
							pl(0,0), pl(1,0), pl(2,0), pl(3,0), pl(4,0), pl(5,0),
							cameras[i].uv[j].first.u, cameras[i].uv[j].first.v,
							cameras[i].uv[j].second.u, cameras[i].uv[j].second.v);

					int ir = tripletListB.size();
					int ic_camera = i * 6;
					tripletListA.emplace_back(ir     , ic_camera    , -jacobian(0,0));
					tripletListA.emplace_back(ir     , ic_camera + 1, -jacobian(0,1));
					tripletListA.emplace_back(ir     , ic_camera + 2, -jacobian(0,2));
					tripletListA.emplace_back(ir     , ic_camera + 3, -jacobian(0,3));
					tripletListA.emplace_back(ir     , ic_camera + 4, -jacobian(0,4));
					tripletListA.emplace_back(ir     , ic_camera + 5, -jacobian(0,5));

					tripletListA.emplace_back(ir + 1 , ic_camera    , -jacobian(1,0));
					tripletListA.emplace_back(ir + 1 , ic_camera + 1, -jacobian(1,1));
					tripletListA.emplace_back(ir + 1 , ic_camera + 2, -jacobian(1,2));
					tripletListA.emplace_back(ir + 1 , ic_camera + 3, -jacobian(1,3));
					tripletListA.emplace_back(ir + 1 , ic_camera + 4, -jacobian(1,4));
					tripletListA.emplace_back(ir + 1 , ic_camera + 5, -jacobian(1,5));

					tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
					tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));

					tripletListB.emplace_back(ir    , 0,  delta(0,0));
					tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
				}
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 6, cameras.size() * 6);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 6, 1);

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
				}
			}

			if(h_x.size() == cameras.size() * 6){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < cameras.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_tait_bryan(pose);
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}

			break;
		}
		case 'r':{
			for(size_t i = 0; i < cameras.size(); i++){
				TaitBryanPose posetb = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
				posetb.px += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.py += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.pz += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.om += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.fi += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				posetb.ka += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00001;
				cameras[i].pose = affine_matrix_from_pose_tait_bryan(posetb);
			}

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras[i].uv.size(); j++){

					RodriguesPose pose = pose_rodrigues_from_affine_matrix(cameras[i].pose);
					Eigen::Matrix<double, 2, 1> delta;
					PLine pl = get_plucker_line(lines[cameras[i].uv[j].first.index_to_line].first, lines[cameras[i].uv[j].first.index_to_line].second);
					observation_equation_perspective_camera_plucker_line_rodrigues_wc(delta, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy, pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
							pl(0,0), pl(1,0), pl(2,0), pl(3,0), pl(4,0), pl(5,0),
							cameras[i].uv[j].first.u, cameras[i].uv[j].first.v,
							cameras[i].uv[j].second.u, cameras[i].uv[j].second.v);

					Eigen::Matrix<double, 2, 6, Eigen::RowMajor> jacobian;
					observation_equation_perspective_camera_plucker_line_rodrigues_wc_jacobian(jacobian, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy, pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
							pl(0,0), pl(1,0), pl(2,0), pl(3,0), pl(4,0), pl(5,0),
							cameras[i].uv[j].first.u, cameras[i].uv[j].first.v,
							cameras[i].uv[j].second.u, cameras[i].uv[j].second.v);

					int ir = tripletListB.size();
					int ic_camera = i * 6;
					tripletListA.emplace_back(ir     , ic_camera    , -jacobian(0,0));
					tripletListA.emplace_back(ir     , ic_camera + 1, -jacobian(0,1));
					tripletListA.emplace_back(ir     , ic_camera + 2, -jacobian(0,2));
					tripletListA.emplace_back(ir     , ic_camera + 3, -jacobian(0,3));
					tripletListA.emplace_back(ir     , ic_camera + 4, -jacobian(0,4));
					tripletListA.emplace_back(ir     , ic_camera + 5, -jacobian(0,5));

					tripletListA.emplace_back(ir + 1 , ic_camera    , -jacobian(1,0));
					tripletListA.emplace_back(ir + 1 , ic_camera + 1, -jacobian(1,1));
					tripletListA.emplace_back(ir + 1 , ic_camera + 2, -jacobian(1,2));
					tripletListA.emplace_back(ir + 1 , ic_camera + 3, -jacobian(1,3));
					tripletListA.emplace_back(ir + 1 , ic_camera + 4, -jacobian(1,4));
					tripletListA.emplace_back(ir + 1 , ic_camera + 5, -jacobian(1,5));

					tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
					tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));

					tripletListB.emplace_back(ir    , 0,  delta(0,0));
					tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
				}
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 6, cameras.size() * 6);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 6, 1);

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
				}
			}

			if(h_x.size() == cameras.size() * 6){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < cameras.size(); i++){
					RodriguesPose pose = pose_rodrigues_from_affine_matrix(cameras[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.sx += h_x[counter++];
					pose.sy += h_x[counter++];
					pose.sz += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_rodrigues(pose);
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

			for(size_t i = 0; i < cameras.size(); i++){
				for(size_t j = 0; j < cameras[i].uv.size(); j++){

					QuaternionPose pose = pose_quaternion_from_affine_matrix(cameras[i].pose);
					Eigen::Matrix<double, 2, 1> delta;
					PLine pl = get_plucker_line(lines[cameras[i].uv[j].first.index_to_line].first, lines[cameras[i].uv[j].first.index_to_line].second);
					observation_equation_perspective_camera_plucker_line_quaternion_wc(delta, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy, pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
							pl(0,0), pl(1,0), pl(2,0), pl(3,0), pl(4,0), pl(5,0),
							cameras[i].uv[j].first.u, cameras[i].uv[j].first.v,
							cameras[i].uv[j].second.u, cameras[i].uv[j].second.v);

					Eigen::Matrix<double, 2, 7, Eigen::RowMajor> jacobian;
					observation_equation_perspective_camera_plucker_line_quaternion_wc_jacobian(jacobian, cam_params.fx, cam_params.fy, cam_params.cx, cam_params.cy, pose.px, pose.py, pose.pz, pose.q0, pose.q1, pose.q2, pose.q3,
							pl(0,0), pl(1,0), pl(2,0), pl(3,0), pl(4,0), pl(5,0),
							cameras[i].uv[j].first.u, cameras[i].uv[j].first.v,
							cameras[i].uv[j].second.u, cameras[i].uv[j].second.v);

					int ir = tripletListB.size();
					int ic_camera = i * 7;
					tripletListA.emplace_back(ir     , ic_camera    , -jacobian(0,0));
					tripletListA.emplace_back(ir     , ic_camera + 1, -jacobian(0,1));
					tripletListA.emplace_back(ir     , ic_camera + 2, -jacobian(0,2));
					tripletListA.emplace_back(ir     , ic_camera + 3, -jacobian(0,3));
					tripletListA.emplace_back(ir     , ic_camera + 4, -jacobian(0,4));
					tripletListA.emplace_back(ir     , ic_camera + 5, -jacobian(0,5));
					tripletListA.emplace_back(ir     , ic_camera + 6, -jacobian(0,6));

					tripletListA.emplace_back(ir + 1 , ic_camera    , -jacobian(1,0));
					tripletListA.emplace_back(ir + 1 , ic_camera + 1, -jacobian(1,1));
					tripletListA.emplace_back(ir + 1 , ic_camera + 2, -jacobian(1,2));
					tripletListA.emplace_back(ir + 1 , ic_camera + 3, -jacobian(1,3));
					tripletListA.emplace_back(ir + 1 , ic_camera + 4, -jacobian(1,4));
					tripletListA.emplace_back(ir + 1 , ic_camera + 5, -jacobian(1,5));
					tripletListA.emplace_back(ir + 1 , ic_camera + 6, -jacobian(1,6));

					tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
					tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));

					tripletListB.emplace_back(ir    , 0,  delta(0,0));
					tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
				}
			}

			for(size_t i = 0 ; i < cameras.size(); i++){
				int ic = i * 7;
				int ir = tripletListB.size();
				QuaternionPose pose = pose_quaternion_from_affine_matrix(cameras[i].pose);

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

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 7);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 7, cameras.size() * 7);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 7, 1);

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
				}
			}

			if(h_x.size() == cameras.size() * 7){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < cameras.size(); i++){
					QuaternionPose pose = pose_quaternion_from_affine_matrix(cameras[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.q0 += h_x[counter++];
					pose.q1 += h_x[counter++];
					pose.q2 += h_x[counter++];
					pose.q3 += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_quaternion(pose);
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
	std::cout << "c: add noise to cameras" << std::endl;
	std::cout << "t: optimize cameras external orientation (Tait-Bryan)" << std::endl;
	std::cout << "r: optimize cameras external orientation (Rodrigues)" << std::endl;
	std::cout << "q: optimize cameras external orientation (Quaternion)" << std::endl;
}

PLine get_plucker_line(const Eigen::Vector3d &from, const Eigen::Vector3d &to)
{
	PLine plucker_line;

	Eigen::Vector3d direction = (to-from);
	direction/=direction.norm();
	Eigen::Vector3d moment = from.cross(direction);

	plucker_line(0,0) = moment.x();
	plucker_line(1,0) = moment.y();
	plucker_line(2,0) = moment.z();
	plucker_line(3,0) = direction.x();
	plucker_line(4,0) = direction.y();
	plucker_line(5,0) = direction.z();

	return plucker_line;
}


