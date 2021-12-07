#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"
#include "example_func_xy_jacobian.h"
#include "constraints_jacobian.h"
#include "relative_pose_tait_bryan_wc_jacobian.h"
#include "point_to_point_source_to_target_tait_bryan_wc_jacobian.h"
#include "smoothness_tait_bryan_wc_jacobian.h"

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -50.0;
float translate_x, translate_y = 0.0;

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

std::vector<Eigen::Vector3d> reference_points;

std::vector<Eigen::Affine3d> vertices;
std::vector<Eigen::Affine3d> vertices_odo;
int number_rows = 0;
int number_cols = 0;
std::vector<std::pair<int, int>> odo_edges;

struct TripletIndexes{
	int index_before;
	int index_curr;
	int index_after;
};

std::vector<TripletIndexes> smoothness_indexes;

bool find_nearest_neighbour(Eigen::Vector3d &p_t, Eigen::Affine3d vertex, const std::vector<Eigen::Vector3d>& reference_points){
	double min_dist_xy = 1000000.0;
	double search_radious = 1.0;
	bool found = false;

	for(size_t i = 0 ; i < reference_points.size(); i++){
		float dist = sqrt ( (vertex(0,3) - reference_points[i].x()) * (vertex(0,3) - reference_points[i].x()) + (vertex(1,3) - reference_points[i].y()) * (vertex(1,3) - reference_points[i].y()) );

		if( dist < search_radious  ) {
			if( dist < min_dist_xy){
				min_dist_xy = dist;
				p_t = reference_points[i];
				p_t.x() = vertex(0,3) + ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
				p_t.y() = vertex(1,3) + ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
				found = true;
			}
		}
	}
	return found;
}

int main(int argc, char *argv[]){

	for(size_t i = 0; i < 10000; i++){
		double x = ((float(rand()%1000000)/1000000.0f - 0.5) * 2.0) * 20.0;
		double y = ((float(rand()%1000000)/1000000.0f - 0.5) * 2.0) * 20.0;
		double z;
		example_func_xy(z, x, y);

		if( (x-4)*(x-4) + (y-5)*(y-5) > 49){
			reference_points.emplace_back(x,y,z*10 + ((float(rand()%1000000)/1000000.0f - 0.5) * 2.0) * 0.1);
		}
	}

	for(int x = -20; x <= 20; x += 1){
		for(int y = -20; y <= 20; y += 1){
			if(number_rows == 0){
				number_cols ++;
			}
			Eigen::Affine3d pose = Eigen::Affine3d::Identity();
			pose(0,3) = x;
			pose(1,3) = y;
			pose(2,3) = -10;
			vertices.push_back(pose);
		}
		number_rows ++;
	}
	vertices_odo = vertices;

	for(int row = 0; row < number_rows; row++){
		for(int col = 0; col < number_cols; col++){
			int index = row + col * number_rows;
			if(row + 1 < number_rows){
				int index_next_row = (row + 1) + col * number_rows;
				odo_edges.emplace_back(index, index_next_row);
			}

			if(col + 1 < number_cols){
				int index_next_col = (row) + (col + 1) * number_rows;
				odo_edges.emplace_back(index, index_next_col);
			}

			if(row + 1 < number_rows && row - 1 >= 0){
				int index_prev_row = (row - 1) + col * number_rows;
				int index_next_row = (row + 1) + col * number_rows;
				TripletIndexes tr;
				tr.index_curr = index;
				tr.index_before = index_prev_row;
				tr.index_after = index_next_row;
				smoothness_indexes.push_back(tr);
			}
			if(col + 1 < number_cols && col - 1 >= 0){
				int index_next_col = (row) + (col + 1) * number_rows;
				int index_prev_col = (row) + (col - 1) * number_rows;
				TripletIndexes tr;
				tr.index_curr = index;
				tr.index_before = index_prev_col;
				tr.index_after = index_next_col;
				smoothness_indexes.push_back(tr);
			}
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
	glutCreateWindow("surface_reconstruction_from_lidar_point_cloud");
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

	glPointSize(3);
	glColor3f(1,0,0);
	glBegin(GL_POINTS);
	for(size_t i = 0; i < reference_points.size(); i++){
		glVertex3f(reference_points[i].x(), reference_points[i].y(), reference_points[i].z());
	}
	glEnd();
	glPointSize(1);

	glPointSize(2);
	glColor3f(0,0,0);
	glBegin(GL_POINTS);
	for(int row = 0; row < number_rows; row++){
		for(int col = 0; col < number_cols; col++){
			int index = row + col * number_rows;
			glVertex3f(vertices[index](0,3), vertices[index](1,3), vertices[index](2,3));
		}
	}
	glEnd();
	glPointSize(1);

	glBegin(GL_LINES);
	for(size_t i = 0; i < odo_edges.size(); i++){
		glVertex3f(vertices[odo_edges[i].first](0,3), vertices[odo_edges[i].first](1,3), vertices[odo_edges[i].first](2,3));
		glVertex3f(vertices[odo_edges[i].second](0,3), vertices[odo_edges[i].second](1,3), vertices[odo_edges[i].second](2,3));
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
			for(size_t i = 0 ; i < vertices.size(); i++){
				TaitBryanPose pose;
				pose.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;

				vertices[i] = vertices[i] * affine_matrix_from_pose_tait_bryan(pose);
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<TaitBryanPose> poses;
			std::vector<TaitBryanPose> poses_odo;

			for(size_t i = 0 ; i < vertices.size(); i++){
				poses.push_back(pose_tait_bryan_from_affine_matrix(vertices[i]));
			}
			for(size_t i = 0 ; i < vertices_odo.size(); i++){
				poses_odo.push_back(pose_tait_bryan_from_affine_matrix(vertices_odo[i]));
			}

			for(size_t i = 0 ; i < odo_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_odo;
				relative_pose_tait_bryan_wc_case1(relative_pose_measurement_odo,
						poses_odo[odo_edges[i].first].px,
						poses_odo[odo_edges[i].first].py,
						poses_odo[odo_edges[i].first].pz,
						poses_odo[odo_edges[i].first].om,
						poses_odo[odo_edges[i].first].fi,
						poses_odo[odo_edges[i].first].ka,
						poses_odo[odo_edges[i].second].px,
						poses_odo[odo_edges[i].second].py,
						poses_odo[odo_edges[i].second].pz,
						poses_odo[odo_edges[i].second].om,
						poses_odo[odo_edges[i].second].fi,
						poses_odo[odo_edges[i].second].ka);

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_tait_bryan_wc_case1(
						delta,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].om,
						poses[odo_edges[i].first].fi,
						poses[odo_edges[i].first].ka,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].om,
						poses[odo_edges[i].second].fi,
						poses[odo_edges[i].second].ka,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_tait_bryan_wc_case1_jacobian(jacobian,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].om,
						poses[odo_edges[i].first].fi,
						poses[odo_edges[i].first].ka,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].om,
						poses[odo_edges[i].second].fi,
						poses[odo_edges[i].second].ka);

				int ir = tripletListB.size();

				int ic_1 = odo_edges[i].first * 6;
				int ic_2 = odo_edges[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
			}

			for(size_t i = 0; i < vertices.size(); i++){
				TaitBryanPose pose_s = pose_tait_bryan_from_affine_matrix(vertices[i]);

				Eigen::Vector3d p_t;//(georeference_data[i].first(0,3), georeference_data[i].first(1,3), georeference_data[i].first(2,3));
				if(!find_nearest_neighbour(p_t, vertices[i], reference_points)){
					continue;
				}

				int ir = tripletListB.size();

				Eigen::Vector3d p_s(0,0,0);
				double delta_x;
				double delta_y;
				double delta_z;
				point_to_point_source_to_target_tait_bryan_wc(delta_x, delta_y, delta_z, pose_s.px, pose_s.py, pose_s.pz, pose_s.om, pose_s.fi, pose_s.ka, p_s.x(), p_s.y(), p_s.z(), p_t.x(), p_t.y(), p_t.z());

				Eigen::Matrix<double, 3, 6, Eigen::RowMajor> jacobian;
				point_to_point_source_to_target_tait_bryan_wc_jacobian(jacobian, pose_s.px, pose_s.py, pose_s.pz, pose_s.om, pose_s.fi, pose_s.ka, p_s.x(), p_s.y(), p_s.z());

				int ic = i * 6;
				tripletListA.emplace_back(ir     , ic + 0, -jacobian(0,0));
				tripletListA.emplace_back(ir + 1 , ic + 1, -jacobian(1,1));
				tripletListA.emplace_back(ir + 2 , ic + 2, -jacobian(2,2));

				tripletListP.emplace_back(ir    , ir    ,  1);
				tripletListP.emplace_back(ir + 1, ir + 1,  1);
				tripletListP.emplace_back(ir + 2, ir + 2,  1);

				tripletListB.emplace_back(ir    , 0,  delta_x);
				tripletListB.emplace_back(ir + 1, 0,  delta_y);
				tripletListB.emplace_back(ir + 2, 0,  delta_z);
			}


			Eigen::SparseMatrix<double> matA(tripletListB.size(), vertices.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(vertices.size() * 6 , vertices.size() * 6);
			Eigen::SparseMatrix<double> AtPB(vertices.size() * 6 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 6 * vertices.size()){
				int counter = 0;

				for(size_t i = 0; i < vertices.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(vertices[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];
					vertices[i] = affine_matrix_from_pose_tait_bryan(pose);
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

			std::vector<TaitBryanPose> poses;
			std::vector<TaitBryanPose> poses_odo;

			for(size_t i = 0 ; i < vertices.size(); i++){
				poses.push_back(pose_tait_bryan_from_affine_matrix(vertices[i]));
			}
			for(size_t i = 0 ; i < vertices_odo.size(); i++){
				poses_odo.push_back(pose_tait_bryan_from_affine_matrix(vertices_odo[i]));
			}

			for(size_t i = 0 ; i < odo_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_odo;
				relative_pose_tait_bryan_wc_case1(relative_pose_measurement_odo,
						poses_odo[odo_edges[i].first].px,
						poses_odo[odo_edges[i].first].py,
						poses_odo[odo_edges[i].first].pz,
						poses_odo[odo_edges[i].first].om,
						poses_odo[odo_edges[i].first].fi,
						poses_odo[odo_edges[i].first].ka,
						poses_odo[odo_edges[i].second].px,
						poses_odo[odo_edges[i].second].py,
						poses_odo[odo_edges[i].second].pz,
						poses_odo[odo_edges[i].second].om,
						poses_odo[odo_edges[i].second].fi,
						poses_odo[odo_edges[i].second].ka);

				Eigen::Matrix<double, 6, 1> delta;
				relative_pose_obs_eq_tait_bryan_wc_case1(
						delta,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].om,
						poses[odo_edges[i].first].fi,
						poses[odo_edges[i].first].ka,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].om,
						poses[odo_edges[i].second].fi,
						poses[odo_edges[i].second].ka,
						relative_pose_measurement_odo(0,0),
						relative_pose_measurement_odo(1,0),
						relative_pose_measurement_odo(2,0),
						relative_pose_measurement_odo(3,0),
						relative_pose_measurement_odo(4,0),
						relative_pose_measurement_odo(5,0));

				Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_tait_bryan_wc_case1_jacobian(jacobian,
						poses[odo_edges[i].first].px,
						poses[odo_edges[i].first].py,
						poses[odo_edges[i].first].pz,
						poses[odo_edges[i].first].om,
						poses[odo_edges[i].first].fi,
						poses[odo_edges[i].first].ka,
						poses[odo_edges[i].second].px,
						poses[odo_edges[i].second].py,
						poses[odo_edges[i].second].pz,
						poses[odo_edges[i].second].om,
						poses[odo_edges[i].second].fi,
						poses[odo_edges[i].second].ka);

				int ir = tripletListB.size();

				int ic_1 = odo_edges[i].first * 6;
				int ic_2 = odo_edges[i].second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
			}

			for(size_t i = 0; i < vertices.size(); i++){
				TaitBryanPose pose_s = pose_tait_bryan_from_affine_matrix(vertices[i]);

				Eigen::Vector3d p_t;
				if(!find_nearest_neighbour(p_t, vertices[i], reference_points)){
					continue;
				}

				int ir = tripletListB.size();
				int ic = i * 6;
				double delta;
				double a = 1;
				Eigen::Matrix<double, 1, 1> jacobian;
				observation_equation_constraint(delta, a, vertices[i](0,3), p_t.x());
				observation_equation_constraint_jacobian(jacobian, a, vertices[i](0,3), p_t.x());

				tripletListA.emplace_back(ir, ic,  -jacobian(0,0));
				tripletListP.emplace_back(ir, ir,  1);
				tripletListB.emplace_back(ir, 0,  delta);

				observation_equation_constraint(delta, a, vertices[i](1,3), p_t.y());
				observation_equation_constraint_jacobian(jacobian, a, vertices[i](1,3), p_t.y());

				tripletListA.emplace_back(ir + 1, ic + 1,  -jacobian(0,0));
				tripletListP.emplace_back(ir + 1, ir + 1,  1);
				tripletListB.emplace_back(ir + 1, 0,  delta);

				observation_equation_constraint(delta, a, vertices[i](2,3), p_t.z());
				observation_equation_constraint_jacobian(jacobian, a, vertices[i](2,3), p_t.z());

				tripletListA.emplace_back(ir + 2, ic + 2,  -jacobian(0,0));
				tripletListP.emplace_back(ir + 2, ir + 2,  1);
				tripletListB.emplace_back(ir + 2, 0,  delta);
			}


			Eigen::SparseMatrix<double> matA(tripletListB.size(), vertices.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(vertices.size() * 6 , vertices.size() * 6);
			Eigen::SparseMatrix<double> AtPB(vertices.size() * 6 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 6 * vertices.size()){
				int counter = 0;

				for(size_t i = 0; i < vertices.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(vertices[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];
					vertices[i] = affine_matrix_from_pose_tait_bryan(pose);
				}
				std::cout << "optimizing with tait bryan finished" << std::endl;
			}else{
				std::cout << "optimizing with tait bryan FAILED" << std::endl;
			}

			break;
		}
		case 'x':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<TaitBryanPose> poses;

			for(size_t i = 0 ; i < vertices.size(); i++){
				poses.push_back(pose_tait_bryan_from_affine_matrix(vertices[i]));
			}

			for(size_t i = 0; i < smoothness_indexes.size(); i++){
				Eigen::Matrix<double, 6, 1> delta;
				smoothness_obs_eq_tait_bryan_wc(delta,
						poses[smoothness_indexes[i].index_before].px,
						poses[smoothness_indexes[i].index_before].py,
						poses[smoothness_indexes[i].index_before].pz,
						poses[smoothness_indexes[i].index_before].om,
						poses[smoothness_indexes[i].index_before].fi,
						poses[smoothness_indexes[i].index_before].ka,
						poses[smoothness_indexes[i].index_curr].px,
						poses[smoothness_indexes[i].index_curr].py,
						poses[smoothness_indexes[i].index_curr].pz,
						poses[smoothness_indexes[i].index_curr].om,
						poses[smoothness_indexes[i].index_curr].fi,
						poses[smoothness_indexes[i].index_curr].ka,
						poses[smoothness_indexes[i].index_after].px,
						poses[smoothness_indexes[i].index_after].py,
						poses[smoothness_indexes[i].index_after].pz,
						poses[smoothness_indexes[i].index_after].om,
						poses[smoothness_indexes[i].index_after].fi,
						poses[smoothness_indexes[i].index_after].ka);

				Eigen::Matrix<double, 6, 18, Eigen::RowMajor> jacobian;
				smoothness_obs_eq_tait_bryan_wc_jacobian(jacobian,
						poses[smoothness_indexes[i].index_before].px,
						poses[smoothness_indexes[i].index_before].py,
						poses[smoothness_indexes[i].index_before].pz,
						poses[smoothness_indexes[i].index_before].om,
						poses[smoothness_indexes[i].index_before].fi,
						poses[smoothness_indexes[i].index_before].ka,
						poses[smoothness_indexes[i].index_curr].px,
						poses[smoothness_indexes[i].index_curr].py,
						poses[smoothness_indexes[i].index_curr].pz,
						poses[smoothness_indexes[i].index_curr].om,
						poses[smoothness_indexes[i].index_curr].fi,
						poses[smoothness_indexes[i].index_curr].ka,
						poses[smoothness_indexes[i].index_after].px,
						poses[smoothness_indexes[i].index_after].py,
						poses[smoothness_indexes[i].index_after].pz,
						poses[smoothness_indexes[i].index_after].om,
						poses[smoothness_indexes[i].index_after].fi,
						poses[smoothness_indexes[i].index_after].ka);

				int ir = tripletListB.size();

				int ic_1 = smoothness_indexes[i].index_before * 6;
				int ic_2 = smoothness_indexes[i].index_curr * 6;
				int ic_3 = smoothness_indexes[i].index_after * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));

					tripletListA.emplace_back(ir + row, ic_3    , -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_3 + 1, -jacobian(row,13));
					tripletListA.emplace_back(ir + row, ic_3 + 2, -jacobian(row,14));
					tripletListA.emplace_back(ir + row, ic_3 + 3, -jacobian(row,15));
					tripletListA.emplace_back(ir + row, ic_3 + 4, -jacobian(row,16));
					tripletListA.emplace_back(ir + row, ic_3 + 5, -jacobian(row,17));
				}
				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
			}

			for(size_t i = 0; i < vertices.size(); i++){
				TaitBryanPose pose_s = pose_tait_bryan_from_affine_matrix(vertices[i]);

				Eigen::Vector3d p_t;//(georeference_data[i].first(0,3), georeference_data[i].first(1,3), georeference_data[i].first(2,3));
				if(!find_nearest_neighbour(p_t, vertices[i], reference_points)){
					continue;
				}

				int ir = tripletListB.size();

				Eigen::Vector3d p_s(0,0,0);
				double delta_x;
				double delta_y;
				double delta_z;
				point_to_point_source_to_target_tait_bryan_wc(delta_x, delta_y, delta_z, pose_s.px, pose_s.py, pose_s.pz, pose_s.om, pose_s.fi, pose_s.ka, p_s.x(), p_s.y(), p_s.z(), p_t.x(), p_t.y(), p_t.z());

				Eigen::Matrix<double, 3, 6, Eigen::RowMajor> jacobian;
				point_to_point_source_to_target_tait_bryan_wc_jacobian(jacobian, pose_s.px, pose_s.py, pose_s.pz, pose_s.om, pose_s.fi, pose_s.ka, p_s.x(), p_s.y(), p_s.z());

				int ic = i * 6;
				tripletListA.emplace_back(ir     , ic + 0, -jacobian(0,0));
				tripletListA.emplace_back(ir + 1 , ic + 1, -jacobian(1,1));
				tripletListA.emplace_back(ir + 2 , ic + 2, -jacobian(2,2));

				tripletListP.emplace_back(ir    , ir    ,  1);
				tripletListP.emplace_back(ir + 1, ir + 1,  1);
				tripletListP.emplace_back(ir + 2, ir + 2,  1);

				tripletListB.emplace_back(ir    , 0,  delta_x);
				tripletListB.emplace_back(ir + 1, 0,  delta_y);
				tripletListB.emplace_back(ir + 2, 0,  delta_z);
			}


			Eigen::SparseMatrix<double> matA(tripletListB.size(), vertices.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(vertices.size() * 6 , vertices.size() * 6);
			Eigen::SparseMatrix<double> AtPB(vertices.size() * 6 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 6 * vertices.size()){
				int counter = 0;

				for(size_t i = 0; i < vertices.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(vertices[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];
					vertices[i] = affine_matrix_from_pose_tait_bryan(pose);
				}
				std::cout << "optimizing with tait bryan finished" << std::endl;
			}else{
				std::cout << "optimizing with tait bryan FAILED" << std::endl;
			}

			break;
		}
		case 'z':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<TaitBryanPose> poses;

			for(size_t i = 0 ; i < vertices.size(); i++){
				poses.push_back(pose_tait_bryan_from_affine_matrix(vertices[i]));
			}

			for(size_t i = 0; i < smoothness_indexes.size(); i++){
				Eigen::Matrix<double, 6, 1> delta;
				smoothness_obs_eq_tait_bryan_wc(delta,
						poses[smoothness_indexes[i].index_before].px,
						poses[smoothness_indexes[i].index_before].py,
						poses[smoothness_indexes[i].index_before].pz,
						poses[smoothness_indexes[i].index_before].om,
						poses[smoothness_indexes[i].index_before].fi,
						poses[smoothness_indexes[i].index_before].ka,
						poses[smoothness_indexes[i].index_curr].px,
						poses[smoothness_indexes[i].index_curr].py,
						poses[smoothness_indexes[i].index_curr].pz,
						poses[smoothness_indexes[i].index_curr].om,
						poses[smoothness_indexes[i].index_curr].fi,
						poses[smoothness_indexes[i].index_curr].ka,
						poses[smoothness_indexes[i].index_after].px,
						poses[smoothness_indexes[i].index_after].py,
						poses[smoothness_indexes[i].index_after].pz,
						poses[smoothness_indexes[i].index_after].om,
						poses[smoothness_indexes[i].index_after].fi,
						poses[smoothness_indexes[i].index_after].ka);

				Eigen::Matrix<double, 6, 18, Eigen::RowMajor> jacobian;
				smoothness_obs_eq_tait_bryan_wc_jacobian(jacobian,
						poses[smoothness_indexes[i].index_before].px,
						poses[smoothness_indexes[i].index_before].py,
						poses[smoothness_indexes[i].index_before].pz,
						poses[smoothness_indexes[i].index_before].om,
						poses[smoothness_indexes[i].index_before].fi,
						poses[smoothness_indexes[i].index_before].ka,
						poses[smoothness_indexes[i].index_curr].px,
						poses[smoothness_indexes[i].index_curr].py,
						poses[smoothness_indexes[i].index_curr].pz,
						poses[smoothness_indexes[i].index_curr].om,
						poses[smoothness_indexes[i].index_curr].fi,
						poses[smoothness_indexes[i].index_curr].ka,
						poses[smoothness_indexes[i].index_after].px,
						poses[smoothness_indexes[i].index_after].py,
						poses[smoothness_indexes[i].index_after].pz,
						poses[smoothness_indexes[i].index_after].om,
						poses[smoothness_indexes[i].index_after].fi,
						poses[smoothness_indexes[i].index_after].ka);

				int ir = tripletListB.size();

				int ic_1 = smoothness_indexes[i].index_before * 6;
				int ic_2 = smoothness_indexes[i].index_curr * 6;
				int ic_3 = smoothness_indexes[i].index_after * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,11));

					tripletListA.emplace_back(ir + row, ic_3    , -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_3 + 1, -jacobian(row,13));
					tripletListA.emplace_back(ir + row, ic_3 + 2, -jacobian(row,14));
					tripletListA.emplace_back(ir + row, ic_3 + 3, -jacobian(row,15));
					tripletListA.emplace_back(ir + row, ic_3 + 4, -jacobian(row,16));
					tripletListA.emplace_back(ir + row, ic_3 + 5, -jacobian(row,17));
				}
				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
			}

			for(size_t i = 0; i < vertices.size(); i++){
				TaitBryanPose pose_s = pose_tait_bryan_from_affine_matrix(vertices[i]);

				Eigen::Vector3d p_t;
				if(!find_nearest_neighbour(p_t, vertices[i], reference_points)){
					continue;
				}

				int ir = tripletListB.size();
				int ic = i * 6;
				double delta;
				double a = 1;
				Eigen::Matrix<double, 1, 1> jacobian;
				observation_equation_constraint(delta, a, vertices[i](0,3), p_t.x());
				observation_equation_constraint_jacobian(jacobian, a, vertices[i](0,3), p_t.x());

				tripletListA.emplace_back(ir, ic,  -jacobian(0,0));
				tripletListP.emplace_back(ir, ir,  1);
				tripletListB.emplace_back(ir, 0,  delta);

				observation_equation_constraint(delta, a, vertices[i](1,3), p_t.y());
				observation_equation_constraint_jacobian(jacobian, a, vertices[i](1,3), p_t.y());

				tripletListA.emplace_back(ir + 1, ic + 1,  -jacobian(0,0));
				tripletListP.emplace_back(ir + 1, ir + 1,  1);
				tripletListB.emplace_back(ir + 1, 0,  delta);

				observation_equation_constraint(delta, a, vertices[i](2,3), p_t.z());
				observation_equation_constraint_jacobian(jacobian, a, vertices[i](2,3), p_t.z());

				tripletListA.emplace_back(ir + 2, ic + 2,  -jacobian(0,0));
				tripletListP.emplace_back(ir + 2, ir + 2,  1);
				tripletListB.emplace_back(ir + 2, 0,  delta);
			}


			Eigen::SparseMatrix<double> matA(tripletListB.size(), vertices.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(vertices.size() * 6 , vertices.size() * 6);
			Eigen::SparseMatrix<double> AtPB(vertices.size() * 6 , 1);

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

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 6 * vertices.size()){
				int counter = 0;

				for(size_t i = 0; i < vertices.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(vertices[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];
					vertices[i] = affine_matrix_from_pose_tait_bryan(pose);
				}
				std::cout << "optimizing with tait bryan finished" << std::endl;
			}else{
				std::cout << "optimizing with tait bryan FAILED" << std::endl;
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
	std::cout << "t: optimize (relative pose constraint and source to target)" << std::endl;
	std::cout << "y: optimize (relative pose constraint and linear function constraint)" << std::endl;
	std::cout << "x: optimize (smoothness constraint and source to target)" << std::endl;
	std::cout << "z: optimize (smoothness constraint and linear function constraint)" << std::endl;
	std::cout << "n: add noise to mesh" << std::endl;
}







