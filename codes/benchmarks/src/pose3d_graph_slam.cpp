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

#include "manif/SE3.h"

using manif::SE3d;
using manif::SE3Tangentd;

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

int main(int argc, char *argv[]){
	if(argc != 2){
		std::cout << "USAGE: " << argv[0] << " file.g2o" << std::endl;
		return 1;
	}

	std::ifstream g2o_file(argv[1]);
    //std::ifstream g2o_file("../data/sphere_bignoise_vertex3.g2o");
	//std::ifstream g2o_file("../data/sphere.g2o");
	//std::ifstream g2o_file("../data/cubicle.g2o");
	//std::ifstream g2o_file("../data/sphere-optimized-g2o.g2o");
	//std::ifstream g2o_file("../data/rim.g2o");

	std::string line;
    while (std::getline(g2o_file, line)) {
        std::stringstream line_stream(line);
        std::string class_element;
        line_stream >> class_element;
        if (class_element == "VERTEX_SE3:QUAT"){
            int id=0;
            QuaternionPose p;
            line_stream>> id;
            line_stream >> p.px;
            line_stream >> p.py;
            line_stream >> p.pz;
            line_stream >> p.q1;
            line_stream >> p.q2;
            line_stream >> p.q3;
            line_stream >> p.q0;
            normalize_quaternion(p);

            Eigen::Affine3d m = affine_matrix_from_pose_quaternion(p);
            m_poses.push_back(m);
        }
        else if(class_element == "EDGE_SE3:QUAT"){
            int pose_id1,pose_id2;
            line_stream>>pose_id1;
            line_stream>>pose_id2;
            QuaternionPose p;
            line_stream >> p.px;
            line_stream >> p.py;
            line_stream >> p.pz;
            line_stream >> p.q1;
            line_stream >> p.q2;
            line_stream >> p.q3;
            line_stream >> p.q0;
            Eigen::Affine3d m = affine_matrix_from_pose_quaternion(p);

            edges_g2o.emplace_back(std::make_tuple(pose_id1,pose_id2, m));
        }
    }
    g2o_file.close();

    //exit(1);

//	for(size_t i = 0 ; i < 100; i++){
//		TaitBryanPose p;
//		p.px = i;
//		p.py = -1;
//		p.pz = 0.0;
//		p.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
//		p.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
//		p.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
//
//		Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(p);
//		m_poses.push_back(m);
//	}
//	for(size_t i = 0 ; i < 100; i++){
//		TaitBryanPose p;
//		p.px = i;
//		p.py = 1;
//		p.pz = 0.0;
//		p.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
//		p.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
//		p.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.01;
//
//		Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(p);
//		m_poses.push_back(m);
//	}
	//m_poses_desired = m_poses;

	//for(size_t i = 1; i < 100; i++){
	//	odo_edges.emplace_back(i-1,i);
	//}
//
//	for(size_t i = 101; i < 200; i++){
//		odo_edges.emplace_back(i-1,i);
//	}
//
//	for(size_t i = 0; i < 100; i+=10){
//		loop_edges.emplace_back(i,i+100);
//	}

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
	glutCreateWindow("relative_pose");
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

				//std::cout << from.px << " " << from.py << " " << from.pz << " " <<
				//		to.px << " " << to.py << " " << to.pz << " " <<
				//		pose_rel.px << " " << pose_rel.py << " " << pose_rel.pz << std::endl;


				//Eigen::Matrix<double, 6, 1> relative_pose_measurement_odo;
				//relative_pose_measurement_odo << pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka;

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

				int ic_1 = first * 6;
				int ic_2 = second * 6;

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

				std::cout << "delta " << delta(0,0) << " " << delta(1,0) << " " << delta(2,0) << " " <<
						delta(3,0) << " " << delta(4,0) << " " << delta(5,0) << std::endl;

				if(abs(first - second) == 1){
					tripletListP.emplace_back(ir ,    ir,     1000000000);
					tripletListP.emplace_back(ir + 1, ir + 1, 1000000000);
					tripletListP.emplace_back(ir + 2, ir + 2, 1000000000);
					tripletListP.emplace_back(ir + 3, ir + 3, 1000000000);
					tripletListP.emplace_back(ir + 4, ir + 4, 1000000000);
					tripletListP.emplace_back(ir + 5, ir + 5, 1000000000);
				}else{
					tripletListP.emplace_back(ir ,    ir,     cauchy(delta(0,0),1)*0.001);
					tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta(1,0),1)*0.001);
					tripletListP.emplace_back(ir + 2, ir + 2, cauchy(delta(2,0),1)*0.001);
					tripletListP.emplace_back(ir + 3, ir + 3, cauchy(delta(3,0),1)*0.001);
					tripletListP.emplace_back(ir + 4, ir + 4, cauchy(delta(4,0),1)*0.001);
					tripletListP.emplace_back(ir + 5, ir + 5, cauchy(delta(5,0),1)*0.001);
				}
			}

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir     , 0, 1);
			tripletListA.emplace_back(ir + 1 , 1, 1);
			tripletListA.emplace_back(ir + 2 , 2, 1);
			tripletListA.emplace_back(ir + 3 , 3, 1);
			tripletListA.emplace_back(ir + 4 , 4, 1);
			tripletListA.emplace_back(ir + 5 , 5, 1);

			tripletListP.emplace_back(ir     , ir,     1000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 1000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 1000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 1000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 1000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 1000000);

			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 6 , m_poses.size() * 6);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 6 , 1);

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
				for(size_t i = 0 ; i < m_poses.size() * 6; i++){
					tripletListA.emplace_back(i, i, 10);
				}
				Eigen::SparseMatrix<double> matLM(m_poses.size() * 6, m_poses.size() * 6);
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

			//for(size_t i = 0 ; i < h_x.size(); i++){
			//	std::cout << h_x[i] << std::endl;
			//}

			if(h_x.size() == 6 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++]*0.1;
					pose.py += h_x[counter++]*0.1;
					pose.pz += h_x[counter++]*0.1;
					pose.om += h_x[counter++]*0.1;
					pose.fi += h_x[counter++]*0.1;
					pose.ka += h_x[counter++]*0.1;
					m_poses[i] = affine_matrix_from_pose_tait_bryan(pose);
				}
				std::cout << "optimizing with tait bryan finished" << std::endl;
			}else{
				std::cout << "optimizing with tait bryan FAILED" << std::endl;
			}

			break;
		}

		case 'r':{

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<RodriguesPose> poses;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_rodrigues_from_affine_matrix(m_poses[i]));
			}

            for(size_t i = 0 ; i < edges_g2o.size(); i++){
                const int first = std::get<0>(edges_g2o[i]);
                const int second = std::get<1>(edges_g2o[i]);
                const Eigen::Affine3d& rel = std::get<2>(edges_g2o[i]);

                //Eigen::Matrix<double, 6, 1> relative_pose_measurement_loop;
                RodriguesPose ps = pose_rodrigues_from_affine_matrix(rel);
                //relative_pose_measurement_loop << ps.px, ps.py, ps.pz, ps.sx, ps.sy,ps.sz;
//                relative_pose_rodrigues_wc(relative_pose_measurement_loop,
//                                           poses_desired[first].px,
//                                           poses_desired[first].py,
//                                           poses_desired[first].pz,
//                                           poses_desired[first].sx,
//                                           poses_desired[first].sy,
//                                           poses_desired[first].sz,
//                                           poses_desired[second].px,
//                                           poses_desired[second].py,
//                                           poses_desired[second].pz,
//                                           poses_desired[second].sx,
//                                           poses_desired[second].sy,
//                                           poses_desired[second].sz);

                Eigen::Matrix<double, 6, 1> delta;
                relative_pose_obs_eq_rodrigues_wc(
                        delta,
                        poses[first].px,
                        poses[first].py,
                        poses[first].pz,
                        poses[first].sx,
                        poses[first].sy,
                        poses[first].sz,
                        poses[second].px,
                        poses[second].py,
                        poses[second].pz,
                        poses[second].sx,
                        poses[second].sy,
                        poses[second].sz,
                        ps.px,
                        ps.py,
                        ps.pz,
                        ps.sx,
                        ps.sy,
                        ps.sz);

                Eigen::Matrix<double, 6, 12, Eigen::RowMajor> jacobian;
                relative_pose_obs_eq_rodrigues_wc_jacobian(jacobian,
					   poses[first].px,
					   poses[first].py,
					   poses[first].pz,
					   poses[first].sx,
					   poses[first].sy,
					   poses[first].sz,
					   poses[second].px,
					   poses[second].py,
					   poses[second].pz,
					   poses[second].sx,
					   poses[second].sy,
					   poses[second].sz);

                int ir = tripletListB.size();

                int ic_1 = first * 6;
                int ic_2 = second * 6;

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


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 6 , m_poses.size() * 6);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 6 , 1);

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

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}
			std::cout << "h_x.size(): " << h_x.size() << std::endl;
			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			if(h_x.size() == 6 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					RodriguesPose pose = pose_rodrigues_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.sx += h_x[counter++];
					pose.sy += h_x[counter++];
					pose.sz += h_x[counter++];

					m_poses[i] = affine_matrix_from_pose_rodrigues(pose);

					Eigen::Vector3d vx(m_poses[i](0,0), m_poses[i](1,0), m_poses[i](2,0));
					Eigen::Vector3d vy(m_poses[i](0,1), m_poses[i](1,1), m_poses[i](2,1));
					Eigen::Vector3d vz(m_poses[i](0,2), m_poses[i](1,2), m_poses[i](2,2));

					//std::cout << std::setprecision(15);
					//std::cout << "norm: "<< vx.norm() << " " << vy.norm() << " " << vz.norm() << " " <<
					//		vx.dot(vy) << " " << vy.dot(vz) << " " << vx.dot(vz) << std::endl;
				}
				std::cout << "optimizing with rodrigues finished" << std::endl;
			}else{
				std::cout << "optimizing with rodrigues FAILED" << std::endl;
			}
			break;
		}
		case 'q':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<QuaternionPose> poses;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_quaternion_from_affine_matrix(m_poses[i]));
			}

			for(size_t i = 0 ; i < edges_g2o.size(); i++){
				const int first = std::get<0>(edges_g2o[i]);
				const int second = std::get<1>(edges_g2o[i]);
				const Eigen::Affine3d& rel = std::get<2>(edges_g2o[i]);
				QuaternionPose pose_rel = pose_quaternion_from_affine_matrix(rel);

				//QuaternionPose from = poses[first];
				//QuaternionPose to = poses[second];

				Eigen::Matrix<double, 7, 1> delta;
				relative_pose_obs_eq_quaternion_wc(
					delta,
					poses[first].px,
					poses[first].py,
					poses[first].pz,
					poses[first].q0,
					poses[first].q1,
					poses[first].q2,
					poses[first].q3,
					poses[second].px,
					poses[second].py,
					poses[second].pz,
					poses[second].q0,
					poses[second].q1,
					poses[second].q2,
					poses[second].q3,
					pose_rel.px,
					pose_rel.py,
					pose_rel.pz,
					pose_rel.q0,
					pose_rel.q1,
					pose_rel.q2,
					pose_rel.q3);

				Eigen::Matrix<double, 7, 14, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_quaternion_wc_jacobian(jacobian,
					poses[first].px,
					poses[first].py,
					poses[first].pz,
					poses[first].q0,
					poses[first].q1,
					poses[first].q2,
					poses[first].q3,
					poses[second].px,
					poses[second].py,
					poses[second].pz,
					poses[second].q0,
					poses[second].q1,
					poses[second].q2,
					poses[second].q3);

				int ir = tripletListB.size();

				int ic_1 = first * 7;
				int ic_2 = second * 7;

				for(size_t row = 0 ; row < 7; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));
					tripletListA.emplace_back(ir + row, ic_1 + 6, -jacobian(row,6));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,11));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_2 + 6, -jacobian(row,13));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));
				tripletListB.emplace_back(ir + 6, 0, delta(6,0));

				tripletListP.emplace_back(ir ,    ir,     1);
				tripletListP.emplace_back(ir + 1, ir + 1, 1);
				tripletListP.emplace_back(ir + 2, ir + 2, 1);
				tripletListP.emplace_back(ir + 3, ir + 3, 1);
				tripletListP.emplace_back(ir + 4, ir + 4, 1);
				tripletListP.emplace_back(ir + 5, ir + 5, 1);
				tripletListP.emplace_back(ir + 6, ir + 6, 1);
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

			for(size_t i = 0 ; i < m_poses.size(); i++){
				int ic = i * 7;
				ir = tripletListB.size();
				QuaternionPose pose = pose_quaternion_from_affine_matrix(m_poses[i]);

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


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 7);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 7 , m_poses.size() * 7);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 7 , 1);

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

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			std::cout << "h_x.size(): " << h_x.size() << std::endl;

			std::cout << "AtPA=AtPB SOLVED" << std::endl;



			if(h_x.size() == 7 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					QuaternionPose pose = pose_quaternion_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.q0 += h_x[counter++];
					pose.q1 += h_x[counter++];
					pose.q2 += h_x[counter++];
					pose.q3 += h_x[counter++];
					m_poses[i] = affine_matrix_from_pose_quaternion(pose);
				}
				std::cout << "optimizing with quaternions finished" << std::endl;
			}else{
				std::cout << "optimizing with quaternions FAILED" << std::endl;
			}
			break;
		}
		case 'x':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < edges_g2o.size(); i++){
				const int first = std::get<0>(edges_g2o[i]);
				const int second = std::get<1>(edges_g2o[i]);
				const Eigen::Affine3d& rel = std::get<2>(edges_g2o[i]);

				Eigen::Matrix<double, 12, 1> delta;
				relative_pose_obs_eq_wc(
						delta,
						m_poses[first](0,3),
						m_poses[first](1,3),
						m_poses[first](2,3),
						m_poses[first](0,0),
						m_poses[first](0,1),
						m_poses[first](0,2),
						m_poses[first](1,0),
						m_poses[first](1,1),
						m_poses[first](1,2),
						m_poses[first](2,0),
						m_poses[first](2,1),
						m_poses[first](2,2),
						m_poses[second](0,3),
						m_poses[second](1,3),
						m_poses[second](2,3),
						m_poses[second](0,0),
						m_poses[second](0,1),
						m_poses[second](0,2),
						m_poses[second](1,0),
						m_poses[second](1,1),
						m_poses[second](1,2),
						m_poses[second](2,0),
						m_poses[second](2,1),
						m_poses[second](2,2),
						rel(0,3),
						rel(1,3),
						rel(2,3),
						rel(0,0),
						rel(0,1),
						rel(0,2),
						rel(1,0),
						rel(1,1),
						rel(1,2),
						rel(2,0),
						rel(2,1),
						rel(2,2));

				Eigen::Matrix<double, 12, 24, Eigen::RowMajor> jacobian;
				relative_pose_obs_eq_wc_jacobian(jacobian,
						m_poses[first](0,3),
						m_poses[first](1,3),
						m_poses[first](2,3),
						m_poses[first](0,0),
						m_poses[first](0,1),
						m_poses[first](0,2),
						m_poses[first](1,0),
						m_poses[first](1,1),
						m_poses[first](1,2),
						m_poses[first](2,0),
						m_poses[first](2,1),
						m_poses[first](2,2),
						m_poses[second](0,3),
						m_poses[second](1,3),
						m_poses[second](2,3),
						m_poses[second](0,0),
						m_poses[second](0,1),
						m_poses[second](0,2),
						m_poses[second](1,0),
						m_poses[second](1,1),
						m_poses[second](1,2),
						m_poses[second](2,0),
						m_poses[second](2,1),
						m_poses[second](2,2));

				int ir = tripletListB.size();
				int ic_1 = first * 12;
				int ic_2 = second * 12;

				for(size_t row = 0 ; row < 12; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -jacobian(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -jacobian(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -jacobian(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -jacobian(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -jacobian(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -jacobian(row,5));
					tripletListA.emplace_back(ir + row, ic_1 + 6, -jacobian(row,6));
					tripletListA.emplace_back(ir + row, ic_1 + 7, -jacobian(row,7));
					tripletListA.emplace_back(ir + row, ic_1 + 8, -jacobian(row,8));
					tripletListA.emplace_back(ir + row, ic_1 + 9, -jacobian(row,9));
					tripletListA.emplace_back(ir + row, ic_1 + 10, -jacobian(row,10));
					tripletListA.emplace_back(ir + row, ic_1 + 11, -jacobian(row,11));

					tripletListA.emplace_back(ir + row, ic_2    , -jacobian(row,12));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -jacobian(row,13));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -jacobian(row,14));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -jacobian(row,15));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -jacobian(row,16));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -jacobian(row,17));
					tripletListA.emplace_back(ir + row, ic_2 + 6, -jacobian(row,18));
					tripletListA.emplace_back(ir + row, ic_2 + 7, -jacobian(row,19));
					tripletListA.emplace_back(ir + row, ic_2 + 8, -jacobian(row,20));
					tripletListA.emplace_back(ir + row, ic_2 + 9, -jacobian(row,21));
					tripletListA.emplace_back(ir + row, ic_2 + 10, -jacobian(row,22));
					tripletListA.emplace_back(ir + row, ic_2 + 11, -jacobian(row,23));
				}

				tripletListB.emplace_back(ir,     0, delta(0,0));
				tripletListB.emplace_back(ir + 1, 0, delta(1,0));
				tripletListB.emplace_back(ir + 2, 0, delta(2,0));
				tripletListB.emplace_back(ir + 3, 0, delta(3,0));
				tripletListB.emplace_back(ir + 4, 0, delta(4,0));
				tripletListB.emplace_back(ir + 5, 0, delta(5,0));
				tripletListB.emplace_back(ir + 6, 0, delta(6,0));
				tripletListB.emplace_back(ir + 7, 0, delta(7,0));
				tripletListB.emplace_back(ir + 8, 0, delta(8,0));
				tripletListB.emplace_back(ir + 9, 0, delta(9,0));
				tripletListB.emplace_back(ir + 10, 0, delta(10,0));
				tripletListB.emplace_back(ir + 11, 0, delta(11,0));

				if(abs(first - second) == 1){
					tripletListP.emplace_back(ir ,    ir,     1000);
					tripletListP.emplace_back(ir + 1, ir + 1, 1000);
					tripletListP.emplace_back(ir + 2, ir + 2, 1000);
					tripletListP.emplace_back(ir + 3, ir + 3, 1000);
					tripletListP.emplace_back(ir + 4, ir + 4, 1000);
					tripletListP.emplace_back(ir + 5, ir + 5, 1000);
					tripletListP.emplace_back(ir + 6, ir + 6, 1000);
					tripletListP.emplace_back(ir + 7, ir + 7, 1000);
					tripletListP.emplace_back(ir + 8, ir + 8, 1000);
					tripletListP.emplace_back(ir + 9, ir + 9, 1000);
					tripletListP.emplace_back(ir + 10, ir + 10, 1000);
					tripletListP.emplace_back(ir + 11, ir + 11, 1000);
				}else{
					tripletListP.emplace_back(ir ,    ir,     1);
					tripletListP.emplace_back(ir + 1, ir + 1, 1);
					tripletListP.emplace_back(ir + 2, ir + 2, 1);
					tripletListP.emplace_back(ir + 3, ir + 3, 1);
					tripletListP.emplace_back(ir + 4, ir + 4, 1);
					tripletListP.emplace_back(ir + 5, ir + 5, 1);
					tripletListP.emplace_back(ir + 6, ir + 6, 1);
					tripletListP.emplace_back(ir + 7, ir + 7, 1);
					tripletListP.emplace_back(ir + 8, ir + 8, 1);
					tripletListP.emplace_back(ir + 9, ir + 9, 1);
					tripletListP.emplace_back(ir + 10, ir + 10, 1);
					tripletListP.emplace_back(ir + 11, ir + 11, 1);
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
			tripletListA.emplace_back(ir + 7 , 7, 1);
			tripletListA.emplace_back(ir + 8 , 8, 1);
			tripletListA.emplace_back(ir + 9 , 9, 1);
			tripletListA.emplace_back(ir + 10 , 10, 1);
			tripletListA.emplace_back(ir + 11 , 11, 1);

			tripletListP.emplace_back(ir     , ir,     10000000000000);
			tripletListP.emplace_back(ir + 1 , ir + 1, 10000000000000);
			tripletListP.emplace_back(ir + 2 , ir + 2, 10000000000000);
			tripletListP.emplace_back(ir + 3 , ir + 3, 10000000000000);
			tripletListP.emplace_back(ir + 4 , ir + 4, 10000000000000);
			tripletListP.emplace_back(ir + 5 , ir + 5, 10000000000000);
			tripletListP.emplace_back(ir + 6 , ir + 6, 10000000000000);
			tripletListP.emplace_back(ir + 7 , ir + 7, 10000000000000);
			tripletListP.emplace_back(ir + 8 , ir + 8, 10000000000000);
			tripletListP.emplace_back(ir + 9 , ir + 9, 10000000000000);
			tripletListP.emplace_back(ir + 10 , ir + 10, 10000000000000);
			tripletListP.emplace_back(ir + 11 , ir + 11, 10000000000000);


			tripletListB.emplace_back(ir     , 0, 0);
			tripletListB.emplace_back(ir + 1 , 0, 0);
			tripletListB.emplace_back(ir + 2 , 0, 0);
			tripletListB.emplace_back(ir + 3 , 0, 0);
			tripletListB.emplace_back(ir + 4 , 0, 0);
			tripletListB.emplace_back(ir + 5 , 0, 0);
			tripletListB.emplace_back(ir + 6 , 0, 0);
			tripletListB.emplace_back(ir + 7 , 0, 0);
			tripletListB.emplace_back(ir + 8 , 0, 0);
			tripletListB.emplace_back(ir + 9 , 0, 0);
			tripletListB.emplace_back(ir + 10 , 0, 0);
			tripletListB.emplace_back(ir + 11 , 0, 0);


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 12);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 12 , m_poses.size() * 12);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 12 , 1);

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
				for(size_t i = 0 ; i < m_poses.size() * 12; i++){
					tripletListA.emplace_back(i, i, 10);
				}
				Eigen::SparseMatrix<double> matLM(m_poses.size() * 12, m_poses.size() * 12);
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

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 12 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					m_poses[i](0,3) += h_x[counter++];
					m_poses[i](1,3) += h_x[counter++];
					m_poses[i](2,3) += h_x[counter++];
					m_poses[i](0,0) += h_x[counter++];
					m_poses[i](0,1) += h_x[counter++];
					m_poses[i](0,2) += h_x[counter++];
					m_poses[i](1,0) += h_x[counter++];
					m_poses[i](1,1) += h_x[counter++];
					m_poses[i](1,2) += h_x[counter++];
					m_poses[i](2,0) += h_x[counter++];
					m_poses[i](2,1) += h_x[counter++];
					m_poses[i](2,2) += h_x[counter++];

					orthogonalize_rotation(m_poses[i]);
				}
				std::cout << "optimizing without rotation matrix parametrization finished" << std::endl;
			}else{
				std::cout << "optimizing without rotation matrix parametrization FAILED" << std::endl;
			}
			break;
		}
		case 'y':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<SE3d> X;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				Eigen::Vector3d t(m_poses[i](0,3), m_poses[i](1,3), m_poses[i](2,3));
				Eigen::Quaterniond q(m_poses[i].rotation());
				X.push_back(SE3d(t,q));
			}

			for(size_t i = 0 ; i < edges_g2o.size(); i++){
				const int first = std::get<0>(edges_g2o[i]);
				const int second = std::get<1>(edges_g2o[i]);
				const Eigen::Affine3d& rel = std::get<2>(edges_g2o[i]);
				Eigen::Vector3d t(rel(0,3), rel(1,3), rel(2,3));
				Eigen::Quaterniond q(rel.rotation());
				SE3d U = SE3d(t,q);

				SE3Tangentd     d;
				SE3Tangentd     u;
				Eigen::Matrix<double, 6, 6>         J_d_xi, J_d_xj;

				SE3d         Xi,
							 Xj;

				Xi = X[first];
				Xj = X[second];

				d  = Xj.rminus(Xi, J_d_xj, J_d_xi);
				u = U.log();

				int ir = tripletListB.size();
				int ic_1 = first * 6;
				int ic_2 = second * 6;

				for(size_t row = 0 ; row < 6; row ++){
					tripletListA.emplace_back(ir + row, ic_1    , -J_d_xi(row,0));
					tripletListA.emplace_back(ir + row, ic_1 + 1, -J_d_xi(row,1));
					tripletListA.emplace_back(ir + row, ic_1 + 2, -J_d_xi(row,2));
					tripletListA.emplace_back(ir + row, ic_1 + 3, -J_d_xi(row,3));
					tripletListA.emplace_back(ir + row, ic_1 + 4, -J_d_xi(row,4));
					tripletListA.emplace_back(ir + row, ic_1 + 5, -J_d_xi(row,5));

					tripletListA.emplace_back(ir + row, ic_2    , -J_d_xj(row,0));
					tripletListA.emplace_back(ir + row, ic_2 + 1, -J_d_xj(row,1));
					tripletListA.emplace_back(ir + row, ic_2 + 2, -J_d_xj(row,2));
					tripletListA.emplace_back(ir + row, ic_2 + 3, -J_d_xj(row,3));
					tripletListA.emplace_back(ir + row, ic_2 + 4, -J_d_xj(row,4));
					tripletListA.emplace_back(ir + row, ic_2 + 5, -J_d_xj(row,5));
				}

				SE3Tangentd delta = d - u;

				tripletListB.emplace_back(ir,     0, delta.coeffs()(0));
				tripletListB.emplace_back(ir + 1, 0, delta.coeffs()(1));
				tripletListB.emplace_back(ir + 2, 0, delta.coeffs()(2));
				tripletListB.emplace_back(ir + 3, 0, delta.coeffs()(3));
				tripletListB.emplace_back(ir + 4, 0, delta.coeffs()(4));
				tripletListB.emplace_back(ir + 5, 0, delta.coeffs()(5));

				if(abs(first - second) == 1){
					tripletListP.emplace_back(ir ,    ir,     1);
					tripletListP.emplace_back(ir + 1, ir + 1, 1);
					tripletListP.emplace_back(ir + 2, ir + 2, 1);
					tripletListP.emplace_back(ir + 3, ir + 3, 1);
					tripletListP.emplace_back(ir + 4, ir + 4, 1);
					tripletListP.emplace_back(ir + 5, ir + 5, 1);
				}else{
					tripletListP.emplace_back(ir ,    ir,     cauchy(delta.coeffs()(0),1));
					tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta.coeffs()(1),1));
					tripletListP.emplace_back(ir + 2, ir + 2, cauchy(delta.coeffs()(2),1));
					tripletListP.emplace_back(ir + 3, ir + 3, cauchy(delta.coeffs()(3),1));
					tripletListP.emplace_back(ir + 4, ir + 4, cauchy(delta.coeffs()(4),1));
					tripletListP.emplace_back(ir + 5, ir + 5, cauchy(delta.coeffs()(5),1));
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


			Eigen::SparseMatrix<double> matA(tripletListB.size(), m_poses.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


			Eigen::SparseMatrix<double> AtPA(m_poses.size() * 6 , m_poses.size() * 6);
			Eigen::SparseMatrix<double> AtPB(m_poses.size() * 6 , 1);

			{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP) * matA;
			AtPB = (AtP) * matB;
			}

			tripletListA.clear();
			tripletListP.clear();
			tripletListB.clear();

			//std::cout << "AtPA.size: " << AtPA.size() << std::endl;
			//std::cout << "AtPB.size: " << AtPB.size() << std::endl;

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



			//for(size_t i = 0 ; i < h_x.size(); i++){
			//	std::cout << h_x[i] << std::endl;
			//}

			if(X.size() * 6 == h_x.size()){

				int counter = 0;
				for(size_t i = 0 ; i < X.size(); i++){
					//int dx_row          = i * 6;
					SE3Tangentd     dx;
					dx.coeffs()(0) = h_x[counter++];
					dx.coeffs()(1) = h_x[counter++];
					dx.coeffs()(2) = h_x[counter++];
					dx.coeffs()(3) = h_x[counter++];
					dx.coeffs()(4) = h_x[counter++];
					dx.coeffs()(5) = h_x[counter++];

					//dx = dX.segment<6>(dx_row);
					//std::cout << "before: " << X[i] << std::endl;
					X[i] = X[i] +  dx;
					//std::cout << "after: " << X[i] << std::endl;
				}

				for (int i = 0 ; i < m_poses.size(); i++){
					m_poses[i](0,0) = X[i].rotation()(0,0);
					m_poses[i](1,0) = X[i].rotation()(1,0);
					m_poses[i](2,0) = X[i].rotation()(2,0);

					m_poses[i](0,1) = X[i].rotation()(0,1);
					m_poses[i](1,1) = X[i].rotation()(1,1);
					m_poses[i](2,1) = X[i].rotation()(2,1);

					m_poses[i](0,2) = X[i].rotation()(0,2);
					m_poses[i](1,2) = X[i].rotation()(1,2);
					m_poses[i](2,2) = X[i].rotation()(2,2);

					m_poses[i](0,3) = X[i].translation()(0);
					m_poses[i](1,3) = X[i].translation()(1);
					m_poses[i](2,3) = X[i].translation()(2);
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
	std::cout << "r: optimize (Rodriguez)" << std::endl;
	std::cout << "q: optimize (Quaternion)" << std::endl;
	std::cout << "x: optimize (Without rotation matrix parametrization)" << std::endl;
	std::cout << "y: optimize (Lie Algebra: manif lib)" << std::endl;
}
