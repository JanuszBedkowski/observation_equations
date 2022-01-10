#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "relative_pose_tait_bryan_wc_jacobian.h"

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
std::vector<Eigen::Affine3d> m_poses_desired;
std::vector<std::pair<int, int>> odo_edges;


void calculate_cov(std::vector<Eigen::Affine3d> m_poses,
				   std::vector<Eigen::Affine3d> m_poses_desired,
				   std::vector<std::pair<int, int>> odo_edges,
				   Eigen::MatrixXd & cov_b);

void draw_ellipse(const Eigen::Matrix3d& covar, Eigen::Vector3d& mean, Eigen::Vector3f color, float nstd  = 3)
{
    Eigen::LLT<Eigen::Matrix<double,3,3> > cholSolver(covar);
    Eigen::Matrix3d transform = cholSolver.matrixL();

    const double pi = 3.141592;
    const double di =0.02;
    const double dj =0.04;
    const double du =di*2*pi;
    const double dv =dj*pi;
    glColor3f(color.x(), color.y(),color.z());

    for (double i = 0; i < 1.0; i+=di)  //horizonal
    {
        for (double j = 0; j < 1.0; j+=dj)  //vertical
        {
            double u = i*2*pi;      //0     to  2pi
            double v = (j-0.5)*pi;  //-pi/2 to pi/2

            const Eigen::Vector3d pp0( cos(v)* cos(u),cos(v) * sin(u),sin(v));
            const Eigen::Vector3d pp1(cos(v) * cos(u + du) ,cos(v) * sin(u + du) ,sin(v));
            const Eigen::Vector3d pp2(cos(v + dv)* cos(u + du) ,cos(v + dv)* sin(u + du) ,sin(v + dv));
            const Eigen::Vector3d pp3( cos(v + dv)* cos(u),cos(v + dv)* sin(u),sin(v + dv));
            Eigen::Vector3d tp0 = transform * (nstd*pp0) + mean;
            Eigen::Vector3d tp1 = transform * (nstd*pp1) + mean;
            Eigen::Vector3d tp2 = transform * (nstd*pp2) + mean;
            Eigen::Vector3d tp3 = transform * (nstd*pp3) + mean;

            glBegin(GL_LINE_LOOP);
            glVertex3dv(tp0.data());
            glVertex3dv(tp1.data());
            glVertex3dv(tp2.data());
            glVertex3dv(tp3.data());
            glEnd();
        }
    }
}

void draw_ellipse2D(const Eigen::Matrix3d& covar, Eigen::Vector3d& mean, Eigen::Vector3f color, float nstd  = 3)
{
    Eigen::LLT<Eigen::Matrix<double,3,3> > cholSolver(covar);
    Eigen::Matrix3d transform = cholSolver.matrixL();

    const double pi = 3.141592;
    const double di = 0.02;
    const double dj = 0.04;
    const double du = di*2*pi;
    const double dv = dj*pi;
    glColor3f(color.x(), color.y(),color.z());

    for (double i = 0; i < 1.0; i+=di) { //horizonal
		double u = i*2*pi;      //0     to  2pi
		const Eigen::Vector3d pp0( cos(u), sin (u),0);
		const Eigen::Vector3d pp1( cos(u+du), sin(u+du),0);
		Eigen::Vector3d tp0 = transform * (nstd*pp0) + mean;
		Eigen::Vector3d tp1 = transform * (nstd*pp1) + mean;
		glBegin(GL_LINE_LOOP);
		glVertex3dv(tp0.data());
		glVertex3dv(tp1.data());
		glEnd();
	}
}


int main(int argc, char *argv[]){
	TaitBryanPose p0;
	p0.px = 0;
	p0.py = 0;
	p0.pz = 0;
	p0.om = 0;
	p0.fi = 0;
	p0.ka = 0;
	Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(p0);
	m_poses.push_back(m);

	TaitBryanPose p_rel_forward_x;
	p_rel_forward_x.px = 1.0;
	p_rel_forward_x.py = 0.0;
	p_rel_forward_x.pz = 0.0;
	p_rel_forward_x.om = 0.0;
	p_rel_forward_x.fi = 0.0;
	p_rel_forward_x.ka = 0.0;

	TaitBryanPose p_rel_rotate;
	p_rel_rotate.px = 0.0;
	p_rel_rotate.py = 0.0;
	p_rel_rotate.pz = 0.0;
	p_rel_rotate.om = 0.0;
	p_rel_rotate.fi = 0.0;
	p_rel_rotate.ka = 90 * M_PI/180.0;

	Eigen::Affine3d m_rel;

	for(int c = 0; c < 4; c++){
		for(size_t i = 0 ; i < 3; i++){
			m_rel = affine_matrix_from_pose_tait_bryan(p_rel_forward_x);
			m = m * m_rel;
			m_poses.push_back(m);
		}
		m_rel = affine_matrix_from_pose_tait_bryan(p_rel_rotate);
		m = m * m_rel;
	}
	m_poses.pop_back();
	m_poses_desired = m_poses;
	m_poses.clear();
	p0.px = 0;
	p0.py = 0;
	p0.pz = 0;
	p0.om = 0;
	p0.fi = 0;
	p0.ka = 30 * M_PI/180;
	m = affine_matrix_from_pose_tait_bryan(p0);
	m_poses.push_back(m);

	p_rel_forward_x.px = 1.1;
	p_rel_forward_x.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
	p_rel_forward_x.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
	p_rel_forward_x.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
	p_rel_forward_x.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
	p_rel_forward_x.ka = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;

	p_rel_rotate.px = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
	p_rel_rotate.py = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
	p_rel_rotate.pz = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
	p_rel_rotate.om = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
	p_rel_rotate.fi = ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
	p_rel_rotate.ka = 80 * M_PI/180.0;


	for(int c = 0; c < 4; c++){
		for(size_t i = 0 ; i < 3; i++){
			m_rel = affine_matrix_from_pose_tait_bryan(p_rel_forward_x);
			m = m * m_rel;
			m_poses.push_back(m);
			p_rel_forward_x.px += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
			p_rel_forward_x.py += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
			p_rel_forward_x.pz += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			p_rel_forward_x.om += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			p_rel_forward_x.fi += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
			p_rel_forward_x.ka += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
		}
		m_rel = affine_matrix_from_pose_tait_bryan(p_rel_rotate);
		m = m * m_rel;
		p_rel_rotate.px += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
		p_rel_rotate.py += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
		p_rel_rotate.pz += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
		p_rel_rotate.om += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
		p_rel_rotate.fi += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
		p_rel_rotate.ka += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.000001;
	}
	m_poses.pop_back();
	for(size_t i = 1; i < m_poses.size(); i++){
		odo_edges.emplace_back(i-1,i);
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
	glutCreateWindow("relative_pose_covariances");
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

	/*glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(1.0f, 0.0f, 0.0f);

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 1.0f, 0.0f);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 1.0f);
	glEnd();*/

	glColor3f(1,0,0);
	glBegin(GL_LINES);
	for(size_t i = 0; i < odo_edges.size(); i++){
		glVertex3f(m_poses[odo_edges[i].first](0,3), m_poses[odo_edges[i].first](1,3), m_poses[odo_edges[i].first](2,3) );
		glVertex3f(m_poses[odo_edges[i].second](0,3), m_poses[odo_edges[i].second](1,3), m_poses[odo_edges[i].second](2,3) );
	}
	glEnd();

	Eigen::MatrixXd cov_b(m_poses.size() * 2, m_poses.size() * 2);
	cov_b = Eigen::MatrixXd::Zero(m_poses.size() * 2, m_poses.size() * 2);
	calculate_cov( m_poses,	m_poses_desired, odo_edges, cov_b);

	for(size_t i = 0; i < m_poses.size(); i++){
		Eigen::Vector3d mean(m_poses[i](0,3), m_poses[i](1,3), m_poses[i](2,3));
		Eigen::Matrix3d cov;
		int r = i * 6;
		int c = i * 6;
		cov(0,0) = cov_b(r+0,c+0);
		cov(0,1) = cov_b(r+0,c+1);
		cov(0,2) = cov_b(r+0,c+2);
		cov(1,0) = cov_b(r+1,c+0);
		cov(1,1) = cov_b(r+1,c+1);
		cov(1,2) = cov_b(r+1,c+2);
		cov(2,0) = cov_b(r+2,c+0);
		cov(2,1) = cov_b(r+2,c+1);
		cov(2,2) = cov_b(r+2,c+2);
		draw_ellipse2D(cov, mean, Eigen::Vector3f(1,0,0),1);
		draw_ellipse2D(cov, mean, Eigen::Vector3f(0,1,0),2);
		draw_ellipse2D(cov, mean, Eigen::Vector3f(0,0,1),3);
	}
	glutSwapBuffers();
}

void keyboard(unsigned char key, int /*x*/, int /*y*/) {
	switch (key) {
		case (27): {
			glutDestroyWindow(glutGetWindow());
			return;
		}
		case 'e':{
			odo_edges.emplace_back(1, m_poses.size()-2);
			break;
		}
		case 'n':{
			for(size_t i = 0 ; i < m_poses.size(); i++){
				TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m_poses[i]);
				pose.px += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.py += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.1;
				pose.pz += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.001;
				pose.om += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.0001;
				pose.fi += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.0001;
				pose.ka += ((float(rand()%1000000))/1000000.0f - 0.5) * 0.0001;
				m_poses[i] = affine_matrix_from_pose_tait_bryan(pose);
			}
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			std::vector<TaitBryanPose> poses;
			std::vector<TaitBryanPose> poses_desired;

			for(size_t i = 0 ; i < m_poses.size(); i++){
				poses.push_back(pose_tait_bryan_from_affine_matrix(m_poses[i]));
			}
			for(size_t i = 0 ; i < m_poses_desired.size(); i++){
				poses_desired.push_back(pose_tait_bryan_from_affine_matrix(m_poses_desired[i]));
			}

			for(size_t i = 0 ; i < odo_edges.size(); i++){
				Eigen::Matrix<double, 6, 1> relative_pose_measurement_odo;
				relative_pose_tait_bryan_wc_case1(relative_pose_measurement_odo,
						poses_desired[odo_edges[i].first].px,
						poses_desired[odo_edges[i].first].py,
						poses_desired[odo_edges[i].first].pz,
						poses_desired[odo_edges[i].first].om,
						poses_desired[odo_edges[i].first].fi,
						poses_desired[odo_edges[i].first].ka,
						poses_desired[odo_edges[i].second].px,
						poses_desired[odo_edges[i].second].py,
						poses_desired[odo_edges[i].second].pz,
						poses_desired[odo_edges[i].second].om,
						poses_desired[odo_edges[i].second].fi,
						poses_desired[odo_edges[i].second].ka);

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

				float angle_diff = delta(3,0);
				if(fabs(angle_diff) > M_PI){
					angle_diff -= 2.0*M_PI;
				}
				if(fabs(angle_diff)< -M_PI){
					angle_diff += 2.0*M_PI;
				}
				tripletListB.emplace_back(ir + 3, 0, angle_diff);

				angle_diff = delta(4,0);
				if(fabs(angle_diff) > M_PI){
					angle_diff -= 2.0*M_PI;
				}
				if(fabs(angle_diff)< -M_PI){
					angle_diff += 2.0*M_PI;
				}
				tripletListB.emplace_back(ir + 4, 0, angle_diff);

				angle_diff = delta(5,0);
				if(fabs(angle_diff) > M_PI){
					angle_diff -= 2.0*M_PI;
				}
				if(fabs(angle_diff)< -M_PI){
					angle_diff += 2.0*M_PI;
				}
				tripletListB.emplace_back(ir + 5, 0, angle_diff);

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

			tripletListP.emplace_back(ir     , ir,     1);
			tripletListP.emplace_back(ir + 1 , ir + 1, 1);
			tripletListP.emplace_back(ir + 2 , ir + 2, 1);
			tripletListP.emplace_back(ir + 3 , ir + 3, 1);
			tripletListP.emplace_back(ir + 4 , ir + 4, 1);
			tripletListP.emplace_back(ir + 5 , ir + 5, 1);

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
			std::cout << "h_x.size(): " << h_x.size() << std::endl;
			std::cout << "AtPA=AtPB SOLVED" << std::endl;

			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}

			if(h_x.size() == 6 * m_poses.size()){
				int counter = 0;

				for(size_t i = 0; i < m_poses.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(m_poses[i]);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];
					m_poses[i] = affine_matrix_from_pose_tait_bryan(pose);
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
	std::cout << "n: add noise to poses" << std::endl;
	std::cout << "t: optimize (Tait-Bryan)" << std::endl;
	std::cout << "e: add edge (loop closure)" << std::endl;
}

void calculate_cov(std::vector<Eigen::Affine3d> m_poses,
				   std::vector<Eigen::Affine3d> m_poses_desired,
				   std::vector<std::pair<int, int>> odo_edges,
				   Eigen::MatrixXd & cov_b)
{
	Eigen::MatrixXd d2sum_dbeta2(m_poses.size() * 6, m_poses.size() * 6);
	d2sum_dbeta2 = Eigen::MatrixXd::Zero(m_poses.size() * 6, m_poses.size() * 6);

	std::vector<TaitBryanPose> poses;
	std::vector<TaitBryanPose> poses_desired;

	for(size_t i = 0 ; i < m_poses.size(); i++){
		poses.push_back(pose_tait_bryan_from_affine_matrix(m_poses[i]));
	}
	for(size_t i = 0 ; i < m_poses_desired.size(); i++){
		poses_desired.push_back(pose_tait_bryan_from_affine_matrix(m_poses_desired[i]));
	}

	for (size_t i = 0; i < odo_edges.size(); i++ ){
		Eigen::Matrix<double, 6, 1> relative_pose_measurement_odo;
		relative_pose_tait_bryan_wc_case1(relative_pose_measurement_odo,
			poses_desired[odo_edges[i].first].px,
			poses_desired[odo_edges[i].first].py,
			poses_desired[odo_edges[i].first].pz,
			poses_desired[odo_edges[i].first].om,
			poses_desired[odo_edges[i].first].fi,
			poses_desired[odo_edges[i].first].ka,
			poses_desired[odo_edges[i].second].px,
			poses_desired[odo_edges[i].second].py,
			poses_desired[odo_edges[i].second].pz,
			poses_desired[odo_edges[i].second].om,
			poses_desired[odo_edges[i].second].fi,
			poses_desired[odo_edges[i].second].ka);

		Eigen::Matrix<double, 12, 12, Eigen::RowMajor> d2sum_dbeta2i;
		relative_pose_obs_eq_tait_bryan_wc_case1_d2sum_dbeta2(
				d2sum_dbeta2i,
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

		int raw = odo_edges[i].first * 6;
		int cal = odo_edges[i].first * 6;
		for(size_t r = 0; r < 6; r++){
			for(size_t c = 0; c < 6; c++){
				d2sum_dbeta2(raw + r, cal + c) =
						d2sum_dbeta2(raw + r, cal + c) + d2sum_dbeta2i(r,c);
			}
		}
		raw = odo_edges[i].second * 6;
		cal = odo_edges[i].second * 6;
		for(size_t r = 0; r < 6; r++){
			for(size_t c = 0; c < 6; c++){
				d2sum_dbeta2(raw + r, cal + c) =
						d2sum_dbeta2(raw + r, cal + c) + d2sum_dbeta2i(r+6,c+6);
			}
		}

		raw = odo_edges[i].first * 6;
		cal = odo_edges[i].second * 6;
		for(size_t r = 0; r < 6; r++){
			for(size_t c = 0; c < 6; c++){
				d2sum_dbeta2(raw + r, cal + c) =
						d2sum_dbeta2(raw + r, cal + c) + d2sum_dbeta2i(r,c+6);
			}
		}
		raw = odo_edges[i].second * 6;
		cal = odo_edges[i].first * 6;
		for(size_t r = 0; r < 6; r++){
			for(size_t c = 0; c < 6; c++){
				d2sum_dbeta2(raw + r, cal + c) =
						d2sum_dbeta2(raw + r, cal + c) + d2sum_dbeta2i(r + 6,c);
			}
		}
	}

	for (size_t i = 0; i < 6; i++ ){
		d2sum_dbeta2(i,i) += 1000000;
	}

	Eigen::MatrixXd d2sum_dbetadx(m_poses.size() * 12, 6 * odo_edges.size());
	d2sum_dbetadx = Eigen::MatrixXd::Zero(m_poses.size() * 12, 6 * odo_edges.size());

	for (int i = 0; i < odo_edges.size() ; i++)
	{
		Eigen::Matrix<double, 6, 1> relative_pose_measurement_odo;
		relative_pose_tait_bryan_wc_case1(
			relative_pose_measurement_odo,
			poses_desired[odo_edges[i].first].px,
			poses_desired[odo_edges[i].first].py,
			poses_desired[odo_edges[i].first].pz,
			poses_desired[odo_edges[i].first].om,
			poses_desired[odo_edges[i].first].fi,
			poses_desired[odo_edges[i].first].ka,
			poses_desired[odo_edges[i].second].px,
			poses_desired[odo_edges[i].second].py,
			poses_desired[odo_edges[i].second].pz,
			poses_desired[odo_edges[i].second].om,
			poses_desired[odo_edges[i].second].fi,
			poses_desired[odo_edges[i].second].ka);
		Eigen::Matrix<double, 12, 6, Eigen::RowMajor> d2sum_dbetadx_temp;
		relative_pose_obs_eq_tait_bryan_wc_case1_d2sum_dbetadx(
			d2sum_dbetadx_temp,
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

		int raw = odo_edges[i].first * 6;
		int cal = i * 6;
		for(size_t r = 0; r < 6; r++){
			for(size_t c = 0; c < 6; c++){
				d2sum_dbetadx(raw + r, cal + c) = d2sum_dbetadx_temp(r,c);
			}
		}
		raw = odo_edges[i].second * 6;
		for(size_t r = 0; r < 6; r++){
			for(size_t c = 0; c < 6; c++){
				d2sum_dbetadx(raw + r, cal + c) = d2sum_dbetadx_temp(r + 6,c);
			}
		}
	}

	Eigen::MatrixXd cov_x(6 * odo_edges.size(), 6 * odo_edges.size());
	cov_x = Eigen::MatrixXd::Zero(6 * odo_edges.size(), 6 * odo_edges.size());

	for(int i = 0 ; i < 6 * odo_edges.size(); i+=6){
		cov_x(i,i)     = 0.005 * 0.005;
		cov_x(i+1,i+1) = 0.05 * 0.05;
		cov_x(i+2,i+2) = 0.000000001 * 0.000000001;
		cov_x(i+3,i+3) = 0.000000001 * 0.000000001;
		cov_x(i+4,i+4) = 0.000000001 * 0.000000001;
		cov_x(i+5,i+5) = 0.01 * 0.01;
	}

	cov_b = d2sum_dbeta2.inverse() * d2sum_dbetadx * cov_x * d2sum_dbetadx.transpose() * d2sum_dbeta2.inverse();
}






