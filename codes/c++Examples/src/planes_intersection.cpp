#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "distance_point_to_plane_jacobian.h"

struct Plane
{
	double a;
	double b;
	double c;
	double d;
};

const unsigned int window_width = 800;
const unsigned int window_height = 600;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -20.0;
float translate_x, translate_y = 0.0;

std::vector<Eigen::Affine3d> sheaf_of_planes;

Eigen::Affine3d m_imu = Eigen::Affine3d::Identity();
Eigen::Affine3d m_imu_volf = Eigen::Affine3d::Identity();

// std::pair<Eigen::Vector3d, Eigen::Vector3d> line;

Eigen::Vector3d intersection(0, 0, 0);

// IMU3 0.8440 -0.0259 -0.9710
// IMU4 0.8452  0.0378 -0.9717
// IMU5 0.8458  0.0376 -1.0700

// 100 0.8563 -0.0382 -0.9612
// 101 0.9838 -0.0404 -0.9605
// 102 0.8438 -0.0373 -0.9710

// IMU1  0.8569 -0.0264 -0.9612
// IMU2  0.8589  0.0499 -0.9614
// IMU10 0.9875  0.0503 -0.9608

//IMU20 0.9824 -0.0408 -0.9661 
//IMU21 0.8464 -0.0417 -0.9683 
//IMU22 0.9022 -0.0401 -1.0360 

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

int main(int argc, char *argv[])
{
	/*for (size_t i = 0; i < 10; i++)
	{
		TaitBryanPose pose;

		pose.px = random(-0.01, 0.01) + 5;
		pose.py = random(-0.01, 0.01) + 5;
		pose.pz = random(-0.01, 0.01) + 5;

		pose.om = random(-90.0, 90.0);
		pose.fi = random(-0.0001, 0.0001);
		pose.ka = random(-0.0001, 0.0001);

		sheaf_of_planes.push_back(affine_matrix_from_pose_tait_bryan(pose));
	}*/
	////////////////////////////////////////////////////////////////
	Eigen::Vector3d IMU3(0.8440, -0.0259, -0.9710);
	Eigen::Vector3d IMU4(0.8452, 0.0378, -0.9717);
	Eigen::Vector3d IMU5(0.8458, 0.0376, -1.0700);

	Eigen::Vector3d a1 = IMU4 - IMU3;
	Eigen::Vector3d a2 = IMU4 - IMU5;
	Eigen::Vector3d a3 = a1.cross(a2);
	Eigen::Vector3d a4 = a1.cross(a3);
	a1 /= a1.norm();
	a3 /= a3.norm();
	a4 /= a4.norm();

	// a1,a3,a4
	Eigen::Affine3d pose1 = Eigen::Affine3d::Identity();
	pose1.translation() = IMU4;
	pose1(0, 0) = a1.x();
	pose1(1, 0) = a1.y();
	pose1(2, 0) = a1.z();

	pose1(0, 1) = a4.x();
	pose1(1, 1) = a4.y();
	pose1(2, 1) = a4.z();

	pose1(0, 2) = a3.x();
	pose1(1, 2) = a3.y();
	pose1(2, 2) = a3.z();

	sheaf_of_planes.push_back(pose1);

	////////////////////////////////////////////////////////////////////
	// Eigen::Vector3d v100(0.8563, -0.0382, -0.9612);
	// Eigen::Vector3d v101(0.9838, -0.0404, -0.9605);
	// Eigen::Vector3d v102(0.8438, -0.0373, -0.9710);
	Eigen::Vector3d v100(0.9824, -0.0408, -0.9661);
	Eigen::Vector3d v101(0.8464, -0.0417, -0.9683);
	Eigen::Vector3d v102(0.9022, -0.0401, -1.0360);

	a1 = v100 - v101;
	a2 = v100 - v102;
	a3 = a1.cross(a2);
	a4 = a1.cross(a3);
	a1 /= a1.norm();
	a3 /= a3.norm();
	a4 /= a4.norm();

	pose1 = Eigen::Affine3d::Identity();
	pose1.translation() = v100;
	pose1(0, 0) = a1.x();
	pose1(1, 0) = a1.y();
	pose1(2, 0) = a1.z();

	pose1(0, 1) = a4.x();
	pose1(1, 1) = a4.y();
	pose1(2, 1) = a4.z();

	pose1(0, 2) = a3.x();
	pose1(1, 2) = a3.y();
	pose1(2, 2) = a3.z();

	sheaf_of_planes.push_back(pose1);

	/////////////////////////////////////////////////////////////////////////
	Eigen::Vector3d IMU1(0.8569, -0.0264, -0.9612);
	Eigen::Vector3d IMU2(0.8589, 0.0499, -0.9614);
	Eigen::Vector3d IMU10(0.9875, 0.0503, -0.9608);

	a1 = IMU1 - IMU2;
	a2 = IMU1 - IMU10;
	a3 = a1.cross(a2);
	a4 = a1.cross(a3);
	a1 /= a1.norm();
	a3 /= a3.norm();
	a4 /= a4.norm();

	// a1,a3,a4
	pose1 = Eigen::Affine3d::Identity();
	pose1.translation() = IMU1;
	pose1(0, 0) = a1.x();
	pose1(1, 0) = a1.y();
	pose1(2, 0) = a1.z();

	pose1(0, 1) = a4.x();
	pose1(1, 1) = a4.y();
	pose1(2, 1) = a4.z();

	pose1(0, 2) = a3.x();
	pose1(1, 2) = a3.y();
	pose1(2, 2) = a3.z();

	sheaf_of_planes.push_back(pose1);

	if (false == initGL(&argc, argv))
	{
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

bool initGL(int *argc, char **argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow("planes_intersection");
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
	gluPerspective(60.0, (GLfloat)window_width / (GLfloat)window_height, 0.01,
				   10000.0);
	glutReshapeFunc(reshape);

	return true;
}

void display()
{
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

	for (size_t i = 0; i < sheaf_of_planes.size(); i++)
	{
		Eigen::Affine3d &m = sheaf_of_planes[i];

		glBegin(GL_LINES);
		glColor3f(m(0, 2), m(1, 2), m(2, 2));
		for (float step = -1.0f; step <= 1.1f; step += 0.1f)
		{
			Eigen::Vector3d v1(-1, step, 0);
			Eigen::Vector3d v2(1, step, 0);

			Eigen::Vector3d v1t = m * v1;
			Eigen::Vector3d v2t = m * v2;

			glVertex3f(v1t.x(), v1t.y(), v1t.z());
			glVertex3f(v2t.x(), v2t.y(), v2t.z());

			v1 = Eigen::Vector3d(step, -1, 0);
			v2 = Eigen::Vector3d(step, 1, 0);

			v1t = m * v1;
			v2t = m * v2;

			glVertex3f(v1t.x(), v1t.y(), v1t.z());
			glVertex3f(v2t.x(), v2t.y(), v2t.z());
		}
		glEnd();
	}

	glColor3f(0, 0.3, 0);
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3f(intersection.x(), intersection.y(), intersection.z());
	glEnd();
	glPointSize(1);

	//m_imu
	glLineWidth(3);
	glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(m_imu(0, 3), m_imu(1, 3), m_imu(2, 3));
		glVertex3f(m_imu(0, 3) + m_imu(0, 0), m_imu(1, 3) + m_imu(1, 0), m_imu(2, 3) + m_imu(2, 0));

		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(m_imu(0, 3), m_imu(1, 3), m_imu(2, 3));
		glVertex3f(m_imu(0, 3) + m_imu(0, 1), m_imu(1, 3) + m_imu(1, 1), m_imu(2, 3) + m_imu(2, 1));

		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(m_imu(0, 3), m_imu(1, 3), m_imu(2, 3));
		glVertex3f(m_imu(0, 3) + m_imu(0, 2), m_imu(1, 3) + m_imu(1, 2), m_imu(2, 3) + m_imu(2, 2));
	glEnd();
	glLineWidth(1);


	//volfram
	glLineWidth(3);
	glBegin(GL_LINES);
		glColor3f(1.0f, 0.5f, 0.5f);
		glVertex3f(m_imu_volf(0, 3), m_imu_volf(1, 3), m_imu_volf(2, 3));
		glVertex3f(m_imu_volf(0, 3) + m_imu_volf(0, 0), m_imu_volf(1, 3) + m_imu_volf(1, 0), m_imu_volf(2, 3) + m_imu_volf(2, 0));

		glColor3f(5.0f, 1.0f, 0.5f);
		glVertex3f(m_imu_volf(0, 3), m_imu_volf(1, 3), m_imu_volf(2, 3));
		glVertex3f(m_imu_volf(0, 3) + m_imu_volf(0, 1), m_imu_volf(1, 3) + m_imu_volf(1, 1), m_imu_volf(2, 3) + m_imu_volf(2, 1));

		glColor3f(0.5f, 0.5f, 1.0f);
		glVertex3f(m_imu_volf(0, 3), m_imu_volf(1, 3), m_imu_volf(2, 3));
		glVertex3f(m_imu_volf(0, 3) + m_imu_volf(0, 2), m_imu_volf(1, 3) + m_imu_volf(1, 2), m_imu_volf(2, 3) + m_imu_volf(2, 2));
	glEnd();
	glLineWidth(1);
	

	glutSwapBuffers();
}

void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
	switch (key)
	{
	case (27):
	{
		glutDestroyWindow(glutGetWindow());
		return;
	}
	case 't':
	{
		intersection.x() += ((rand() % 1000000000000 / 1000000000000.0) - 0.5) * 0.000000000001;
		intersection.y() += ((rand() % 1000000000000 / 1000000000000.0) - 0.5) * 0.000000000001;
		intersection.z() += ((rand() % 1000000000000 / 1000000000000.0) - 0.5) * 0.000000000001;

		std::vector<Eigen::Triplet<double>> tripletListA;
		std::vector<Eigen::Triplet<double>> tripletListP;
		std::vector<Eigen::Triplet<double>> tripletListB;

		for (size_t i = 0; i < sheaf_of_planes.size(); i++)
		{
			Plane plane;
			plane.a = sheaf_of_planes[i](0, 2);
			plane.b = sheaf_of_planes[i](1, 2);
			plane.c = sheaf_of_planes[i](2, 2);
			plane.d = -plane.a * sheaf_of_planes[i](0, 3) - plane.b * sheaf_of_planes[i](1, 3) - plane.c * sheaf_of_planes[i](2, 3);

			Eigen::Matrix<double, 1, 1> delta;
			delta_distance_point_to_plane(delta, intersection.x(), intersection.y(), intersection.z(), plane.a, plane.b, plane.c, plane.d);

			Eigen::Matrix<double, 1, 3> jacobian;
			delta_distance_point_to_plane_jacobian(jacobian, intersection.x(), intersection.y(), intersection.z(), plane.a, plane.b, plane.c, plane.d);

			int ir = tripletListB.size();

			for (size_t j = 0; j < 1; j++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					tripletListA.emplace_back(ir + j, k, -jacobian(j, k));
				}
			}

			tripletListP.emplace_back(ir, ir, 1);
			// tripletListP.emplace_back(ir + 1, ir + 1, 1);
			// tripletListP.emplace_back(ir + 2, ir + 2, 1);

			tripletListB.emplace_back(ir, 0, delta(0, 0));
			// tripletListB.emplace_back(ir, 1, delta(1, 0));
			// tripletListB.emplace_back(ir, 2, delta(2, 0));
		}

		Eigen::SparseMatrix<double> matA(tripletListB.size(), 3);
		Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
		Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

		matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
		matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
		matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

		Eigen::SparseMatrix<double> AtPA(3, 3);
		Eigen::SparseMatrix<double> AtPB(3, 1);

		{
			Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
			AtPA = (AtP)*matA;
			AtPB = (AtP)*matB;
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

		for (int k = 0; k < x.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(x, k); it; ++it)
			{
				h_x.push_back(it.value());
				// std::cout << it.value() << std::endl;
			}
		}

		std::cout << "h_x.size(): " << h_x.size() << std::endl;
		for (size_t i = 0; i < h_x.size(); i++)
		{
			std::cout << h_x[i] << std::endl;
		}

		if (h_x.size() == 3)
		{
			std::cout << "AtPA=AtPB SOLVED" << std::endl;
			std::cout << "update" << std::endl;

			intersection.x() += h_x[0];
			intersection.y() += h_x[1];
			intersection.z() += h_x[2];

			std::cout << "intersection: " << intersection.x() << " " << intersection.y() << " " << intersection.z() << std::endl;

			for (size_t i = 0; i < sheaf_of_planes.size(); i++)
			{
				Plane plane;
				plane.a = sheaf_of_planes[i](0, 2);
				plane.b = sheaf_of_planes[i](1, 2);
				plane.c = sheaf_of_planes[i](2, 2);
				plane.d = -plane.a * sheaf_of_planes[i](0, 3) - plane.b * sheaf_of_planes[i](1, 3) - plane.c * sheaf_of_planes[i](2, 3);

				Eigen::Matrix<double, 1, 1> delta;
				delta_distance_point_to_plane(delta, intersection.x(), intersection.y(), intersection.z(), plane.a, plane.b, plane.c, plane.d);

				std::cout << "distance: " << delta << std::endl;
			}

			// orthogonalize_rotation(Eigen::Affine3d & m)

			Eigen::Affine3d m = Eigen::Affine3d::Identity();
			m.translation() = intersection;
			m(0, 0) = sheaf_of_planes[0](0, 2);
			m(1, 0) = sheaf_of_planes[0](1, 2);
			m(2, 0) = sheaf_of_planes[0](2, 2);

			m(0, 1) = sheaf_of_planes[1](0, 2);
			m(1, 1) = sheaf_of_planes[1](1, 2);
			m(2, 1) = sheaf_of_planes[1](2, 2);

			m(0, 2) = -sheaf_of_planes[2](0, 2);
			m(1, 2) = -sheaf_of_planes[2](1, 2);
			m(2, 2) = -sheaf_of_planes[2](2, 2);

			std::cout << "m before" << std::endl;
			std::cout << m.matrix() << std::endl;

			{
			m_imu = m;

			// check angles between vetors
			Eigen::Vector3d vx(m_imu(0, 0), m_imu(1, 0), m_imu(2, 0));
			Eigen::Vector3d vy(m_imu(0, 1), m_imu(1, 1), m_imu(2, 1));
			Eigen::Vector3d vz(m_imu(0, 2), m_imu(1, 2), m_imu(2, 2));

			std::cout << "lengths vx: " << vx.norm() << " vy: " << vy.norm() << " vz: " << vz.norm() << std::endl;
			double a1 = acos(vx.dot(vy) / (vx.norm() * vy.norm()));
			std::cout << "angle_deg xy: " << a1 * 180.0 / M_PI << std::endl;
			double a2 = acos(vy.dot(vz) / (vy.norm() * vz.norm()));
			std::cout << "angle_deg yz: " << a2 * 180.0 / M_PI << std::endl;
			double a3 = acos(vx.dot(vz) / (vx.norm() * vz.norm()));
			std::cout << "angle_deg xz: " << a3 * 180.0 / M_PI << std::endl;
			}


			orthogonalize_rotation(m);

			m(0, 0) *= -1.0;
			m(1, 0) *= -1.0;
			m(2, 0) *= -1.0;

			m(0, 1) *= -1.0;
			m(1, 1) *= -1.0;
			m(2, 1) *= -1.0;

			//m(0, 2) *= -1.0;
			//m(1, 2) *= -1.0;
			//m(2, 2) *= -1.0;

			std::cout << "m after" << std::endl;
			std::cout << m.matrix() << std::endl;

			m_imu = m;

			//check angles between vetors
			Eigen::Vector3d vx(m_imu(0, 0), m_imu(1, 0), m_imu(2, 0));
			Eigen::Vector3d vy(m_imu(0, 1), m_imu(1, 1), m_imu(2, 1));
			Eigen::Vector3d vz(m_imu(0, 2), m_imu(1, 2), m_imu(2, 2));

			std::cout << "lengths vx: " << vx.norm() << " vy: " << vy.norm() << " vz: " << vz.norm() << std::endl;
			double a1 = acos(vx.dot(vy) / (vx.norm() * vy.norm()));
			std::cout << "angle_deg xy: " << a1*180.0/M_PI << std::endl;
			double a2 = acos(vy.dot(vz) / (vy.norm() * vz.norm()));
			std::cout << "angle_deg yz: " << a2 * 180.0 / M_PI << std::endl;
			double a3 = acos(vx.dot(vz) / (vx.norm() * vz.norm()));
			std::cout << "angle_deg xz: " << a3 * 180.0 / M_PI << std::endl;

			//(0.999965 | 0.00690769 | -0.00467485 0.00692001 | -0.999973 | 0.00262451 0.00465659 | 0.00265677 | 0.999986)

			{
			m_imu = m;

			// check angles between vetors
			Eigen::Vector3d vx(0.999965, 0.00690769, -0.00467485);
			Eigen::Vector3d vy(0.00692001, -0.999973, 0.00262451);
			Eigen::Vector3d vz(0.00465659, 0.00265677, 0.999986);

			std::cout << "lengths vx: " << vx.norm() << " vy: " << vy.norm() << " vz: " << vz.norm() << std::endl;
			double a1 = acos(vx.dot(vy) / (vx.norm() * vy.norm()));
			std::cout << "angle_deg xy: " << a1 * 180.0 / M_PI << std::endl;
			double a2 = acos(vy.dot(vz) / (vy.norm() * vz.norm()));
			std::cout << "angle_deg yz: " << a2 * 180.0 / M_PI << std::endl;
			double a3 = acos(vx.dot(vz) / (vx.norm() * vz.norm()));
			std::cout << "angle_deg xz: " << a3 * 180.0 / M_PI << std::endl;

			m_imu_volf = m;
			m_imu_volf(0, 0) = vx.x();
			m_imu_volf(1, 0) = vx.y();
			m_imu_volf(2, 0) = vx.z();

			m_imu_volf(0, 1) = vy.x();
			m_imu_volf(1, 1) = vy.y();
			m_imu_volf(2, 1) = vy.z();

			m_imu_volf(0, 2) = vz.x();
			m_imu_volf(1, 2) = vz.y();
			m_imu_volf(2, 2) = vz.z();
			}
		}
		else
		{
			std::cout << "AtPA=AtPB FAILED" << std::endl;
		}

		break;
	}
	case 'n':
	{
		/*
		TaitBryanPose pose;

		pose.px = random(-2.0, 2.0);
		pose.py = random(-2.0, 2.0);
		pose.pz = random(-2.0, 2.0);

		pose.om = random(-2.0, 2.0);
		pose.fi = random(-2.0, 2.0);
		pose.ka = random(-2.0, 2.0);

		for (size_t i = 0; i < sheaf_of_planes.size(); i++)
		{
			TaitBryanPose pose_omfika;

			pose_omfika.px = 0;
			pose_omfika.py = 0;
			pose_omfika.pz = 0;

			pose_omfika.om = random(-0.01, 0.01);
			pose_omfika.fi = random(-0.01, 0.01);
			pose_omfika.ka = random(-0.01, 0.01);

			sheaf_of_planes[i] = affine_matrix_from_pose_tait_bryan(pose_omfika) * affine_matrix_from_pose_tait_bryan(pose) * sheaf_of_planes[i] * affine_matrix_from_pose_tait_bryan(pose_omfika).inverse();
		}
		*/
		break;
	}
	}

	printHelp();
	glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
		mouse_buttons |= 1 << button;
	}
	else if (state == GLUT_UP)
	{
		mouse_buttons = 0;
	}

	mouse_old_x = x;
	mouse_old_y = y;
}

void motion(int x, int y)
{
	float dx, dy;
	dx = (float)(x - mouse_old_x);
	dy = (float)(y - mouse_old_y);

	if (mouse_buttons & 1)
	{
		rotate_x += dy * 0.2f;
		rotate_y += dx * 0.2f;
	}
	else if (mouse_buttons & 4)
	{
		translate_z += dy * 0.05f;
	}
	else if (mouse_buttons & 3)
	{
		translate_x += dx * 0.05f;
		translate_y -= dy * 0.05f;
	}

	mouse_old_x = x;
	mouse_old_y = y;

	glutPostRedisplay();
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (GLfloat)w / (GLfloat)h, 0.01, 10000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void printHelp()
{
	std::cout << "-------help-------" << std::endl;
	// std::cout << "n: modify planes" << std::endl;
	std::cout << "t: optimize" << std::endl;
}
