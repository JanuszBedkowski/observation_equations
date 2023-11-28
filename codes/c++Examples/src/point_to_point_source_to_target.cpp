#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"
#include "point_to_point_source_to_target_tait_bryan_wc_jacobian.h"
#include "point_to_point_source_to_target_rodrigues_wc_jacobian.h"
#include "point_to_point_source_to_target_quaternion_wc_jacobian.h"
#include "quaternion_constraint_jacobian.h"

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

Eigen::Affine3d pose_source;
Eigen::Affine3d m_imu;
std::vector<Eigen::Vector3d> points_target_global;
std::vector<Eigen::Vector3d> points_source_local;

int main(int argc, char *argv[])
{

	// for(size_t i = 0 ; i < 100; i++){
	//	Eigen::Vector3d p;
	//	p.x() = random(-100.0, 100.0);
	//	p.y() = random(-100.0, 100.0);
	//	p.z() = random(-100.0, 100.0);
	//	points_target_global.push_back(p);
	// }

	points_target_global.emplace_back(-2.5974, 2.9859, 0.4050);
	points_target_global.emplace_back(-4.4065, 9.4782, -0.7643);
	points_target_global.emplace_back(0.0281, 9.4712, 1.3241);
	points_target_global.emplace_back(2.5421, 2.6989, 0.5996);
	points_target_global.emplace_back(7.2702, -0.7880, 0.1246);
	points_target_global.emplace_back(2.1284, -3.5989, 1.1592);

	points_target_global.emplace_back(-2.4885, -3.5796, 0.7323);
	points_target_global.emplace_back(-7.9781, -3.4132, 1.1302);
	points_target_global.emplace_back(-6.0474, -2.2724, -1.0379);
	points_target_global.emplace_back(-3.7999, -0.0691, -1.2115);
	points_target_global.emplace_back(-7.9691, 1.5899, 1.4859);
	points_target_global.emplace_back(-7.8522, 2.9080, -1.1242);
	/*
	T1	-2.5974	2.9859	0.4050
	T2	-4.4065	9.4782	-0.7643
	T3	0.0281	9.4712	1.3241
	T4	2.5421	2.6989	0.5996
	T5	7.2702	-0.7880	0.1246
	T6	2.1284	-3.5989	1.1592

	T7	-2.4885	-3.5796	0.7323
	T8	-7.9781	-3.4132	1.1302
	T9	-6.0474	-2.2724	-1.0379
	T10	-3.7999	-0.0691	-1.2115
	T11	-7.9691	1.5899	1.4859
	T12	-7.8522	2.9080	-1.1242

	IMU1	-0.0264	0.8569	-0.9612
	IMU2	0.0499	0.8589	-0.9614
	IMU3	-0.0259	0.8440	-0.9710
	IMU4	0.0378	0.8452	-0.9717
	IMU5	0.0376	0.8458	-1.0700
	IMU6	0.0296	0.9937	-1.0143
	IMU7	-0.0226	0.9929	-1.0091
	IMU8	0.0084	0.9934	-1.0763
	IMU9 	-0.0324 0.9842 -0.9605
	IMU10 	0.0503 0.9875 -0.9608
	IMU20 -0.0408  0.9824 -0.9661
	IMU21 -0.0417 0.8464 -0.9683
	IMU22 -0.0401 0.9022 -1.0360
	zero_tachimetr 0 0 0*/

	points_source_local.emplace_back(-1.9274, 2.6313, 0.9558);
	points_source_local.emplace_back(-8.4013, 4.5254, -0.1893);
	points_source_local.emplace_back(-8.4443, 0.0885, 1.8957);
	points_source_local.emplace_back(-1.7062, -2.5128, 1.1451);
	points_source_local.emplace_back(1.7165, -7.2857, 0.6531);
	points_source_local.emplace_back(4.6001, -2.1822, 1.6830);

	points_source_local.emplace_back(4.6389, 2.4359, 1.2592);
	points_source_local.emplace_back(4.5459, 7.9286, 1.6627);
	points_source_local.emplace_back(3.3712, 6.0130, -0.5035);
	points_source_local.emplace_back(1.1376, 3.7936, -0.6709);
	points_source_local.emplace_back(-0.4571, 7.9836, 2.0364);
	points_source_local.emplace_back(-1.7863, 7.8864, -0.5695);

	for (auto &p : points_source_local)
	{
		p.y() *= -1.0;
	}

	/*
	T1 -1.9274 2.6313 0.9558
	T2 -8.4013 4.5254 -0.1893
	T3 -8.4443 0.0885 1.8957
	T4 -1.7062 -2.5128 1.1451
	T5 1.7165 -7.2857 0.6531
	T6 4.6001 -2.1822 1.6830

	T7 4.6389 2.4359 1.2592
	T8 4.5459 7.9286 1.6627
	T9 3.3712 6.0130 -0.5035
	T10 1.1376 3.7936 -0.6709
	T11 -0.4571 7.9836 2.0364
	T12 -1.7863 7.8864 -0.5695
	zero_Imager 0 0 0
	*/

	TaitBryanPose pose;
	// pose.px = -4;
	// pose.py = 0.4;
	// pose.pz = -0.2;
	// pose.om = 0.1;
	// pose.fi = 0.2;
	// pose.ka = 0.3;
	pose.px = 0;
	pose.py = 0;
	pose.pz = 0;
	pose.om = 0;
	pose.fi = 0;
	pose.ka = 0;
	pose_source = affine_matrix_from_pose_tait_bryan(pose);

	/*pose_source(0, 0) = -0.013042502105;
	pose_source(0, 1) = -0.999914586544;
	pose_source(0, 2) = -0.000834811886;
	pose_source(0, 3) = 0.008754968643;
	pose_source(1, 0) = -0.999908566475;
	pose_source(1, 1) = 0.013039439917;
	pose_source(1, 2) = 0.003574120346;
	pose_source(1, 3) = 1.021424531937;
	pose_source(2, 0) = 0.003562929574;
	pose_source(2, 1) = -0.000881351007;
	pose_source(2, 2) = 0.999993264675;
	pose_source(2, 3) = -0.541456818581;
	pose_source(3, 0) = 0.000000000000;
	pose_source(3, 1) = 0.000000000000;
	pose_source(3, 2) = 0.000000000000;
	pose_source(3, 3) = 1.000000000000;*/

	// Eigen::Vector3d vx(pose_source(0, 0), pose_source(1, 0), pose_source(2, 0));
	// Eigen::Vector3d vy(pose_source(0, 1), pose_source(1, 1), pose_source(2, 1));
	// Eigen::Vector3d vz = vx.cross(vy);

	// std::cout << "vz " << vz << std::endl;

	// Eigen::Affine3d m_inv = pose_source.inverse();
	// for(size_t j = 0 ; j < points_target_global.size(); j++){
	//	Eigen::Vector3d vt = m_inv * points_target_global[j];
	//	points_source_local.push_back(vt);
	// }

	m_imu(0, 0) = 0.999967;
	m_imu(0, 1) = -0.00593569;
	m_imu(0, 2) = 0.00551016;

	// m_imu(0, 3) = 0.843641;
	m_imu(0, 3) = -0.0418461;

	m_imu(1, 0) = -0.00597756;
	m_imu(1, 1) = -0.999953;
	m_imu(1, 2) = 0.0076135;

	// m_imu(1, 3) = -0.0418461;
	m_imu(1, 3) = 0.843641;

	m_imu(2, 0) = 0.00546471;
	m_imu(2, 1) = -0.00764619;
	m_imu(2, 2) = -0.999956;
	m_imu(2, 3) = -0.96122;
	m_imu(3, 0) = 0;
	m_imu(3, 1) = 0;
	m_imu(3, 2) = 0;
	m_imu(3, 3) = 1;

	auto mm = m_imu;
	mm(0, 0) = m_imu(0, 1);
	mm(1, 0) = m_imu(1, 1);
	mm(2, 0) = m_imu(2, 1);

	mm(0, 1) = -m_imu(0, 0);
	mm(1, 1) = -m_imu(1, 0);
	mm(2, 1) = -m_imu(2, 0);

	Eigen::Affine3d m2 = Eigen::Affine3d::Identity();
	m2(0, 0) = -1;
	m2(1, 1) = 1;
	m2(2, 2) = -1;

	m_imu = mm;

	Eigen::Vector3d imu_offset(0.058, 0.046, 0.071);

	m_imu = m_imu * m2;

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
	glutCreateWindow("point to point source to target");
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
	gluPerspective(60.0, (GLfloat)window_width / (GLfloat)window_height, 0.01, 10000.0);
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

	Eigen::Affine3d m = pose_source;

	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(m(0, 3), m(1, 3), m(2, 3));
	glVertex3f(m(0, 3) + m(0, 0), m(1, 3) + m(1, 0), m(2, 3) + m(2, 0));

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(m(0, 3), m(1, 3), m(2, 3));
	glVertex3f(m(0, 3) + m(0, 1), m(1, 3) + m(1, 1), m(2, 3) + m(2, 1));

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(m(0, 3), m(1, 3), m(2, 3));
	glVertex3f(m(0, 3) + m(0, 2), m(1, 3) + m(1, 2), m(2, 3) + m(2, 2));
	glEnd();

	glColor3f(1, 0, 0);
	glPointSize(3);
	glBegin(GL_POINTS);
	for (const auto &p : points_source_local)
	{
		Eigen::Vector3d vt;
		vt = m * p;
		glVertex3f(vt.x(), vt.y(), vt.z());
	}
	glEnd();
	glPointSize(1);

	glColor3f(0, 0, 1);
	glPointSize(3);
	glBegin(GL_POINTS);

	for (auto &p : points_target_global)
	{
		glVertex3f(p.x(), p.y(), p.z());
	}
	glEnd();
	glPointSize(1);

	glColor3f(0, 1, 0);
	glBegin(GL_LINES);
	for (int i = 0; i < points_source_local.size(); i++)
	{
		Eigen::Vector3d v(points_source_local[i].x(), points_source_local[i].y(), points_source_local[i].z());
		Eigen::Vector3d vt;
		vt = m * v;
		glVertex3f(vt.x(), vt.y(), vt.z());
		glVertex3f(points_target_global[i].x(), points_target_global[i].y(), points_target_global[i].z());
	}
	glEnd();

	m = m_imu;
	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(m(0, 3), m(1, 3), m(2, 3));
	glVertex3f(m(0, 3) + m(0, 0) * 0.1, m(1, 3) + m(1, 0) * 0.1, m(2, 3) + m(2, 0) * 0.1);

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(m(0, 3), m(1, 3), m(2, 3));
	glVertex3f(m(0, 3) + m(0, 1) * 0.1, m(1, 3) + m(1, 1) * 0.1, m(2, 3) + m(2, 1) * 0.1);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(m(0, 3), m(1, 3), m(2, 3));
	glVertex3f(m(0, 3) + m(0, 2) * 0.1, m(1, 3) + m(1, 2) * 0.1, m(2, 3) + m(2, 2) * 0.1);
	glEnd();

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
	case 'n':
	{
		TaitBryanPose pose;
		pose.px = random(-1.0, 1.0);
		pose.py = random(-1.0, 1.0);
		pose.pz = random(-1.0, 1.0);
		pose.om = random(-0.5, 0.5);
		pose.fi = random(-0.5, 0.5);
		pose.ka = random(-0.5, 0.5);

		pose_source = pose_source * affine_matrix_from_pose_tait_bryan(pose);
		break;
	}
	case 't':
	{
		std::vector<Eigen::Triplet<double>> tripletListA;
		std::vector<Eigen::Triplet<double>> tripletListP;
		std::vector<Eigen::Triplet<double>> tripletListB;

		TaitBryanPose pose_s = pose_tait_bryan_from_affine_matrix(pose_source);

		for (size_t i = 0; i < points_source_local.size(); i++)
		{
			Eigen::Vector3d &p_t = points_target_global[i];
			Eigen::Vector3d &p_s = points_source_local[i];
			double delta_x;
			double delta_y;
			double delta_z;
			point_to_point_source_to_target_tait_bryan_wc(delta_x, delta_y, delta_z,
														  pose_s.px, pose_s.py, pose_s.pz, pose_s.om, pose_s.fi, pose_s.ka,
														  p_s.x(), p_s.y(), p_s.z(), p_t.x(), p_t.y(), p_t.z());

			std::cout << "delta_x " << delta_x << " delta_y " << delta_y << " delta_z " << delta_z << std::endl;

			Eigen::Matrix<double, 3, 6, Eigen::RowMajor>
				jacobian;
			point_to_point_source_to_target_tait_bryan_wc_jacobian(jacobian, pose_s.px, pose_s.py, pose_s.pz, pose_s.om, pose_s.fi, pose_s.ka, p_s.x(), p_s.y(), p_s.z());

			int ir = tripletListB.size();

			for (int row = 0; row < 3; row++)
			{
				for (int col = 0; col < 6; col++)
				{
					if (jacobian(row, col) != 0.0)
					{
						tripletListA.emplace_back(ir + row, col, -jacobian(row, col));
					}
				}
			}

			tripletListP.emplace_back(ir, ir, 1);
			tripletListP.emplace_back(ir + 1, ir + 1, 1);
			tripletListP.emplace_back(ir + 2, ir + 2, 1);

			tripletListB.emplace_back(ir, 0, delta_x);
			tripletListB.emplace_back(ir + 1, 0, delta_y);
			tripletListB.emplace_back(ir + 2, 0, delta_z);
		}

		Eigen::SparseMatrix<double> matA(tripletListB.size(), 6);
		Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
		Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

		matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
		matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
		matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

		Eigen::SparseMatrix<double> AtPA(6, 6);
		Eigen::SparseMatrix<double> AtPB(6, 1);

		Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
		AtPA = AtP * matA;
		AtPB = AtP * matB;

		tripletListA.clear();
		tripletListP.clear();
		tripletListB.clear();

		Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(AtPA);
		Eigen::SparseMatrix<double> x = solver.solve(AtPB);

		std::vector<double> h_x;

		for (int k = 0; k < x.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(x, k); it; ++it)
			{
				h_x.push_back(it.value());
			}
		}

		if (h_x.size() == 6)
		{
			for (size_t i = 0; i < h_x.size(); i++)
			{
				std::cout << h_x[i] << std::endl;
			}
			std::cout << "AtPA=AtPB SOLVED" << std::endl;
			// std::cout << "update" << std::endl;

			int counter = 0;

			TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(pose_source);
			pose.px += h_x[counter++];
			pose.py += h_x[counter++];
			pose.pz += h_x[counter++];
			pose.om += h_x[counter++];
			pose.fi += h_x[counter++];
			pose.ka += h_x[counter++];

			pose_source = affine_matrix_from_pose_tait_bryan(pose);

			std::cout << pose_source.matrix() << std::endl;

			std::cout << "m_imu.inv() * pose_source" << std::endl;
			std::cout << (m_imu.inverse() * pose_source).matrix() << std::endl; 
		}
		else
		{
			std::cout << "AtPA=AtPB FAILED" << std::endl;
		}
		break;
	}
	case 'r':
	{

		TaitBryanPose posetb = pose_tait_bryan_from_affine_matrix(pose_source);
		posetb.om += random(-0.000001, 0.000001);
		posetb.fi += random(-0.000001, 0.000001);
		posetb.ka += random(-0.000001, 0.000001);
		pose_source = affine_matrix_from_pose_tait_bryan(posetb);

		std::vector<Eigen::Triplet<double>> tripletListA;
		std::vector<Eigen::Triplet<double>> tripletListP;
		std::vector<Eigen::Triplet<double>> tripletListB;

		RodriguesPose pose_s = pose_rodrigues_from_affine_matrix(pose_source);

		for (size_t i = 0; i < points_source_local.size(); i++)
		{
			Eigen::Vector3d &p_t = points_target_global[i];
			Eigen::Vector3d &p_s = points_source_local[i];
			double delta_x;
			double delta_y;
			double delta_z;
			point_to_point_source_to_target_rodrigues_wc(delta_x, delta_y, delta_z, pose_s.px, pose_s.py, pose_s.pz, pose_s.sx, pose_s.sy, pose_s.sz, p_s.x(), p_s.y(), p_s.z(), p_t.x(), p_t.y(), p_t.z());

			Eigen::Matrix<double, 3, 6, Eigen::RowMajor> jacobian;
			point_to_point_source_to_target_rodrigues_wc_jacobian(jacobian, pose_s.px, pose_s.py, pose_s.pz, pose_s.sx, pose_s.sy, pose_s.sz, p_s.x(), p_s.y(), p_s.z());

			int ir = tripletListB.size();

			tripletListA.emplace_back(ir, 0, -jacobian(0, 0));
			tripletListA.emplace_back(ir, 1, -jacobian(0, 1));
			tripletListA.emplace_back(ir, 2, -jacobian(0, 2));
			tripletListA.emplace_back(ir, 3, -jacobian(0, 3));
			tripletListA.emplace_back(ir, 4, -jacobian(0, 4));
			tripletListA.emplace_back(ir, 5, -jacobian(0, 5));

			tripletListA.emplace_back(ir + 1, 0, -jacobian(1, 0));
			tripletListA.emplace_back(ir + 1, 1, -jacobian(1, 1));
			tripletListA.emplace_back(ir + 1, 2, -jacobian(1, 2));
			tripletListA.emplace_back(ir + 1, 3, -jacobian(1, 3));
			tripletListA.emplace_back(ir + 1, 4, -jacobian(1, 4));
			tripletListA.emplace_back(ir + 1, 5, -jacobian(1, 5));

			tripletListA.emplace_back(ir + 2, 0, -jacobian(2, 0));
			tripletListA.emplace_back(ir + 2, 1, -jacobian(2, 1));
			tripletListA.emplace_back(ir + 2, 2, -jacobian(2, 2));
			tripletListA.emplace_back(ir + 2, 3, -jacobian(2, 3));
			tripletListA.emplace_back(ir + 2, 4, -jacobian(2, 4));
			tripletListA.emplace_back(ir + 2, 5, -jacobian(2, 5));

			tripletListP.emplace_back(ir, ir, cauchy(delta_x, 1));
			tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta_y, 1));
			tripletListP.emplace_back(ir + 2, ir + 2, cauchy(delta_z, 1));

			tripletListB.emplace_back(ir, 0, delta_x);
			tripletListB.emplace_back(ir + 1, 0, delta_y);
			tripletListB.emplace_back(ir + 2, 0, delta_z);
		}

		Eigen::SparseMatrix<double> matA(tripletListB.size(), 6);
		Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
		Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

		matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
		matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
		matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

		Eigen::SparseMatrix<double> AtPA(6, 6);
		Eigen::SparseMatrix<double> AtPB(6, 1);

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
			}
		}

		if (h_x.size() == 6)
		{
			for (size_t i = 0; i < h_x.size(); i++)
			{
				std::cout << h_x[i] << std::endl;
			}
			std::cout << "AtPA=AtPB SOLVED" << std::endl;
			std::cout << "update" << std::endl;

			int counter = 0;

			RodriguesPose pose = pose_rodrigues_from_affine_matrix(pose_source);
			pose.px += h_x[counter++];
			pose.py += h_x[counter++];
			pose.pz += h_x[counter++];
			pose.sx += h_x[counter++];
			pose.sy += h_x[counter++];
			pose.sz += h_x[counter++];

			pose_source = affine_matrix_from_pose_rodrigues(pose);
		}
		else
		{
			std::cout << "AtPA=AtPB FAILED" << std::endl;
		}
		break;
	}
	case 'q':
	{
		std::vector<Eigen::Triplet<double>> tripletListA;
		std::vector<Eigen::Triplet<double>> tripletListP;
		std::vector<Eigen::Triplet<double>> tripletListB;

		QuaternionPose pose_s = pose_quaternion_from_affine_matrix(pose_source);

		for (size_t i = 0; i < points_source_local.size(); i++)
		{
			Eigen::Vector3d &p_t = points_target_global[i];
			Eigen::Vector3d &p_s = points_source_local[i];
			double delta_x;
			double delta_y;
			double delta_z;
			point_to_point_source_to_target_quaternion_wc(delta_x, delta_y, delta_z, pose_s.px, pose_s.py, pose_s.pz, pose_s.q0, pose_s.q1, pose_s.q2, pose_s.q3, p_s.x(), p_s.y(), p_s.z(), p_t.x(), p_t.y(), p_t.z());

			Eigen::Matrix<double, 3, 7, Eigen::RowMajor> jacobian;
			point_to_point_source_to_target_quaternion_wc_jacobian(jacobian, pose_s.px, pose_s.py, pose_s.pz, pose_s.q0, pose_s.q1, pose_s.q2, pose_s.q3, p_s.x(), p_s.y(), p_s.z());

			int ir = tripletListB.size();

			tripletListA.emplace_back(ir, 0, -jacobian(0, 0));
			tripletListA.emplace_back(ir, 1, -jacobian(0, 1));
			tripletListA.emplace_back(ir, 2, -jacobian(0, 2));
			tripletListA.emplace_back(ir, 3, -jacobian(0, 3));
			tripletListA.emplace_back(ir, 4, -jacobian(0, 4));
			tripletListA.emplace_back(ir, 5, -jacobian(0, 5));
			tripletListA.emplace_back(ir, 6, -jacobian(0, 6));

			tripletListA.emplace_back(ir + 1, 0, -jacobian(1, 0));
			tripletListA.emplace_back(ir + 1, 1, -jacobian(1, 1));
			tripletListA.emplace_back(ir + 1, 2, -jacobian(1, 2));
			tripletListA.emplace_back(ir + 1, 3, -jacobian(1, 3));
			tripletListA.emplace_back(ir + 1, 4, -jacobian(1, 4));
			tripletListA.emplace_back(ir + 1, 5, -jacobian(1, 5));
			tripletListA.emplace_back(ir + 1, 6, -jacobian(1, 6));

			tripletListA.emplace_back(ir + 2, 0, -jacobian(2, 0));
			tripletListA.emplace_back(ir + 2, 1, -jacobian(2, 1));
			tripletListA.emplace_back(ir + 2, 2, -jacobian(2, 2));
			tripletListA.emplace_back(ir + 2, 3, -jacobian(2, 3));
			tripletListA.emplace_back(ir + 2, 4, -jacobian(2, 4));
			tripletListA.emplace_back(ir + 2, 5, -jacobian(2, 5));
			tripletListA.emplace_back(ir + 2, 6, -jacobian(2, 6));

			tripletListP.emplace_back(ir, ir, cauchy(delta_x, 1));
			tripletListP.emplace_back(ir + 1, ir + 1, cauchy(delta_y, 1));
			tripletListP.emplace_back(ir + 2, ir + 2, cauchy(delta_z, 1));

			tripletListB.emplace_back(ir, 0, delta_x);
			tripletListB.emplace_back(ir + 1, 0, delta_y);
			tripletListB.emplace_back(ir + 2, 0, delta_z);
		}

		int ir = tripletListB.size();
		QuaternionPose pose = pose_quaternion_from_affine_matrix(pose_source);

		double delta;
		quaternion_constraint(delta, pose.q0, pose.q1, pose.q2, pose.q3);

		Eigen::Matrix<double, 1, 4> jacobian;
		quaternion_constraint_jacobian(jacobian, pose.q0, pose.q1, pose.q2, pose.q3);

		tripletListA.emplace_back(ir, 3, -jacobian(0, 0));
		tripletListA.emplace_back(ir, 4, -jacobian(0, 1));
		tripletListA.emplace_back(ir, 5, -jacobian(0, 2));
		tripletListA.emplace_back(ir, 6, -jacobian(0, 3));

		tripletListP.emplace_back(ir, ir, 1000000.0);

		tripletListB.emplace_back(ir, 0, delta);

		Eigen::SparseMatrix<double> matA(tripletListB.size(), 7);
		Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
		Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

		matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
		matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
		matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

		Eigen::SparseMatrix<double> AtPA(7, 7);
		Eigen::SparseMatrix<double> AtPB(7, 1);

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
			}
		}

		if (h_x.size() == 7)
		{
			for (size_t i = 0; i < h_x.size(); i++)
			{
				std::cout << h_x[i] << std::endl;
			}
			std::cout << "AtPA=AtPB SOLVED" << std::endl;
			std::cout << "update" << std::endl;

			int counter = 0;

			QuaternionPose pose = pose_quaternion_from_affine_matrix(pose_source);
			pose.px += h_x[counter++];
			pose.py += h_x[counter++];
			pose.pz += h_x[counter++];
			pose.q0 += h_x[counter++];
			pose.q1 += h_x[counter++];
			pose.q2 += h_x[counter++];
			pose.q3 += h_x[counter++];

			pose_source = affine_matrix_from_pose_quaternion(pose);
		}
		else
		{
			std::cout << "AtPA=AtPB FAILED" << std::endl;
		}
		break;
	}
	case 'u':
	{
		Eigen::MatrixXd d2sum_dbeta2(6, 6);
		d2sum_dbeta2 = Eigen::MatrixXd::Zero(6, 6);

		Eigen::MatrixXd d2sum_dbetadx(6, 6 * points_source_local.size());
		d2sum_dbetadx = Eigen::MatrixXd::Zero(6, 6 * points_source_local.size());

		TaitBryanPose pose_s = pose_tait_bryan_from_affine_matrix(pose_source);

		for (size_t i = 0; i < points_source_local.size(); i++)
		{
			Eigen::Vector3d &p_s = points_source_local[i];
			Eigen::Vector3d &p_t = points_target_global[i];

			Eigen::Matrix<double, 6, 6, Eigen::RowMajor> d2sum_dbeta2i;
			point_to_point_source_to_target_tait_bryan_wc_d2sum_dbeta2(
				d2sum_dbeta2i,
				pose_s.px, pose_s.py, pose_s.pz, pose_s.om, pose_s.fi, pose_s.ka,
				p_s.x(), p_s.y(), p_s.z(), p_t.x(), p_t.y(), p_t.z());

			for (size_t r = 0; r < 6; r++)
			{
				for (size_t c = 0; c < 6; c++)
				{
					d2sum_dbeta2(r, c) =
						d2sum_dbeta2(r, c) + d2sum_dbeta2i(r, c);
				}
			}
		}

		for (size_t i = 0; i < points_source_local.size(); i++)
		{
			Eigen::Vector3d &p_s = points_source_local[i];
			Eigen::Vector3d &p_t = points_target_global[i];

			Eigen::Matrix<double, 6, 6, Eigen::RowMajor> d2sum_dbetadxi(6, 6);

			point_to_point_source_to_target_tait_bryan_wc_d2sum_dbetadx(
				d2sum_dbetadxi,
				pose_s.px, pose_s.py, pose_s.pz, pose_s.om, pose_s.fi, pose_s.ka,
				p_s.x(), p_s.y(), p_s.z(), p_t.x(), p_t.y(), p_t.z());

			int cal = i * 6;
			for (size_t r = 0; r < 6; r++)
			{
				for (size_t c = 0; c < 6; c++)
				{
					d2sum_dbetadx(r, cal + c) = d2sum_dbetadxi(r, c);
				}
			}
		}
		Eigen::MatrixXd cov_x(6 * points_source_local.size(), 6 * points_source_local.size());
		cov_x = Eigen::MatrixXd::Zero(6 * points_source_local.size(), 6 * points_source_local.size());

		for (int i = 0; i < 6 * points_source_local.size(); i += 6)
		{
			// source local
			cov_x(i, i) = 0.01 * 0.01;
			cov_x(i + 1, i + 1) = 0.01 * 0.01;
			cov_x(i + 2, i + 2) = 0.01 * 0.01;
			// target global
			cov_x(i + 3, i + 3) = 0.01 * 0.01;
			cov_x(i + 4, i + 4) = 0.01 * 0.01;
			cov_x(i + 5, i + 5) = 0.01 * 0.01;
		}

		auto cov_b = d2sum_dbeta2.inverse() * d2sum_dbetadx * cov_x * d2sum_dbetadx.transpose() * d2sum_dbeta2.inverse();
		std::cout << "cov_b" << std::endl;
		std::cout << cov_b << std::endl;

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
	std::cout << "n: add noise to pose source" << std::endl;
	std::cout << "t: optimize (Tait-Bryan)" << std::endl;
	std::cout << "r: optimize (Rodrigues)" << std::endl;
	std::cout << "q: optimize (Quaternion)" << std::endl;
	std::cout << "u: uncertainty (Tait-Bryan)" << std::endl;
}
