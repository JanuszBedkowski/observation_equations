#include <Eigen/Eigen>
#include <GL/freeglut.h>
#include <iostream>
#include <cmath>   

#include <slerp_point_to_point_source_to_target_quaternion_wc_jacobian.h>
#include <structures.h>
#include <transformations.h>

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

class PoseInterpolation {
public:
	PoseInterpolation() {
		time_stamp_1 = 1000000;
		time_stamp_2 = 2000000;
		time_stamp_middle = 1200000;

		Eigen::AngleAxisd xAngle1(0.1, Eigen::Vector3d::UnitX());
		Eigen::AngleAxisd yAngle1(0.2, Eigen::Vector3d::UnitY());
		Eigen::AngleAxisd zAngle1(0.3, Eigen::Vector3d::UnitZ());
		Eigen::Quaternion<double> q1 = xAngle1 * yAngle1 * zAngle1;
		Eigen::Matrix3d rotationMatrix1 = q1.matrix();

		m1.linear() = rotationMatrix1;
		m1.translation() = Eigen::Vector3d(1, 1, 1);

		Eigen::AngleAxisd xAngle2(-1.6, Eigen::Vector3d::UnitX());
		Eigen::AngleAxisd yAngle2(0.1, Eigen::Vector3d::UnitY());
		Eigen::AngleAxisd zAngle2(-0.3, Eigen::Vector3d::UnitZ());
		Eigen::Quaternion<double> q2 = xAngle2 * yAngle2 * zAngle2;
		Eigen::Matrix3d rotationMatrix2 = q2.matrix();

		m1.linear() = rotationMatrix1;
		m1.translation() = Eigen::Vector3d(1, 1, 1);

		m2.linear() = rotationMatrix2;
		m2.translation() = Eigen::Vector3d(10, -1, 5);
	
	};
	~PoseInterpolation() { ; };

	Eigen::Affine3d pose_interpolation(double t, double t1, double t2,
		Eigen::Affine3d const& aff1, Eigen::Affine3d const& aff2);

	Eigen::Affine3d pose_interpolation_Eigen(double t, double t1, double t2,
		Eigen::Affine3d const& aff1, Eigen::Affine3d const& aff2);

	float time_stamp_1;
	float time_stamp_2;
	float time_stamp_middle;
	Eigen::Affine3d m1 = Eigen::Affine3d::Identity();
	Eigen::Affine3d m2 = Eigen::Affine3d::Identity();
};

Eigen::Affine3d PoseInterpolation::pose_interpolation(double t, double t1, double t2,
    Eigen::Affine3d const& aff1, Eigen::Affine3d const& aff2) {
    // assume here t1 <= t <= t2
    //double alpha = 0.0;
    //if (t2 != t1)
    //   alpha = (t - t1) / (t2 - t1);

    //Eigen::Quaternion<double> rot1(aff1.linear());
    //Eigen::Quaternion<double> rot2(aff2.linear());

    //Eigen::Vector3d trans1 = aff1.translation();
    //Eigen::Vector3d trans2 = aff2.translation();

    Eigen::Affine3d result = Eigen::Affine3d::Identity();

    QuaternionPose pq_0 = pose_quaternion_from_affine_matrix(aff1);
    QuaternionPose pq_1 = pose_quaternion_from_affine_matrix(aff2);

    double px;
    double py;
    double pz;
    double q0;
    double q1;
    double q2;
    double q3;

    double px_0 = pq_0.px;
    double py_0 = pq_0.py;
    double pz_0 = pq_0.pz;
    double q0_0 = pq_0.q0;
    double q1_0 = pq_0.q1;
    double q2_0 = pq_0.q2;
    double q3_0 = pq_0.q3;

    double px_1 = pq_1.px;
	double py_1 = pq_1.py;
	double pz_1 = pq_1.pz;
	double q0_1 = pq_1.q0;
	double q1_1 = pq_1.q1;
	double q2_1 = pq_1.q2;
	double q3_1 = pq_1.q3;

    slerp_point_to_point_source_to_target_quaternion_wc_interpolate(
    		px, py, pz, q0, q1, q2, q3,
    		px_0, py_0, pz_0, q0_0, q1_0, q2_0, q3_0,
    		px_1, py_1, pz_1, q0_1, q1_1, q2_1, q3_1,
    		t1, t2, t);

    QuaternionPose pq;
    pq.px = px;
    pq.py = py;
    pq.pz = pz;
    pq.q0 = q0;
    pq.q1 = q1;
    pq.q2 = q2;
    pq.q3 = q3;

    result = affine_matrix_from_pose_quaternion(pq);

    return result;
}

Eigen::Affine3d PoseInterpolation::pose_interpolation_Eigen(double query_time, double t1, double t2,
    Eigen::Affine3d const& aff1, Eigen::Affine3d const& aff2) {
	
	Eigen::Affine3d m_interpolated;
	m_interpolated = Eigen::Matrix4d::Identity();
    double res = (query_time - t1) / (t2 - t1);

	const Eigen::Vector3d diff_translation = aff2.translation() - aff1.translation();
        
	m_interpolated.translation() = aff2.translation() + diff_translation * res;

	// Eigen::Matrix3d r1 = aff1.linear();
	// Eigen::Matrix3d r2 = aff2.linear();
	Eigen::Matrix3d r1 = aff1.rotation();
	Eigen::Matrix3d r2 = aff2.rotation();
	Eigen::Quaterniond q1(r1);
	Eigen::Quaterniond q2(r2);
	Eigen::Quaterniond qt = q1.slerp(res, q2);
	
	m_interpolated.matrix().topLeftCorner(3, 3) = qt.toRotationMatrix();
	return m_interpolated;
}

PoseInterpolation pi;

int main(int argc, char *argv[]){
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
	glutCreateWindow("pose_interpolation");
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
	gluPerspective(60.0, (GLfloat) window_width / (GLfloat) window_height, 0.01, 10000.0);
	glutReshapeFunc(reshape);

	return true;
}

void draw_axes(Eigen::Affine3d m){
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

void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(translate_x, translate_y, translate_z);
	glRotatef(rotate_x, 1.0, 0.0, 0.0);
	glRotatef(rotate_y, 0.0, 0.0, 1.0);

	draw_axes(pi.m1);
	draw_axes(pi.m2);

	//Eigen::Affine3d m3 = pi.pose_interpolation(pi.time_stamp_middle, pi.time_stamp_1, pi.time_stamp_2, pi.m1, pi.m2);
	Eigen::Affine3d m3 = pi.pose_interpolation_Eigen(pi.time_stamp_middle, pi.time_stamp_1, pi.time_stamp_2, pi.m1, pi.m2);
	
	draw_axes(m3);

	glutSwapBuffers();
}


void keyboard(unsigned char key, int /*x*/, int /*y*/) {
	switch (key) {
		case (27): {
			glutDestroyWindow(glutGetWindow());
			return;
		}
		case '-':{
			pi.time_stamp_middle -= 10000;
			break;
		}
		case '=':{
			pi.time_stamp_middle += 10000;
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
	std::cout << "=: move pose forward" << std::endl;
	std::cout << "-: move pose backward" << std::endl;
}

