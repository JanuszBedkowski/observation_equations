#include <iostream>
#include <GL/freeglut.h>
#include <deque>
#include <vector>
#include <Eigen/Eigen>

#include "structures.h"
#include "transformations.h"
#include "point_to_point_source_to_target_mirror_tait_bryan_wc_jacobian.h"

struct Mirror{
	std::vector<Eigen::Vector3d> vertexes_local;
	Eigen::Affine3d pose;
	double a;
	double b;
	double c;
	double d;
};

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

Mirror mirror;
TaitBryanPose origindir;
Eigen::Vector3d origin(0,0,0);
Eigen::Vector3d dir(0,0,1);
double dir_length = 10;

void draw_mirror(const Mirror &mirror);
bool check_if_ray_intersect_mirror(const Eigen::Vector3d& origin,const Eigen::Vector3d& dir, const Mirror &mirror);
void update_mirror_abcd_from_pose(Mirror &mirror);

int main(int argc, char *argv[]){
	mirror.pose = Eigen::Affine3d::Identity();
	mirror.vertexes_local.emplace_back(0,1,0);
	mirror.vertexes_local.emplace_back(0,0,0);
	mirror.vertexes_local.emplace_back(1,0,0);
	mirror.pose(2,3) = 0.5;

	update_mirror_abcd_from_pose(mirror);

	origindir.om = origindir.fi = origindir.ka = 0;
	origindir.px = origindir.py = origindir.pz = 0;

	Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(origindir);
	origin.x() = m(0,3);
	origin.y() = m(1,3);
	origin.z() = m(2,3);

	dir.x() = m(0,2);
	dir.y() = m(1,2);
	dir.z() = m(2,2);

	dir*=dir_length;

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
	glutCreateWindow("mirror-sim");
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

	draw_mirror(mirror);

	glColor3f(1,0,0);
	glBegin(GL_LINES);
		glVertex3f(origin.x(), origin.y(), origin.z());
		glVertex3f(origin.x() + dir.x(), origin.y() + dir.y(), origin.z() + dir.z());
	glEnd();

	bool is_intersection =
			check_if_ray_intersect_mirror(origin, dir, mirror);

	Eigen::Matrix<double, 3, 1> intersection;
	double ray_length = (dir-origin).norm();
	Eigen::Vector3d dirn = dir/dir.norm();
	point_to_point_source_to_target_mirror_tait_bryan_wc_get_intersection(
			intersection, 0, 0, 0, 0, 0, 0,
			dirn.x(), dirn.y(), dirn.z(), ray_length,
			mirror.a,
			mirror.b,
			mirror.c,
			mirror.d);

	if(is_intersection){
		glPointSize(5);
		glColor3f(1,0,0);
		glBegin(GL_POINTS);
			glVertex3f(intersection(0,0), intersection(1,0), intersection(2,0));
		glEnd();
		glPointSize(1);

		double x = 0;
		double y = 0;
		double z = 0;
		transform_point_mirror_tait_bryan_wc(
				x, y, z,
				0, 0, 0, 0, 0, 0,
				dirn.x(), dirn.y(), dirn.z(), ray_length,
				mirror.a,
				mirror.b,
				mirror.c,
				mirror.d);
		glLineWidth(2);
		glColor3f(0,0,0);
		glBegin(GL_LINES);
			glVertex3f(intersection(0,0), intersection(1,0), intersection(2,0));
			glVertex3f(x, y, z);
		glEnd();
		glLineWidth(1);
	}
	glutSwapBuffers();
}


void keyboard(unsigned char key, int /*x*/, int /*y*/) {
	switch (key) {
		case (27): {
			glutDestroyWindow(glutGetWindow());
			return;
		}
		case 'a':{
			origindir.om += 0.01;

			Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(origindir);
			origin.x() = m(0,3);
			origin.y() = m(1,3);
			origin.z() = m(2,3);

			dir.x() = m(0,2);
			dir.y() = m(1,2);
			dir.z() = m(2,2);
			dir*=dir_length;
			break;
		}
		case 'd':{
			origindir.om -= 0.01;

			Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(origindir);
			origin.x() = m(0,3);
			origin.y() = m(1,3);
			origin.z() = m(2,3);

			dir.x() = m(0,2);
			dir.y() = m(1,2);
			dir.z() = m(2,2);
			dir*=dir_length;
			break;
		}
		case 'w':{
			origindir.fi += 0.01;

			Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(origindir);
			origin.x() = m(0,3);
			origin.y() = m(1,3);
			origin.z() = m(2,3);

			dir.x() = m(0,2);
			dir.y() = m(1,2);
			dir.z() = m(2,2);
			dir*=dir_length;
			break;
		}
		case 's':{
			origindir.fi -= 0.01;

			Eigen::Affine3d m = affine_matrix_from_pose_tait_bryan(origindir);
			origin.x() = m(0,3);
			origin.y() = m(1,3);
			origin.z() = m(2,3);

			dir.x() = m(0,2);
			dir.y() = m(1,2);
			dir.z() = m(2,2);
			dir*=dir_length;
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
	std::cout << "adws: change ray direction" << std::endl;
}

void draw_mirror(const Mirror &mirror)
{
	glColor3f(0.5,0.3,0.8);
	glBegin(GL_LINE_STRIP);
		for(const auto &v:mirror.vertexes_local){
			Eigen::Vector3d vertex = mirror.pose * v;
			glVertex3f(vertex.x(), vertex.y(), vertex.z());
//			std::cout << vertex << " ----" << std::endl;
		}
		Eigen::Vector3d vertex = mirror.pose * mirror.vertexes_local[0];
		glVertex3f(vertex.x(), vertex.y(), vertex.z());
	glEnd();
}

bool check_if_ray_intersect_mirror(const Eigen::Vector3d& origin,const Eigen::Vector3d& dir, const Mirror &mirror)
{
    if (dir == origin) return false;

    // project on plane
    Eigen::Vector3d dirn = dir/dir.norm();

    Eigen::Vector3d intersection;
    double a = mirror.a * dirn.x() + mirror.b * dirn.y() + mirror.c * dirn.z();
    double dist = mirror.a * origin.x() + mirror.b * origin.y() + mirror.c * origin.z() + mirror.d;

    intersection = origin - dirn * (dist/a);

    // project intersection on 2D
    Eigen::Vector3d intersection_2D = mirror.pose.inverse() * intersection;
    Eigen::Vector2d point_2D = intersection_2D.head<2>();

    double angle = 0;
    for (int i = 0; i < mirror.vertexes_local.size()-1; i++)
    {
    	Eigen::Vector2d vertex1 = mirror.vertexes_local[i].head<2>();
    	Eigen::Vector2d vertex2 = mirror.vertexes_local[i+1].head<2>();

        Eigen::Vector2d v1 = vertex1 - point_2D;
        Eigen::Vector2d v2 = vertex2 - point_2D;
        angle += acos(v1.dot(v2)/(v1.norm()*v2.norm()));
    }
    Eigen::Vector2d vertex1 = mirror.vertexes_local.back().head<2>();
    Eigen::Vector2d vertex2 = mirror.vertexes_local.front().head<2>();
    Eigen::Vector2d v1 = vertex1 - point_2D;
    Eigen::Vector2d v2 = vertex2 - point_2D;
    angle += acos(v1.dot(v2)/(v1.norm()*v2.norm()));

    if(std::abs(angle - 2*M_PI) < 1e-5){
        return true;
    }
    return false;
}

void update_mirror_abcd_from_pose(Mirror &mirror){
	mirror.a = mirror.pose(0,2);
	mirror.b = mirror.pose(1,2);
	mirror.c = mirror.pose(2,2);
	mirror.d = -mirror.a * mirror.pose(0,3) - mirror.b * mirror.pose(1,3) - mirror.c * mirror.pose(2,3);
}
