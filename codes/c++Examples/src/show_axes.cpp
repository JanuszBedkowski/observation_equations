#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

#include "structures.h"
#include "transformations.h"
#include "distance_point_to_plane_jacobian.h"


const unsigned int window_width = 800;
const unsigned int window_height = 600;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -20.0;
float translate_x, translate_y = 0.0;

Eigen::Affine3d m1 = Eigen::Affine3d::Identity();
Eigen::Affine3d m2 = Eigen::Affine3d::Identity();


bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

int main(int argc, char *argv[])
{
	m1(0, 0) = 0.00593569;
	m1(0,1) = 0.99995300;
	m1(0,2) = 0.00764619;
	m1(0,3) = 0.09058620;
	m1(1, 0) = - 0.99996700;
	m1(1, 1) = 0.00597756;
	m1(1, 2) = - 0.00546471;
	m1(1,3) = 0.89661100;
	m1(2,0) = - 0.00551016;
	m1(2,1) = - 0.00761350;
	m1(2,2) = 0.99995600;
	m1(2,3) = 1.03661000;
	m1(3,0) = 0;
	m1(3,1) = 0;
	m1(3,2) = 0;
	m1(3,3) = 1;

	m2(0,0) = 0.00709691;
	m2(0,1) = -0.99991192;
	m2(0,2) = 0.01121514;
	m2(0,3) = 1.10787461;
	m2(1,0) = 0.99996456;
	m2(1,1) = 0.00704560;
	m2(1,2) =-0.00460851;
	m2(1,3) = 0.89691993;
	m2(2,0) = 0.00452908;
	m2(2,1) = 0.01124745;
	m2(2,2) = 0.99992649;
	m2(2,3) =0.48735229;
	m2(3,0) = 0;
	m2(3,1) = 0;
	m2(3,2) = 0;
	m2(3,3) = 1;

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
	glutCreateWindow("show_axes");
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
	glVertex3f(10.0f, 0.0f, 0.0f);

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 10.0f, 0.0f);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 10.0f);
	glEnd();

	
	//m_imu
	glLineWidth(1);
	glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(m1(0, 3), m1(1, 3), m1(2, 3));
		glVertex3f(m1(0, 3) + m1(0, 0), m1(1, 3) + m1(1, 0), m1(2, 3) + m1(2, 0));

		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(m1(0, 3), m1(1, 3), m1(2, 3));
		glVertex3f(m1(0, 3) + m1(0, 1), m1(1, 3) + m1(1, 1), m1(2, 3) + m1(2, 1));

		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(m1(0, 3), m1(1, 3), m1(2, 3));
		glVertex3f(m1(0, 3) + m1(0, 2), m1(1, 3) + m1(1, 2), m1(2, 3) + m1(2, 2));
	glEnd();
	glLineWidth(1);

	glLineWidth(4);
	glBegin(GL_LINES);
		glColor3f(1.0f, 0.5f, 0.5f);
		glVertex3f(m2(0, 3), m2(1, 3), m2(2, 3));
		glVertex3f(m2(0, 3) + m2(0, 0), m2(1, 3) + m2(1, 0), m2(2, 3) + m2(2, 0));

		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(m2(0, 3), m2(1, 3), m2(2, 3));
		glVertex3f(m2(0, 3) + m2(0, 1), m2(1, 3) + m2(1, 1), m2(2, 3) + m2(2, 1));

		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(m2(0, 3), m2(1, 3), m2(2, 3));
		glVertex3f(m2(0, 3) + m2(0, 2), m2(1, 3) + m2(1, 2), m2(2, 3) + m2(2, 2));
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
		

		break;
	}
	case 'n':
	{
		
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
