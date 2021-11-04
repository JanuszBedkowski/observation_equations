#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

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

std::vector<Eigen::Affine3d> vertices;
std::vector<Eigen::Affine3d> vertices_odo;
int number_rows = 0;
int number_cols = 0;
std::vector<std::pair<int, int>> odo_edges;

#define RENDER_TYPE_CIRCLE 0
#define RENDER_TYPE_ROTATEDELIPSOID 1
#define RENDER_TYPE_ROTATEDBOUNDINGBOX 2
int render_type = RENDER_TYPE_CIRCLE;

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
	for(int x = -20; x <= 20; x += 1){
		for(int y = -20; y <= 20; y += 1){
			if(number_rows == 0){
				number_cols ++;
			}
			Eigen::Affine3d pose = Eigen::Affine3d::Identity();
			pose(0,3) = x;
			pose(1,3) = y;
			pose(2,3) = 0;
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
		}
	}


	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow("view_nns");
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

	switch(render_type){
		case RENDER_TYPE_CIRCLE:{
			glColor3f(1,0,0);
			glPointSize(4);
			glBegin(GL_POINTS);
			for(size_t i = 0 ; i < 10000; i++){
				float x = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 20.0;
				float y = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 20.0;

				if( (x-4)*(x-4) + (y-5)*(y-5) < 5*5 ){
					glColor3f(0,1,0);
				}else{
					glColor3f(1,0,0);
				}
				glVertex3f(x,y,0);
			}
			glEnd();
			glPointSize(1);

			glColor3f(0,0,1);
			glPointSize(20);
			glBegin(GL_POINTS);
			glVertex3f(4,5,0);
			glEnd();
			glPointSize(1);
			break;
		}
		case RENDER_TYPE_ROTATEDELIPSOID:{
			Eigen::Affine2d m = Eigen::Affine2d::Identity();
			double theta = 0.3;
			m(0,0) = cos(theta);
			m(0,1) = -sin(theta);
			m(1,0) = sin(theta);
			m(1,1) = cos(theta);
			m(0,2) = 4;
			m(1,2) = 5;
			glColor3f(1,0,0);
			glPointSize(4);
			glBegin(GL_POINTS);
			for(size_t i = 0 ; i < 10000; i++){
				float x = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 20.0;
				float y = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 20.0;

				Eigen::Vector2d v(x,y);
				Eigen::Vector2d vt = m.inverse() * v;


				if( (vt.x()*vt.x())/100 + (vt.y()*vt.y())/9 < 1 ){
					glColor3f(0,1,0);
				}else{
					glColor3f(1,0,0);
				}


				glVertex3f(x,y,0);
			}
			glEnd();
			glPointSize(1);

			glLineWidth(3);
			glColor3f(0,0,0);
			glBegin(GL_LINES);
				Eigen::Vector2d v(0,0);
				Eigen::Vector2d vt = m * v;
				glVertex3f(vt.x(),vt.y(),0);

				v = Eigen::Vector2d (10,0);
				vt = m * v;
				glVertex3f(vt.x(),vt.y(),0);
			glEnd();
			glLineWidth(1);


			glColor3f(0,0,1);
			glPointSize(20);
			glBegin(GL_POINTS);
			glVertex3f(4,5,0);
			glEnd();
			glPointSize(1);
			break;
		}

		case RENDER_TYPE_ROTATEDBOUNDINGBOX:{
			Eigen::Affine2d m = Eigen::Affine2d::Identity();
			double theta = 0.3;
			m(0,0) = cos(theta);
			m(0,1) = -sin(theta);
			m(1,0) = sin(theta);
			m(1,1) = cos(theta);
			m(0,2) = 4;
			m(1,2) = 5;
			glColor3f(1,0,0);
			glPointSize(4);
			glBegin(GL_POINTS);
			for(size_t i = 0 ; i < 10000; i++){
				float x = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 20.0;
				float y = ((float(rand()%1000000))/1000000.0f - 0.5) * 2.0 * 20.0;

				Eigen::Vector2d v(x,y);
				Eigen::Vector2d vt = m.inverse() * v;

				if(fabs(vt.x()) < 10 && fabs(vt.y()) < 3){
					glColor3f(0,1,0);
				}else{
					glColor3f(1,0,0);
				}
				glVertex3f(x,y,0);
			}
			glEnd();
			glPointSize(1);

			glLineWidth(3);
			glColor3f(0,0,0);
			glBegin(GL_LINES);
				Eigen::Vector2d v(0,0);
				Eigen::Vector2d vt = m * v;
				glVertex3f(vt.x(),vt.y(),0);

				v = Eigen::Vector2d (10,0);
				vt = m * v;
				glVertex3f(vt.x(),vt.y(),0);
			glEnd();
			glLineWidth(1);


			glColor3f(0,0,1);
			glPointSize(20);
			glBegin(GL_POINTS);
			glVertex3f(4,5,0);
			glEnd();
			glPointSize(1);
			break;
		}

	}



	glutSwapBuffers();
}



void keyboard(unsigned char key, int /*x*/, int /*y*/) {
	switch (key) {
		case (27): {
			glutDestroyWindow(glutGetWindow());
			return;
		}
		case '1':{
			render_type = RENDER_TYPE_CIRCLE;
			break;
		}
		case '2':{
			render_type = RENDER_TYPE_ROTATEDELIPSOID;
			break;
		}
		case '3':{
			render_type = RENDER_TYPE_ROTATEDBOUNDINGBOX;
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
	std::cout << "1: RENDER_TYPE_CIRCLE" << std::endl;
	std::cout << "2: RENDER_TYPE_ROTATEDELIPSOID" << std::endl;
	std::cout << "3: RENDER_TYPE_ROTATEDBOUNDINGBOX" << std::endl;
}







