#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>
#include <fstream>

struct ScanPose{
	Eigen::Affine3d m;
	std::vector<Eigen::Vector3d> pc;
	std::string scan_file_name;
};

std::vector<ScanPose> scan_poses;

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

void split(std::string &str, char delim, std::vector<std::string> &out)
{
	size_t start;
	size_t end = 0;

	while ((start = str.find_first_not_of(delim, end)) != std::string::npos)
	{
		end = str.find(delim, start);
		out.push_back(str.substr(start, end - start));
	}
}

int main(int argc, char *argv[]){
	if (false == initGL(&argc, argv)) {
		return 4;
	}

	std::ifstream f;
	f.open("../data/poses.csv");
	if(f.good()) {
		std::string s;
		while(!f.eof())	{
			getline(f,s);
			std::vector<std::string> strs;
			split(s, ',' , strs);
			ScanPose sp;
			if(strs.size() >= 13){
				sp.m = Eigen::Affine3d::Identity();
				std::istringstream(strs[0]) >> sp.scan_file_name;
				std::istringstream(strs[1]) >> sp.m(0,0);
				std::istringstream(strs[2]) >> sp.m(1,0);
				std::istringstream(strs[3]) >> sp.m(2,0);
				std::istringstream(strs[4]) >> sp.m(0,1);
				std::istringstream(strs[5]) >> sp.m(1,1);
				std::istringstream(strs[6]) >> sp.m(2,1);
				std::istringstream(strs[7]) >> sp.m(0,2);
				std::istringstream(strs[8]) >> sp.m(1,2);
				std::istringstream(strs[9]) >> sp.m(2,2);
				std::istringstream(strs[10]) >> sp.m(0,3);
				std::istringstream(strs[11]) >> sp.m(1,3);
				std::istringstream(strs[12]) >> sp.m(2,3);
				std::cout << sp.m(0,3) << " " << sp.m(1,3) << " " << sp.m(2,3) << std::endl;
				scan_poses.push_back(sp);
			}

		}
	}else{
		std::cout << "problem with opening file: ../data/poses.csv" << std::endl;
	}
	f.close();

	for(size_t i = 0; i < scan_poses.size(); i++){
		std::cout << scan_poses[i].scan_file_name << std::endl;
		std::ifstream f;
		f.open("../data/" + scan_poses[i].scan_file_name);
		if(f.good()) {
			std::string s;
			while(!f.eof())	{
				getline(f,s);
				std::vector<std::string> strs;
				split(s, ',' , strs);
				Eigen::Vector3d point;
				if(strs.size() == 3){
					std::istringstream(strs[0]) >> point.x();
					std::istringstream(strs[1]) >> point.y();
					std::istringstream(strs[2]) >> point.z();
					scan_poses[i].pc.push_back(point);
				}
			}
		}else{
			std::cout << "problem with opening file: ../data/poses.csv" << std::endl;
		}
		f.close();
	}

	std::cout << "loading done" << std::endl;

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
	glutCreateWindow("particle_filer_demo");
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

	glColor3f(1,0,0);
	glBegin(GL_POINTS);
	for(size_t i = 0 ; i < scan_poses.size(); i++){
		for(size_t j = 0 ; j < scan_poses[i].pc.size(); j++){
			Eigen::Vector3d p = scan_poses[i].m * scan_poses[i].pc[j];
			glVertex3f(p.x(), p.y(), p.z());
		}
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

}







