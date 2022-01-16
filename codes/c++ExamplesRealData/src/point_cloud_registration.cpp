#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <thread>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>

#include "structures.h"
#include "transformations.h"
#include "cauchy.h"

#include "point_to_point_tait_bryan_wc_jacobian.h"
#include "point_to_point_source_to_target_tait_bryan_wc_jacobian.h"

#include "rgd.h"

struct ScanPose{
	Eigen::Affine3d m;
	pcl::PointCloud<pcl::PointXYZ> pc;
};

pcl::PointCloud<pcl::PointXYZ> pc_ground_truth;
std::vector<ScanPose> scan_poses;
int current_scan_index = 0;

float sradius = 1.0;
bool show_ground_truth = true;

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
void set_initial_guess(std::vector<ScanPose>& scan_poses);

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

std::vector<std::pair<int,int>> nns(ScanPose &sp1, ScanPose &sp2, float radius);

std::vector<std::pair<int,int>> pairs_temp;
std::vector<Bucket> buckets_render;
bool show_ndt_covariances = true;

int main(int argc, char *argv[]){

	std::vector<std::string> pcd_file_names;
	pcd_file_names.push_back("../data/pcd/scan000.pcd");
	pcd_file_names.push_back("../data/pcd/scan001.pcd");
	pcd_file_names.push_back("../data/pcd/scan002.pcd");
	pcd_file_names.push_back("../data/pcd/scan003.pcd");
	pcd_file_names.push_back("../data/pcd/scan004.pcd");
	pcd_file_names.push_back("../data/pcd/scan005.pcd");
	pcd_file_names.push_back("../data/pcd/scan006.pcd");
	pcd_file_names.push_back("../data/pcd/scan007.pcd");
	pcd_file_names.push_back("../data/pcd/scan008.pcd");
	pcd_file_names.push_back("../data/pcd/scan009.pcd");
	pcd_file_names.push_back("../data/pcd/scan010.pcd");
	pcd_file_names.push_back("../data/pcd/scan011.pcd");
	pcd_file_names.push_back("../data/pcd/scan012.pcd");

	for(size_t i = 0; i < pcd_file_names.size(); i++){
		std::cout << "loading file: " << pcd_file_names[i] << std::endl;
		ScanPose sp;
		sp.m = Eigen::Affine3d::Identity();
		pcl::PointCloud<pcl::PointXYZ> pc;
		if (pcl::io::loadPCDFile(pcd_file_names[i], pc) == -1) {
			std::cout << "PROBLEM WITH LODAING pcd: " << pcd_file_names[i] << std::endl;
			return 1;
		}else{
			sp.pc = pc;
			scan_poses.push_back(sp);
		}
	}

	if (pcl::io::loadPCDFile("../data/pcd/ground_truth.pcd", pc_ground_truth) == -1) {
		std::cout << "PROBLEM WITH LODAING pcd: ../data/pcd/ground_truth.pcd" << std::endl;
		return 1;
	}

	set_initial_guess(scan_poses);

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
	glutCreateWindow("point_cloud_registration");
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

	glBegin(GL_POINTS);
	glColor3f(1.0, 0.0, 0.0);
	for(size_t i = 0; i < scan_poses.size(); i++){
		if(i == current_scan_index){
			glColor3f(0.0, 1.0, 0.0);
		}else{
			glColor3f(1.0, 0.0, 0.0);
		}
		if(i+1 == current_scan_index){
			glColor3f(0.0, 0.0, 1.0);
		}
		for(size_t j = 0; j < scan_poses[i].pc.size(); j++){
			Eigen::Vector3d v(scan_poses[i].pc[j].x, scan_poses[i].pc[j].y, scan_poses[i].pc[j].z);
			Eigen::Vector3d vt = scan_poses[i].m * v;
			glVertex3f(vt.x(), vt.y(), vt.z());
		}
	}
	glEnd();

	glColor3f(0,0,0);
	glBegin(GL_LINES);
	for(size_t i = 0 ; i < pairs_temp.size(); i++){
		pcl::PointXYZ p1 = scan_poses[0].pc[pairs_temp[i].first];
		pcl::PointXYZ p2 = scan_poses[1].pc[pairs_temp[i].second];

		Eigen::Vector3d v1(p1.x, p1.y, p1.z);
		Eigen::Vector3d v1t = scan_poses[0].m * v1;

		glVertex3f(v1t.x(), v1t.y(), v1t.z());

		Eigen::Vector3d v2(p2.x, p2.y, p2.z);
		Eigen::Vector3d v2t = scan_poses[1].m * v2;

		glVertex3f(v2t.x(), v2t.y(), v2t.z());

	}
	glEnd();

	if(show_ground_truth){
		glColor3f(0.7, 0.7, 0.7);
		glBegin(GL_POINTS);
		for(size_t i = 0; i < pc_ground_truth.size(); i++){
			glVertex3f(pc_ground_truth[i].x, pc_ground_truth[i].y, pc_ground_truth[i].z);
		}
		glEnd();
	}

	if(show_ndt_covariances){
		for(size_t i = 0 ; i < buckets_render.size(); i++){
			if(buckets_render[i].number_of_points > 10){
				draw_ellipse(buckets_render[i].cov, buckets_render[i].mean, Eigen::Vector3f(0.0, 0.0, 1.0), 1);
			}
		}
	}

	glutSwapBuffers();
}

void ndt_job(int i, Job* job, std::vector<Bucket>* buckets, Eigen::SparseMatrix<double > *AtPA,
		Eigen::SparseMatrix<double > *AtPB, std::vector<PointBucketIndexPair> *index_pair_internal, std::vector<Point3D> *pp,
		std::vector<TaitBryanPose> *poses, std::vector<Eigen::Affine3d> *mposes_inv, size_t trajectory_size) {

	std::vector<Eigen::Triplet<double>> tripletListA;
	std::vector<Eigen::Triplet<double>> tripletListP;
	std::vector<Eigen::Triplet<double>> tripletListB;

	for (size_t ii = job->index_begin_inclusive; ii < job->index_end_exclusive; ii++) {
		Bucket& b = (*buckets)[ii];
		if (b.number_of_points < 5)continue;

		Eigen::Vector3d mean(0, 0, 0);
		Eigen::Matrix3d cov;
		cov.setZero();

		for (int index = b.index_begin; index < b.index_end; index++) {
			const auto& p = (*pp)[(*index_pair_internal)[index].index_of_point];
			mean += Eigen::Vector3d(p.x, p.y, p.z);
		}
		mean /= b.number_of_points;

		for (int index = b.index_begin; index < b.index_end; index++) {
			const auto& p = (*pp)[(*index_pair_internal)[index].index_of_point];
			cov(0, 0) += (mean.x() - p.x) * (mean.x() - p.x);
			cov(0, 1) += (mean.x() - p.x) * (mean.y() - p.y);
			cov(0, 2) += (mean.x() - p.x) * (mean.z() - p.z);
			cov(1, 0) += (mean.y() - p.y) * (mean.x() - p.x);
			cov(1, 1) += (mean.y() - p.y) * (mean.y() - p.y);
			cov(1, 2) += (mean.y() - p.y) * (mean.z() - p.z);
			cov(2, 0) += (mean.z() - p.z) * (mean.x() - p.x);
			cov(2, 1) += (mean.z() - p.z) * (mean.y() - p.y);
			cov(2, 2) += (mean.z() - p.z) * (mean.z() - p.z);
		}
		cov /= b.number_of_points;

		(*buckets)[ii].mean = mean;
		(*buckets)[ii].cov = cov;


		Eigen::Matrix3d infm = cov.inverse();

		if (!(infm(0, 0) == infm(0, 0)))continue;
		if (!(infm(0, 1) == infm(0, 1)))continue;
		if (!(infm(0, 2) == infm(0, 2)))continue;

		if (!(infm(1, 0) == infm(1, 0)))continue;
		if (!(infm(1, 1) == infm(1, 1)))continue;
		if (!(infm(1, 2) == infm(1, 2)))continue;

		if (!(infm(2, 0) == infm(2, 0)))continue;
		if (!(infm(2, 1) == infm(2, 1)))continue;
		if (!(infm(2, 2) == infm(2, 2)))continue;



		for (int index = b.index_begin; index < b.index_end; index++) {
			const auto& p = (*pp)[(*index_pair_internal)[index].index_of_point];

			Eigen::Vector3d point_local(p.x, p.y, p.z);
			point_local = (*mposes_inv)[p.index_pose] * point_local;


			TaitBryanPose pose_s = (*poses)[p.index_pose];
			double delta_x;
			double delta_y;
			double delta_z;

			point_to_point_source_to_target_tait_bryan_wc(delta_x, delta_y, delta_z,
				pose_s.px, pose_s.py, pose_s.pz, pose_s.om, pose_s.fi, pose_s.ka,
				point_local.x(), point_local.y(), point_local.z(), mean.x(), mean.y(), mean.z());

			Eigen::Matrix<double, 3, 6, Eigen::RowMajor> jacobian;
			point_to_point_source_to_target_tait_bryan_wc_jacobian(jacobian,
				pose_s.px, pose_s.py, pose_s.pz, pose_s.om, pose_s.fi, pose_s.ka,
				point_local.x(), point_local.y(), point_local.z());


			int ir = tripletListB.size();
			int c = p.index_pose * 6;

			for (int row = 0; row < 3; row++) {
				for (int col = 0; col < 6; col++) {
					if (jacobian(row, col) != 0.0) {
						tripletListA.emplace_back(ir + row, c + col, -jacobian(row, col));
					}
				}
			}

			tripletListP.emplace_back(ir, ir, infm(0, 0));
			tripletListP.emplace_back(ir, ir + 1, infm(0, 1));
			tripletListP.emplace_back(ir, ir + 2, infm(0, 2));
			tripletListP.emplace_back(ir + 1, ir, infm(1, 0));
			tripletListP.emplace_back(ir + 1, ir + 1, infm(1, 1));
			tripletListP.emplace_back(ir + 1, ir + 2, infm(1, 2));
			tripletListP.emplace_back(ir + 2, ir, infm(2, 0));
			tripletListP.emplace_back(ir + 2, ir + 1, infm(2, 1));
			tripletListP.emplace_back(ir + 2, ir + 2, infm(2, 2));

			tripletListB.emplace_back(ir, 0, delta_x);
			tripletListB.emplace_back(ir + 1, 0, delta_y);
			tripletListB.emplace_back(ir + 2, 0, delta_z);
		}
	}

	Eigen::SparseMatrix<double> matA(tripletListB.size(), trajectory_size * 6);
	Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
	Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

	matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
	matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
	matB.setFromTriplets(tripletListB.begin(), tripletListB.end());


	Eigen::SparseMatrix<double> AtPAt(trajectory_size * 6, trajectory_size * 6);
	Eigen::SparseMatrix<double> AtPBt(trajectory_size * 6, 1);

	{
		Eigen::SparseMatrix<double> AtP = matA.transpose() * matP;
		AtPAt = AtP * matA;
		AtPBt = AtP * matB;

		(*AtPA) = AtPAt;
		(*AtPB) = AtPBt;
	}
}

void keyboard(unsigned char key, int /*x*/, int /*y*/) {
	switch (key) {
		case (27): {
			glutDestroyWindow(glutGetWindow());
			return;
		}
		case '-':{
			current_scan_index --;
			if(current_scan_index < 0)current_scan_index = 0;
			break;
		}
		case '=':{
			current_scan_index ++;
			if(current_scan_index >= scan_poses.size())current_scan_index = scan_poses.size() - 1;
			break;
		}
		case 'a':{
			TaitBryanPose pose;
			pose.px = 0;
			pose.py = -0.1;
			pose.pz = 0;
			pose.om = 0;
			pose.fi = 0;
			pose.ka = 0;
			scan_poses[current_scan_index].m = scan_poses[current_scan_index].m * affine_matrix_from_pose_tait_bryan(pose);
			break;
		}
		case 'd':{
			TaitBryanPose pose;
			pose.px = 0;
			pose.py = 0.1;
			pose.pz = 0;
			pose.om = 0;
			pose.fi = 0;
			pose.ka = 0;
			scan_poses[current_scan_index].m = scan_poses[current_scan_index].m * affine_matrix_from_pose_tait_bryan(pose);
			break;
		}
		case 'w':{
			TaitBryanPose pose;
			pose.px = 0.1;
			pose.py = 0;
			pose.pz = 0;
			pose.om = 0;
			pose.fi = 0;
			pose.ka = 0;
			scan_poses[current_scan_index].m = scan_poses[current_scan_index].m * affine_matrix_from_pose_tait_bryan(pose);
			break;
		}
		case 's':{
			TaitBryanPose pose;
			pose.px = -0.1;
			pose.py = 0;
			pose.pz = 0;
			pose.om = 0;
			pose.fi = 0;
			pose.ka = 0;
			scan_poses[current_scan_index].m = scan_poses[current_scan_index].m * affine_matrix_from_pose_tait_bryan(pose);
			break;
		}
		case 'z':{
			TaitBryanPose pose;
			pose.px = 0;
			pose.py = 0;
			pose.pz = 0;
			pose.om = 0;
			pose.fi = 0;
			pose.ka = -0.01;
			scan_poses[current_scan_index].m = scan_poses[current_scan_index].m * affine_matrix_from_pose_tait_bryan(pose);
			break;
		}
		case 'x':{
			TaitBryanPose pose;
			pose.px = 0;
			pose.py = 0;
			pose.pz = 0;
			pose.om = 0;
			pose.fi = 0;
			pose.ka = 0.01;
			scan_poses[current_scan_index].m = scan_poses[current_scan_index].m * affine_matrix_from_pose_tait_bryan(pose);
			break;
		}
		case 'p':{
			for(size_t i = 0; i < scan_poses.size(); i++){
				std::cout << "scan: " << i << std::endl;
				std::cout << scan_poses[i].m.matrix() << std::endl;
			}
			break;
		}
		case 'n':{
			pairs_temp = nns(scan_poses[0], scan_poses[1], sradius);
			break;
		}
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < scan_poses.size() ; i++){
				for(size_t j = i+1 ; j < scan_poses.size() ; j++){
					//if(i == j)continue;
					std::vector<std::pair<int,int>> nn = nns(scan_poses[i], scan_poses[j], sradius);
					std::cout << nn.size() << "," << scan_poses[i].pc.size() << "," << scan_poses[j].pc.size() << std::endl;

					TaitBryanPose pose_1 = pose_tait_bryan_from_affine_matrix(scan_poses[i].m);
					TaitBryanPose pose_2 = pose_tait_bryan_from_affine_matrix(scan_poses[j].m);

					for(size_t k = 0 ; k < nn.size(); k+=1){
						pcl::PointXYZ &p_1 = scan_poses[i].pc[nn[k].first];
						pcl::PointXYZ &p_2 = scan_poses[j].pc[nn[k].second];
						double delta_x;
						double delta_y;
						double delta_z;
						point_to_point_tait_bryan_wc(delta_x, delta_y, delta_z, pose_1.px, pose_1.py, pose_1.pz, pose_1.om, pose_1.fi, pose_1.ka, pose_2.px, pose_2.py, pose_2.pz, pose_2.om, pose_2.fi, pose_2.ka, p_1.x, p_1.y, p_1.z, p_2.x, p_2.y, p_2.z);

						Eigen::Matrix<double, 3, 12, Eigen::RowMajor> jacobian;
						point_to_point_tait_bryan_wc_jacobian(jacobian, pose_1.px, pose_1.py, pose_1.pz, pose_1.om, pose_1.fi, pose_1.ka, pose_2.px, pose_2.py, pose_2.pz, pose_2.om, pose_2.fi, pose_2.ka, p_1.x, p_1.y, p_1.z, p_2.x, p_2.y, p_2.z);

						int ir = tripletListB.size();
						int ic_1 = i * 6;
						int ic_2 = j * 6;

						tripletListA.emplace_back(ir     , ic_1    , -jacobian(0,0));
						tripletListA.emplace_back(ir     , ic_1 + 1, -jacobian(0,1));
						tripletListA.emplace_back(ir     , ic_1 + 2, -jacobian(0,2));
						tripletListA.emplace_back(ir     , ic_1 + 3, -jacobian(0,3));
						tripletListA.emplace_back(ir     , ic_1 + 4, -jacobian(0,4));
						tripletListA.emplace_back(ir     , ic_1 + 5, -jacobian(0,5));

						tripletListA.emplace_back(ir     , ic_2    , -jacobian(0,6));
						tripletListA.emplace_back(ir     , ic_2 + 1, -jacobian(0,7));
						tripletListA.emplace_back(ir     , ic_2 + 2, -jacobian(0,8));
						tripletListA.emplace_back(ir     , ic_2 + 3, -jacobian(0,9));
						tripletListA.emplace_back(ir     , ic_2 + 4, -jacobian(0,10));
						tripletListA.emplace_back(ir     , ic_2 + 5, -jacobian(0,11));

						tripletListA.emplace_back(ir + 1 , ic_1    , -jacobian(1,0));
						tripletListA.emplace_back(ir + 1 , ic_1 + 1, -jacobian(1,1));
						tripletListA.emplace_back(ir + 1 , ic_1 + 2, -jacobian(1,2));
						tripletListA.emplace_back(ir + 1 , ic_1 + 3, -jacobian(1,3));
						tripletListA.emplace_back(ir + 1 , ic_1 + 4, -jacobian(1,4));
						tripletListA.emplace_back(ir + 1 , ic_1 + 5, -jacobian(1,5));

						tripletListA.emplace_back(ir + 1 , ic_2    , -jacobian(1,6));
						tripletListA.emplace_back(ir + 1 , ic_2 + 1, -jacobian(1,7));
						tripletListA.emplace_back(ir + 1 , ic_2 + 2, -jacobian(1,8));
						tripletListA.emplace_back(ir + 1 , ic_2 + 3, -jacobian(1,9));
						tripletListA.emplace_back(ir + 1 , ic_2 + 4, -jacobian(1,10));
						tripletListA.emplace_back(ir + 1 , ic_2 + 5, -jacobian(1,11));

						tripletListA.emplace_back(ir + 2 , ic_1    , -jacobian(2,0));
						tripletListA.emplace_back(ir + 2 , ic_1 + 1, -jacobian(2,1));
						tripletListA.emplace_back(ir + 2 , ic_1 + 2, -jacobian(2,2));
						tripletListA.emplace_back(ir + 2 , ic_1 + 3, -jacobian(2,3));
						tripletListA.emplace_back(ir + 2 , ic_1 + 4, -jacobian(2,4));
						tripletListA.emplace_back(ir + 2 , ic_1 + 5, -jacobian(2,5));

						tripletListA.emplace_back(ir + 2 , ic_2    , -jacobian(2,6));
						tripletListA.emplace_back(ir + 2 , ic_2 + 1, -jacobian(2,7));
						tripletListA.emplace_back(ir + 2 , ic_2 + 2, -jacobian(2,8));
						tripletListA.emplace_back(ir + 2 , ic_2 + 3, -jacobian(2,9));
						tripletListA.emplace_back(ir + 2 , ic_2 + 4, -jacobian(2,10));
						tripletListA.emplace_back(ir + 2 , ic_2 + 5, -jacobian(2,11));

						tripletListP.emplace_back(ir    , ir    ,  1);//cauchy(delta_x, 1));
						tripletListP.emplace_back(ir + 1, ir + 1,  1);//cauchy(delta_y, 1));
						tripletListP.emplace_back(ir + 2, ir + 2,  1);//cauchy(delta_z, 1));

						tripletListB.emplace_back(ir    , 0,  delta_x);
						tripletListB.emplace_back(ir + 1, 0,  delta_y);
						tripletListB.emplace_back(ir + 2, 0,  delta_z);
					}
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

			Eigen::SparseMatrix<double> matA(tripletListB.size(), scan_poses.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(scan_poses.size() * 6, scan_poses.size() * 6);
			Eigen::SparseMatrix<double> AtPB(scan_poses.size() * 6, 1);

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

			if(h_x.size() == scan_poses.size() * 6){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < scan_poses.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(scan_poses[i].m);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];

					scan_poses[i].m = affine_matrix_from_pose_tait_bryan(pose);
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'l':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0 ; i < scan_poses.size() ; i++){
				Eigen::Affine3d pose_source = scan_poses[i].m;

				for(size_t j = 0 ; j < scan_poses.size() ; j++){
					if(i == j)continue;
					std::vector<std::pair<int,int>> nn = nns(scan_poses[i], scan_poses[j], sradius);
					std::cout << nn.size() << "," << scan_poses[i].pc.size() << "," << scan_poses[j].pc.size() << std::endl;

					for(size_t k = 0 ; k < nn.size(); k+=1){
						pcl::PointXYZ &p_1 = scan_poses[i].pc[nn[k].first];
						pcl::PointXYZ &p_2 = scan_poses[j].pc[nn[k].second];

						Eigen::Vector3d p_t(p_2.x, p_2.y, p_2.z);//
						Eigen::Vector3d p_s(p_1.x, p_1.y, p_1.z);//

						int ir = tripletListB.size();
						int ic_1 = i * 6;
						Eigen::Matrix3d px;
						px(0,0) = 0;
						px(0,1) = -p_s.z();
						px(0,2) =  p_s.y();
						px(1,0) = p_s.z();
						px(1,1) = 0;
						px(1,2) = -p_s.x();
						px(2,0) = -p_s.y();
						px(2,1) = p_s.x();
						px(2,2) = 0;

						Eigen::Matrix3d R = pose_source.rotation();
						Eigen::Matrix3d Rpx = R*px;

						tripletListA.emplace_back(ir     ,ic_1 + 0, R(0,0));
						tripletListA.emplace_back(ir     ,ic_1 + 1, R(0,1));
						tripletListA.emplace_back(ir     ,ic_1 + 2, R(0,2));
						tripletListA.emplace_back(ir     ,ic_1 + 3, -Rpx(0,0));
						tripletListA.emplace_back(ir     ,ic_1 + 4, -Rpx(0,1));
						tripletListA.emplace_back(ir     ,ic_1 + 5, -Rpx(0,2));

						tripletListA.emplace_back(ir + 1 ,ic_1 + 0, R(1,0));
						tripletListA.emplace_back(ir + 1 ,ic_1 + 1, R(1,1));
						tripletListA.emplace_back(ir + 1 ,ic_1 + 2, R(1,2));
						tripletListA.emplace_back(ir + 1 ,ic_1 + 3, -Rpx(1,0));
						tripletListA.emplace_back(ir + 1 ,ic_1 + 4, -Rpx(1,1));
						tripletListA.emplace_back(ir + 1 ,ic_1 + 5, -Rpx(1,2));

						tripletListA.emplace_back(ir + 2 ,ic_1 + 0, R(2,0));
						tripletListA.emplace_back(ir + 2 ,ic_1 + 1, R(2,1));
						tripletListA.emplace_back(ir + 2 ,ic_1 + 2, R(2,2));
						tripletListA.emplace_back(ir + 2 ,ic_1 + 3, -Rpx(2,0));
						tripletListA.emplace_back(ir + 2 ,ic_1 + 4, -Rpx(2,1));
						tripletListA.emplace_back(ir + 2 ,ic_1 + 5, -Rpx(2,2));

						tripletListP.emplace_back(ir    , ir    ,  1);
						tripletListP.emplace_back(ir + 1, ir + 1,  1);
						tripletListP.emplace_back(ir + 2, ir + 2,  1);

						Eigen::Vector3d target = scan_poses[j].m * p_t;
						Eigen::Vector3d source = scan_poses[i].m * p_s;

						tripletListB.emplace_back(ir    , 0,  target.x() - source.x());
						tripletListB.emplace_back(ir + 1, 0,  target.y() - source.y());
						tripletListB.emplace_back(ir + 2, 0,  target.z() - source.z());
					}
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

			Eigen::SparseMatrix<double> matA(tripletListB.size(), scan_poses.size() * 6);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(scan_poses.size() * 6, scan_poses.size() * 6);
			Eigen::SparseMatrix<double> AtPB(scan_poses.size() * 6, 1);

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

			if(h_x.size() == scan_poses.size() * 6){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < scan_poses.size(); i++){
					RodriguesPose pose_update;
					pose_update.px = h_x[counter++];
					pose_update.py = h_x[counter++];
					pose_update.pz = h_x[counter++];
					pose_update.sx = h_x[counter++];
					pose_update.sy = h_x[counter++];
					pose_update.sz = h_x[counter++];

					scan_poses[i].m = scan_poses[i].m * affine_matrix_from_pose_rodrigues(pose_update);
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case '1':{
			sradius -= 0.01;
			std::cout << "sradius: " << sradius << std::endl;
			if(sradius < 0)sradius = 0.01;
			break;		
		} 
		case '2':{
			sradius += 0.01;
			std::cout << "sradius: " << sradius << std::endl;
			break;		
		} 
		case '3':{
			for(size_t i = 0; i < scan_poses.size(); i++){
				pcl::PointCloud<pcl::PointXYZ> pc;
				for(size_t j = 0; j < scan_poses[i].pc.size(); j++){
					Eigen::Vector3d v(scan_poses[i].pc[j].x, scan_poses[i].pc[j].y, scan_poses[i].pc[j].z);
					Eigen::Vector3d vt = scan_poses[i].m * v;
					pcl::PointXYZ p;
					p.x = vt.x();
					p.y = vt.y();
					p.z = vt.z();
					pc.push_back(p);
				}
				pcl::io::savePCDFileBinary(std::to_string(i) + ".pcd", pc);
				std::cout << "file: " << std::to_string(i) + ".pcd" << std::endl;
 			}
			break;
		}
		case 'g':{
			show_ground_truth =! show_ground_truth;
			break;
		}
		case 'o':{
			GridParameters rgd_params;
			rgd_params.resolution_X = sradius;
			rgd_params.resolution_Y = sradius;
			rgd_params.resolution_Z = sradius;
			rgd_params.bounding_box_extension = sradius;

			std::vector<Point3D> points_global;
			for(size_t i = 0; i < scan_poses.size(); i++){
				for(size_t j = 0; j < scan_poses[i].pc.size(); j++){
					Eigen::Vector3d v(scan_poses[i].pc[j].x, scan_poses[i].pc[j].y, scan_poses[i].pc[j].z);
					Eigen::Vector3d vt = scan_poses[i].m * v;
					Point3D p;
					p.x = vt.x();
					p.y = vt.y();
					p.z = vt.z();
					p.index_pose = i;
					points_global.push_back(p);
				}
			}

			std::vector<PointBucketIndexPair> index_pair;
			std::vector<Bucket> buckets;

			grid_calculate_params(points_global, rgd_params);
			build_rgd(points_global, index_pair, buckets, rgd_params, 8);

			std::vector<Job> jobs = get_jobs(buckets.size(), 8);

			std::vector<std::thread> threads;

			std::vector<Eigen::SparseMatrix<double>> AtPAtmp(jobs.size());
			std::vector<Eigen::SparseMatrix<double>> AtPBtmp(jobs.size());

			for (size_t i = 0; i < jobs.size(); i++) {
				AtPAtmp[i] = Eigen::SparseMatrix<double>(scan_poses.size() * 6, scan_poses.size() * 6);
				AtPBtmp[i] = Eigen::SparseMatrix<double>(scan_poses.size() * 6, 1);
			}

			std::vector<TaitBryanPose> poses;
			std::vector<Eigen::Affine3d> mposes_inv;
			for(size_t i = 0; i < scan_poses.size(); i++){
				poses.push_back(pose_tait_bryan_from_affine_matrix(scan_poses[i].m));
				mposes_inv.push_back(scan_poses[i].m.inverse());
			}

			for (size_t k = 0; k < jobs.size(); k++) {
				threads.push_back(std::thread(ndt_job, k, &jobs[k], &buckets, &(AtPAtmp[k]), &(AtPBtmp[k]), &index_pair, &points_global, &poses, &mposes_inv, scan_poses.size()));
			}

			for (size_t j = 0; j < threads.size(); j++) {
				threads[j].join();
			}

			bool init = false;
			Eigen::SparseMatrix<double> AtPA_ndt(scan_poses.size() * 6, scan_poses.size() * 6);
			Eigen::SparseMatrix<double> AtPB_ndt(scan_poses.size() * 6, 1);

			for (size_t k = 0; k < jobs.size(); k++) {
				if (!init) {
					if (AtPBtmp[k].size() > 0) {
						AtPA_ndt = AtPAtmp[k];
						AtPB_ndt = AtPBtmp[k];
						init = true;
					}
				}
				else {
					if (AtPBtmp[k].size() > 0) {

						AtPA_ndt += AtPAtmp[k];
						AtPB_ndt += AtPBtmp[k];
					}
				}
			}
			std::cout << "start solving AtPA=AtPB" << std::endl;
			Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(AtPA_ndt);

			std::cout << "x = solver.solve(AtPB)" << std::endl;
			Eigen::SparseMatrix<double> x = solver.solve(AtPB_ndt);

			std::vector<double> h_x;

			for (int k = 0; k < x.outerSize(); ++k) {
				for (Eigen::SparseMatrix<double>::InnerIterator it(x, k); it; ++it) {
					if (it.value() == it.value()) {
						h_x.push_back(it.value());
						std::cout << it.value() << std::endl;
					}
				}
			}

			if(h_x.size() == scan_poses.size() * 6){
				buckets_render = buckets;
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				int counter = 0;
				for (size_t i = 0; i < scan_poses.size(); i++) {

					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(scan_poses[i].m);

					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];

					scan_poses[i].m = affine_matrix_from_pose_tait_bryan(pose);
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}

			break;
		}
		case '4':{
			show_ndt_covariances = !show_ndt_covariances;
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
	std::cout << "-: current_scan_index--" << std::endl;
	std::cout << "=: current_scan_index++" << std::endl;
	std::cout << "awsdzx: move current_scan (green)" << std::endl;
	std::cout << "p: print poses" << std::endl;
	std::cout << "t: optimize (Tait-Bryan)" << std::endl;
	std::cout << "l: optimize (Lie algebra)" << std::endl;
	std::cout << "1: sradius -= 0.01" << std::endl;
	std::cout << "2: sradius += 0.01" << std::endl;
	std::cout << "3: save current point clouds" << std::endl;
	std::cout << "g: show_ground_truth =! show_ground_truth" << std::endl;
	std::cout << "o: optimize (NDT)" << std::endl;
	std::cout << "n: nns" << std::endl;
	std::cout << "4: show ndt covariances on/off" << std::endl;
}

void set_initial_guess(std::vector<ScanPose>& scan_poses){
	scan_poses[0].m(0,0) = 0.380925;
	scan_poses[0].m(0,1) = 0.924606;
	scan_poses[0].m(0,3) = -6.06303;

	scan_poses[0].m(1,0) = -0.924606;
	scan_poses[0].m(1,1) = 0.380925;
	scan_poses[0].m(1,3) = -8.2396;

	scan_poses[1].m(0,0) = 0.710913;
	scan_poses[1].m(0,1) = 0.70328;
	scan_poses[1].m(0,3) = -1.03452;

	scan_poses[1].m(1,0) = -0.70328;
	scan_poses[1].m(1,1) = 0.710913;
	scan_poses[1].m(1,3) = -7.77656;

	scan_poses[2].m(0,0) = 0.613745;
	scan_poses[2].m(0,1) = 0.789504;
	scan_poses[2].m(0,3) = 3.99851;

	scan_poses[2].m(1,0) = -0.789504;
	scan_poses[2].m(1,1) = 0.613745;
	scan_poses[2].m(1,3) = -8.18709;

	scan_poses[3].m(0,0) = 0.751805;
	scan_poses[3].m(0,1) = 0.659385;
	scan_poses[3].m(0,3) = 9.99224;

	scan_poses[3].m(1,0) = -0.659385;
	scan_poses[3].m(1,1) = 0.751805;
	scan_poses[3].m(1,3) = -8.54882;

	scan_poses[4].m(0,0) = -0.996673;
	scan_poses[4].m(0,1) = 0.0815022;
	scan_poses[4].m(0,3) = 12.6963;

	scan_poses[4].m(1,0) = -0.0815022;
	scan_poses[4].m(1,1) = -0.996673;
	scan_poses[4].m(1,3) = -8.73861;

	scan_poses[5].m(0,0) = 0.639602;
	scan_poses[5].m(0,1) = 0.768705;
	scan_poses[5].m(0,3) = 17.8592;

	scan_poses[5].m(1,0) = -0.768705;
	scan_poses[5].m(1,1) = 0.639602;
	scan_poses[5].m(1,3) = -2.78423;

	scan_poses[6].m(0,0) = 0.745174;
	scan_poses[6].m(0,1) = 0.66687;
	scan_poses[6].m(0,3) = 22.9195;

	scan_poses[6].m(1,0) = -0.66687;
	scan_poses[6].m(1,1) = 0.745174;
	scan_poses[6].m(1,3) = -1.98392;

	scan_poses[7].m(0,0) = -0.0491839;
	scan_poses[7].m(0,1) = 0.998789;
	scan_poses[7].m(0,3) = 31.7827;

	scan_poses[7].m(1,0) = -0.998789;
	scan_poses[7].m(1,1) = -0.0491839;
	scan_poses[7].m(1,3) = -2.2143;

	scan_poses[8].m(0,0) = -0.128844;
	scan_poses[8].m(0,1) = 0.991665;
	scan_poses[8].m(0,3) = 39.0272;

	scan_poses[8].m(1,0) = -0.991665;
	scan_poses[8].m(1,1) = -0.128844;
	scan_poses[8].m(1,3) = -2.29705;

	scan_poses[9].m(0,0) = -0.34215;
	scan_poses[9].m(0,1) = 0.939646;
	scan_poses[9].m(0,3) = 48.1018;

	scan_poses[9].m(1,0) = -0.939646;
	scan_poses[9].m(1,1) = -0.34215;
	scan_poses[9].m(1,3) = -1.94245;

	scan_poses[10].m(0,0) = -0.158532;
	scan_poses[10].m(0,1) = 0.987354;
	scan_poses[10].m(0,3) = 54.2044;

	scan_poses[10].m(1,0) = -0.987354;
	scan_poses[10].m(1,1) = -0.158532;
	scan_poses[10].m(1,3) = -7.96743;

	scan_poses[11].m(0,0) = -0.197888;
	scan_poses[11].m(0,1) = 0.980225;
	scan_poses[11].m(0,3) = 65.5777;

	scan_poses[11].m(1,0) = -0.980225;
	scan_poses[11].m(1,1) = -0.197888;
	scan_poses[11].m(1,3) = -8.39231;

	scan_poses[12].m(0,0) = -0.360872;
	scan_poses[12].m(0,1) = 0.932615;
	scan_poses[12].m(0,3) = 78.1712;

	scan_poses[12].m(1,0) = -0.932615;
	scan_poses[12].m(1,1) = -0.360872;
	scan_poses[12].m(1,3) = -7.76261;
}

std::vector<std::pair<int,int>> nns(ScanPose &sp1, ScanPose &sp2, float radius)
{
	pcl::PointCloud<pcl::PointXYZ> pc1;
	pcl::PointCloud<pcl::PointXYZ> pc2;

	for(size_t i = 0; i < sp1.pc.size(); i++){
		Eigen::Vector3d v(sp1.pc[i].x, sp1.pc[i].y, sp1.pc[i].z);
		Eigen::Vector3d vt = sp1.m * v;
		pc1.push_back(pcl::PointXYZ(vt.x(), vt.y(), vt.z()));
	}

	for(size_t i = 0 ; i < sp2.pc.size(); i++){
		Eigen::Vector3d v(sp2.pc[i].x, sp2.pc[i].y, sp2.pc[i].z);
		Eigen::Vector3d vt = sp2.m * v;
		pc2.push_back(pcl::PointXYZ(vt.x(), vt.y(), vt.z()));
	}


	std::vector<std::pair<int,int>> result;

	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	for(size_t i = 0; i < pc1.size(); i++){
		cloud->push_back(pc1[i]);
	}
	int K = 1;
	
	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;

	kdtree.setInputCloud (cloud);
	for(size_t k = 0; k < pc2.size(); k++){
		if ( kdtree.radiusSearch (pc2[k], radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0 ){
			for (std::size_t i = 0; i < pointIdxRadiusSearch.size (); ++i){
				result.emplace_back(pointIdxRadiusSearch[i], k);
				break;
			}
		}
	}

	return result;
}

