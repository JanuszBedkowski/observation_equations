#include <GL/freeglut.h>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "../../metric_camera_colinearity_tait_bryan_wc_jacobian.h"
#include "../../metric_camera_colinearity_rodrigues_wc_jacobian.h"
#include "../../c++Examples/include/structures.h"
#include "../../c++Examples/include/transformations.h"
#include "../../ray_intersection_observation_equation_jacobian.h"
#include "../../c++Examples/include/cauchy.h"
#include "../../constraint_fixed_parameter_jacobian.h"

const unsigned int window_width = 1920;
const unsigned int window_height = 1080;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = -50.7999;
float rotate_y = 170.8;
float translate_z = -4345;
float translate_x = 1560;
float translate_y = -410;

bool initGL(int *argc, char **argv);
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void printHelp();

struct MetricCameraIntrinsicParameters{
	double ksi_0;
	double eta_0;
	double c;
};

struct Camera{
	Eigen::Affine3d pose;
};

struct ColinearityObservation{
	double ksi;
	double eta;
	int index_camera;
	int index_tie_point;
};

struct TiePoint{
	bool is_control_point;
	Eigen::Vector3d coordinates;
};

MetricCameraIntrinsicParameters cam_params;
std::vector<ColinearityObservation> colinearity_observations;
std::vector<Camera> cameras;
std::vector<TiePoint> tie_points;

double offset_x = 307442.51206;
double offset_y = 645119.05660;
double offset_z = 223.10;
Eigen::Vector3d inters;

Eigen::Vector3d bundle_of_rays_intersection(
		MetricCameraIntrinsicParameters cam_params,
		std::vector<ColinearityObservation> colinearity_observations,
		std::vector<Camera> cameras){
	Eigen::Vector3d intersection(0.1,0.1,0.1);

	for(int iter = 0 ; iter < 10; iter++){
		std::vector<Eigen::Triplet<double>> tripletListA;
		std::vector<Eigen::Triplet<double>> tripletListP;
		std::vector<Eigen::Triplet<double>> tripletListB;

		for(size_t i = 0; i < colinearity_observations.size(); i++){
			Eigen::Vector3d p(colinearity_observations[i].ksi - cam_params.ksi_0, colinearity_observations[i].eta - cam_params.eta_0, -cam_params.c);
			Eigen::Vector3d pt = cameras[colinearity_observations[i].index_camera].pose.rotation() * p;

			Eigen::Vector3d vx = pt.cross(Eigen::Vector3d (1,0,0));
			Eigen::Vector3d vy = pt.cross(vx);

			Eigen::Matrix<double, 2, 1> delta;
			ray_intersection_observation_equation(delta,
					intersection.x(), intersection.y(), intersection.z(),
					cameras[colinearity_observations[i].index_camera].pose(0,3), cameras[colinearity_observations[i].index_camera].pose(1,3), cameras[colinearity_observations[i].index_camera].pose(2,3),
					vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z());
			Eigen::Matrix<double, 2, 3> delta_jacobian;
			ray_intersection_observation_equation_jacobian(delta_jacobian,
					intersection.x(), intersection.y(), intersection.z(),
					cameras[colinearity_observations[i].index_camera].pose(0,3), cameras[colinearity_observations[i].index_camera].pose(1,3), cameras[colinearity_observations[i].index_camera].pose(2,3),
					vx.x(), vx.y(), vx.z(), vy.x(), vy.y(), vy.z());

			int ir = tripletListB.size();
			tripletListA.emplace_back(ir, 0, -delta_jacobian(0,0));
			tripletListA.emplace_back(ir, 1, -delta_jacobian(0,1));
			tripletListA.emplace_back(ir, 2, -delta_jacobian(0,2));
			tripletListA.emplace_back(ir + 1, 0, -delta_jacobian(1,0));
			tripletListA.emplace_back(ir + 1, 1, -delta_jacobian(1,1));
			tripletListA.emplace_back(ir + 1, 2, -delta_jacobian(1,2));

			tripletListP.emplace_back(ir    , ir    ,  1);
			tripletListP.emplace_back(ir + 1, ir + 1,  1);

			tripletListB.emplace_back(ir    , 0,  delta(0,0));
			tripletListB.emplace_back(ir + 1, 0,  delta(1,0));
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

		if(h_x.size() == 3){
			for(size_t i = 0 ; i < h_x.size(); i++){
				std::cout << h_x[i] << std::endl;
			}
			std::cout << "AtPA=AtPB SOLVED" << std::endl;
			std::cout << "update" << std::endl;
			intersection.x() += h_x[0];
			intersection.y() += h_x[1];
			intersection.z() += h_x[2];
		}else{
			std::cout << "AtPA=AtPB FAILED" << std::endl;
		}
	}

	return intersection;
}

int main(int argc, char *argv[]){
	cam_params.ksi_0 = 0;
	cam_params.eta_0 = 0;
	cam_params.c = 10000;

	//30056 0
	Camera cam;
	cam.pose = Eigen::Affine3d::Identity();
	cam.pose(0,0) = -0.624666076092;
	cam.pose(1,0) = -0.780886032111;
	cam.pose(2,0) =  0.003049300557;

	cam.pose(0,1) =  0.780890878954;
	cam.pose(1,1) = -0.624667002649;
	cam.pose(2,1) =  0.000755623701;

	cam.pose(0,2) =  0.001314741446;
	cam.pose(1,2) =  0.002853183485;
	cam.pose(2,2) =  0.999995065387;

	cam.pose(1,3) =  645119.05660;
	cam.pose(0,3) =  307442.51206;
	cam.pose(2,3) =  2193.61072;
	cameras.push_back(cam);

	//30055 1
	cam.pose = Eigen::Affine3d::Identity();
	cam.pose(0,0) = -0.634915760062;
	cam.pose(1,0) = -0.772493945357;
	cam.pose(2,0) =  -0.011622478706;

	cam.pose(0,1) =  0.772572805540;
	cam.pose(1,1) = -0.634907717713;
	cam.pose(2,1) =  -0.004842533230;

	cam.pose(0,2) =  -0.003638373829;
	cam.pose(1,2) =  -0.012053811648;
	cam.pose(2,2) =  0.999920730789;

	cam.pose(1,3) =  644710.53021;
	cam.pose(0,3) =  307898.39698;
	cam.pose(2,3) =  2193.27122;
	cameras.push_back(cam);

	//30054 2
	cam.pose = Eigen::Affine3d::Identity();
	cam.pose(0,0) = -0.652006792301;
	cam.pose(1,0) = -0.756112836475;
	cam.pose(2,0) =  -0.056396110778;

	cam.pose(0,1) =  0.757597712077;
	cam.pose(1,1) = -0.652669168106;
	cam.pose(2,1) = -0.008286353851;

	cam.pose(0,2) =  -0.030542584191 ;
	cam.pose(1,2) =  -0.048128323490;
	cam.pose(2,2) =  0.998374085716;

	cam.pose(1,3) =  644305.76605;
	cam.pose(0,3) =  308362.60224;
	cam.pose(2,3) =  2189.45394;
	cameras.push_back(cam);

	//30053 3
	cam.pose = Eigen::Affine3d::Identity();
	cam.pose(0,0) = -0.674360242193;
	cam.pose(1,0) = -0.738139907585;
	cam.pose(2,0) =  -0.019691129447;

	cam.pose(0,1) =  0.738395320937;
	cam.pose(1,1) = -0.674231513386;
	cam.pose(2,1) = -0.013572633300;

	cam.pose(0,2) =  -0.003257877718 ;
	cam.pose(1,2) =  -0.023692682127;
	cam.pose(2,2) =  0.999713980620;

	cam.pose(1,3) =  643895.92663;
	cam.pose(0,3) =  308809.51927;
	cam.pose(2,3) =  2186.58433;
	cameras.push_back(cam);

	//20049 4
	cam.pose = Eigen::Affine3d::Identity();
	cam.pose(0,0) = -0.681968279802;
	cam.pose(1,0) = -0.731323293151;
	cam.pose(2,0) =  -0.009246958318;

	cam.pose(0,1) =  0.731330032327 ;
	cam.pose(1,1) = -0.682015904523;
	cam.pose(2,1) = 0.003269525247;

	cam.pose(0,2) =  -0.008697652611;
	cam.pose(1,2) =  -0.004532865817;
	cam.pose(2,2) =  0.999951900827 ;

	cam.pose(1,3) =  646572.85818;
	cam.pose(0,3) =  308693.56340;
	cam.pose(2,3) =  2193.17111;
	cameras.push_back(cam);

	//20050 5
	cam.pose = Eigen::Affine3d::Identity();
	cam.pose(0,0) = -0.682691500788 ;
	cam.pose(1,0) = -0.730529323495;
	cam.pose(2,0) =  0.016100380945;

	cam.pose(0,1) =  0.730593360233 ;
	cam.pose(1,1) = -0.682807980473;
	cam.pose(2,1) = -0.002569783183 ;

	cam.pose(0,2) =  0.012870770569;
	cam.pose(1,2) = 0.010008462278;
	cam.pose(2,2) =   0.999867078140 ;

	cam.pose(1,3) =  646165.77814;
	cam.pose(0,3) =   309153.36786;
	cam.pose(2,3) =   2190.80044;
	cameras.push_back(cam);

	//20051 6
	cam.pose = Eigen::Affine3d::Identity();
	cam.pose(0,0) = -0.681043157946 ;
	cam.pose(1,0) = -0.732074427290 ;
	cam.pose(2,0) =   0.015724182713 ;

	cam.pose(0,1) =  0.732125432608;
	cam.pose(1,1) = -0.681161780346;
	cam.pose(2,1) = -0.003313596919  ;

	cam.pose(0,2) =  0.013136511858;
	cam.pose(1,2) = 0.009255371562;
	cam.pose(2,2) =   0.999870876740  ;

	cam.pose(1,3) =  645763.45239;
	cam.pose(0,3) =   309620.77521;
	cam.pose(2,3) =   2188.92747;
	cameras.push_back(cam);

	//20052 7
	cam.pose = Eigen::Affine3d::Identity();
	cam.pose(0,0) = -0.681855577878  ;
	cam.pose(1,0) = -0.731463292215 ;
	cam.pose(2,0) =   -0.005867116749  ;

	cam.pose(0,1) =  0.731472377219 ;
	cam.pose(1,1) = -0.681870586424;
	cam.pose(2,1) = 0.000815313520   ;

	cam.pose(0,2) =  -0.004596986250;
	cam.pose(1,2) = -0.003735707765;
	cam.pose(2,2) =   0.999982455949   ;

	cam.pose(1,3) = 645356.06683;
	cam.pose(0,3) =   310081.03728;
	cam.pose(2,3) =   2190.84287;
	cameras.push_back(cam);


	//controlpoints
	TiePoint tp;
	tp.is_control_point = true;
	tp.coordinates.x() = 308084.47;
	tp.coordinates.y() = 645766.50;
	tp.coordinates.z() = 223.10;
	tie_points.push_back(tp);

	tp.coordinates.x() = 306849.35;
	tp.coordinates.y() = 644397.25;
	tp.coordinates.z() = 255.23;
	tie_points.push_back(tp);

	tp.coordinates.x() = 307756.61;
	tp.coordinates.y() = 643682.93;
	tp.coordinates.z() = 255.70;
	tie_points.push_back(tp);

	tp.coordinates.x() = 309337.60;
	tp.coordinates.y() = 647227.36;
	tp.coordinates.z() = 258.59;
	tie_points.push_back(tp);

	tp.coordinates.x() = 310237.72;
	tp.coordinates.y() = 646409.65;
	tp.coordinates.z() = 241.63;
	tie_points.push_back(tp);






	/////////////////////////////////////////////////////////////////
	//301	484.75	-4637.83	308084.47	645766.50	223.10	30056
	ColinearityObservation obs;
	obs.index_tie_point = 0;
	obs.eta = 484.75;
	obs.ksi = -4637.83;
	obs.index_camera = 0;
	colinearity_observations.push_back(obs);

	//3561	-70.92	4781.83	306849.35	644397.25	255.23	30056
	obs.index_tie_point = 1;
	obs.eta = -70.92;
	obs.ksi = 4781.83;
	obs.index_camera = 0;
	colinearity_observations.push_back(obs);

	//301	-2607.00	-4593.17	308084.47	645766.50	223.10	30055
	obs.index_tie_point = 0;
	obs.eta = -2607.00;
	obs.ksi = -4593.17;
	obs.index_camera = 1;
	colinearity_observations.push_back(obs);

	//3541	2869.50	4710.00	307756.61	643682.93	255.70	30055
	obs.index_tie_point = 2;
	obs.eta = 2869.50;
	obs.ksi = 4710.00;
	obs.index_camera = 1;
	colinearity_observations.push_back(obs);

	//3561	-3119.50	4820.83	306849.35	644397.25	255.23	30055
	obs.index_tie_point = 1;
	obs.eta = -3119.50;
	obs.ksi = 4820.83;
	obs.index_camera = 1;
	colinearity_observations.push_back(obs);

	//3541	-194.33	5180.92	307756.61	643682.93	255.70	30054
	obs.index_tie_point = 2;
	obs.eta = -194.33;
	obs.ksi = 5180.92;
	obs.index_camera = 2;
	colinearity_observations.push_back(obs);

	//3541	-3161.92	4710.42	307756.61	643682.93	255.70	30053
	obs.index_tie_point = 2;
	obs.eta = -3161.92;
	obs.ksi = 4710.42;
	obs.index_camera = 3;
	colinearity_observations.push_back(obs);

	//2492	94.25	-4632.08	309337.60	647227.36	258.59	20049
	obs.index_tie_point = 3;
	obs.eta = 94.25;
	obs.ksi = -4632.08;
	obs.index_camera = 4;
	colinearity_observations.push_back(obs);

	//301	499.83	5218.17	308084.47	645766.50	223.10	20049
	obs.index_tie_point = 0;
	obs.eta = 499.83;
	obs.ksi = 5218.17;
	obs.index_camera = 4;
	colinearity_observations.push_back(obs);

	//2492	-3049.83	-4859.33	309337.60	647227.36	258.59	20050
	obs.index_tie_point = 3;
	obs.eta = -3049.83;
	obs.ksi = -4859.33;
	obs.index_camera = 5;
	colinearity_observations.push_back(obs);

	//2512	3263.08	-4914.75	310237.72	646409.65	241.63	20050
	obs.index_tie_point = 4;
	obs.eta = 3263.08;
	obs.ksi = -4914.75;
	obs.index_camera = 5;
	colinearity_observations.push_back(obs);

	//301	-2534.92	4985.42	308084.47	645766.50	223.10	20050
	obs.index_tie_point = 0;
	obs.eta = -2534.92;
	obs.ksi = 4985.42;
	obs.index_camera = 5;
	colinearity_observations.push_back(obs);

	/*//2512	92.75	-4779.42	309337.60	647227.36	258.59	20051
	obs.tie_point.x() = 309337.60;
	obs.tie_point.y() = 647227.36;
	obs.tie_point.z() = 258.59;
	obs.eta = 92.75;
	obs.ksi = -4779.42;
	obs.index_to_camera = 6;
	colinearity_observations.push_back(obs);*/

	//2512	-3098.67	-4432.83	310237.72	646409.65	241.63	20052
	obs.index_tie_point = 4;
	obs.eta = -3098.67;
	obs.ksi = -4432.83;
	obs.index_camera = 7;
	colinearity_observations.push_back(obs);

	for(size_t i = 0 ; i < cameras.size(); i ++){
		cameras[i].pose(0,3) -= offset_x;
		cameras[i].pose(1,3) -= offset_y;
		cameras[i].pose(2,3) -= offset_z;
	}

	for(size_t i = 0 ; i < tie_points.size(); i++){
		tie_points[i].coordinates.x() -= offset_x;
		tie_points[i].coordinates.y() -= offset_y;
		tie_points[i].coordinates.z() -= offset_z;
	}


	//
	std::cout << "colinearity_observations.size(): " << colinearity_observations.size() << std::endl;

	std::vector<ColinearityObservation> colinearity_observations_manual;
///////////////////////////////////////////////////////////////////////////////////

	//30056 0
	//30055 1
	//30054 2
	//30053 3
	//20049 4
	//20050 5
	//20051 6
	//20052 7

	//30054
	obs.eta = 4175.5;
	obs.ksi = 12549.6;
	obs.index_camera = 2;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//30053
	obs.eta = 984.4;
	obs.ksi = 13037.6;
	obs.index_camera = 3;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//30052
	obs.eta = 1063.7;
	obs.ksi = 3310.2;
	obs.index_camera = 7;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//30051
	obs.eta = 4254.7;
	obs.ksi = 3628.8;
	obs.index_camera = 6;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//20050
	obs.eta = 7398.4;
	obs.ksi = 3761.9;
	obs.index_camera = 5;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations_manual[i].eta = colinearity_observations_manual[i].eta - 7680*0.5;
		colinearity_observations_manual[i].ksi = -colinearity_observations_manual[i].ksi + 13824*0.5;
	}

	inters = bundle_of_rays_intersection(
				cam_params,
				colinearity_observations_manual,
				cameras);

	tp.is_control_point = false;
	tp.coordinates = inters;
	tie_points.push_back(tp);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations.push_back(colinearity_observations_manual[i]);
	}
	colinearity_observations_manual.clear();
	///////////////////////////////////////////////////////////////////

	//2
	//
	//30054
	obs.eta = 4405.3;
	obs.ksi = 12372.6;
	obs.index_camera = 2;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//30052
	obs.eta = 1311.6;
	obs.ksi = 3127.6;
	obs.index_camera = 7;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//30051
	obs.eta = 4500.4;
	obs.ksi = 3448.4;
	obs.index_camera = 6;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//20050
	obs.eta = 7644.4;
	obs.ksi = 3582.4;
	obs.index_camera = 5;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//30055
	obs.eta = 7662.8;
	obs.ksi = 12860.8;
	obs.index_camera = 1;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations_manual[i].eta = colinearity_observations_manual[i].eta - 7680*0.5;
		colinearity_observations_manual[i].ksi = -colinearity_observations_manual[i].ksi + 13824*0.5;
	}

	inters = bundle_of_rays_intersection(
				cam_params,
				colinearity_observations_manual,
				cameras);

	tp.is_control_point = false;
	tp.coordinates = inters;
	tie_points.push_back(tp);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations.push_back(colinearity_observations_manual[i]);
	}
	colinearity_observations_manual.clear();

	/////////////////////

	//3
	//30056
	obs.eta = 4754.2;
	obs.ksi = 11142.0;
	obs.index_camera = 0;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//30055
	obs.eta = 1646.7;
	obs.ksi = 11109.8;
	obs.index_camera = 1;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//30049
	obs.eta = 4803.0;
	obs.ksi = 1265.6;
	obs.index_camera = 4;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//20050
	obs.eta = 1747.9;
	obs.ksi = 1507.2;
	obs.index_camera = 5;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations_manual[i].eta = colinearity_observations_manual[i].eta - 7680*0.5;
		colinearity_observations_manual[i].ksi = -colinearity_observations_manual[i].ksi + 13824*0.5;
	}

	inters = bundle_of_rays_intersection(
				cam_params,
				colinearity_observations_manual,
				cameras);
	tp.is_control_point = false;
	tp.coordinates = inters;
	tie_points.push_back(tp);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations.push_back(colinearity_observations_manual[i]);
	}
	colinearity_observations_manual.clear();

	//5
	//30049
	obs.eta = 5831.6;
	obs.ksi = 5690.2;
	obs.index_camera = 4;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//20050
	obs.eta = 2785.8;
	obs.ksi = 5863.4;
	obs.index_camera = 5;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations_manual[i].eta = colinearity_observations_manual[i].eta - 7680*0.5;
		colinearity_observations_manual[i].ksi = -colinearity_observations_manual[i].ksi + 13824*0.5;
	}

	inters = bundle_of_rays_intersection(
				cam_params,
				colinearity_observations_manual,
				cameras);
	tp.is_control_point = false;
	tp.coordinates = inters;
	tie_points.push_back(tp);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations.push_back(colinearity_observations_manual[i]);
	}
	colinearity_observations_manual.clear();

	//6
	//30051
	obs.eta = 5460.6;
	obs.ksi = 13278.1;
	obs.index_camera = 6;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//30052
	obs.eta = 2257.2;
	obs.ksi = 12891.2;
	obs.index_camera = 7;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations_manual[i].eta = colinearity_observations_manual[i].eta - 7680*0.5;
		colinearity_observations_manual[i].ksi = -colinearity_observations_manual[i].ksi + 13824*0.5;
	}

	inters = bundle_of_rays_intersection(
				cam_params,
				colinearity_observations_manual,
				cameras);
	tp.is_control_point = false;
	tp.coordinates = inters;
	tie_points.push_back(tp);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations.push_back(colinearity_observations_manual[i]);
	}
	colinearity_observations_manual.clear();

	//7
	//30056
	obs.eta = 5530.7;
	obs.ksi = 6720.0;
	obs.index_camera = 0;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//30055
	obs.eta = 2444.1;
	obs.ksi = 6730.1;
	obs.index_camera = 1;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations_manual[i].eta = colinearity_observations_manual[i].eta - 7680*0.5;
		colinearity_observations_manual[i].ksi = -colinearity_observations_manual[i].ksi + 13824*0.5;
	}

	inters = bundle_of_rays_intersection(
				cam_params,
				colinearity_observations_manual,
				cameras);
	tp.is_control_point = false;
	tp.coordinates = inters;
	tie_points.push_back(tp);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations.push_back(colinearity_observations_manual[i]);
	}
	colinearity_observations_manual.clear();


	///////////////for students///////////////////////
	//30056 0
	//30055 1
	//30054 2
	//30053 3
	//20049 4
	//20050 5
	//20051 6
	//20052 7

	//8
	//20049
	obs.eta = 6272.0;
	obs.ksi = 3725.3;
	obs.index_camera = 4;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	//20050
	obs.eta = 3238.0;
	obs.ksi = 3916.5;
	obs.index_camera = 5;
	obs.index_tie_point = tie_points.size();
	colinearity_observations_manual.push_back(obs);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations_manual[i].eta = colinearity_observations_manual[i].eta - 7680*0.5;
		colinearity_observations_manual[i].ksi = -colinearity_observations_manual[i].ksi + 13824*0.5;
	}

	inters = bundle_of_rays_intersection(
				cam_params,
				colinearity_observations_manual,
				cameras);
	tp.is_control_point = false;
	tp.coordinates = inters;
	tie_points.push_back(tp);

	for(size_t i = 0 ; i < colinearity_observations_manual.size(); i++){
		colinearity_observations.push_back(colinearity_observations_manual[i]);
	}
	colinearity_observations_manual.clear();
	///////////////for students///////////////////////



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
	glutCreateWindow("photogrammetry");
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
	for(size_t i = 0 ; i < cameras.size(); i ++){
		float scale = 100;

		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(cameras[i].pose(0,3), cameras[i].pose(1,3), cameras[i].pose(2,3));
		glVertex3f(cameras[i].pose(0,3) + cameras[i].pose(0,0) * scale, cameras[i].pose(1,3) + cameras[i].pose(0,1) * scale, cameras[i].pose(2,3) + cameras[i].pose(0,2) * scale);

		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(cameras[i].pose(0,3), cameras[i].pose(1,3), cameras[i].pose(2,3));
		glVertex3f(cameras[i].pose(0,3) + cameras[i].pose(1,0) * scale, cameras[i].pose(1,3) + cameras[i].pose(1,1) * scale, cameras[i].pose(2,3) + cameras[i].pose(1,2) * scale);

		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(cameras[i].pose(0,3), cameras[i].pose(1,3), cameras[i].pose(2,3));
		glVertex3f(cameras[i].pose(0,3) + cameras[i].pose(2,0) * scale, cameras[i].pose(1,3) + cameras[i].pose(2,1) * scale, cameras[i].pose(2,3) + cameras[i].pose(2,2) * scale);

	}
	glEnd();

	glColor3f(0,1,0);
	glBegin(GL_LINES);
	for(size_t i = 0 ; i < colinearity_observations.size(); i++){
		if(tie_points[colinearity_observations[i].index_tie_point].is_control_point){
			glColor3f(0,1,0);
		}else{
			glColor3f(0.5,0.5,0.5);
		}

		glVertex3f(cameras[colinearity_observations[i].index_camera].pose(0,3), cameras[colinearity_observations[i].index_camera].pose(1,3), cameras[colinearity_observations[i].index_camera].pose(2,3));
		glVertex3f(tie_points[colinearity_observations[i].index_tie_point].coordinates.x(), tie_points[colinearity_observations[i].index_tie_point].coordinates.y(), tie_points[colinearity_observations[i].index_tie_point].coordinates.z());
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
		case 't':{
			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < colinearity_observations.size(); i++){

				TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[colinearity_observations[i].index_camera].pose);

				Eigen::Matrix<double, 2, 1> delta;
				observation_equation_metric_camera_colinearity_tait_bryan_wc(delta, cam_params.ksi_0, cam_params.eta_0, cam_params.c, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
						tie_points[colinearity_observations[i].index_tie_point].coordinates.x(), tie_points[colinearity_observations[i].index_tie_point].coordinates.y(), tie_points[colinearity_observations[i].index_tie_point].coordinates.z(),
						colinearity_observations[i].ksi, colinearity_observations[i].eta);

				Eigen::Matrix<double, 2, 9, Eigen::RowMajor> jacobian;
				observation_equation_metric_camera_colinearity_tait_bryan_wc_jacobian(jacobian, cam_params.ksi_0, cam_params.eta_0, cam_params.c, pose.px, pose.py, pose.pz, pose.om, pose.fi, pose.ka,
						tie_points[colinearity_observations[i].index_tie_point].coordinates.x(), tie_points[colinearity_observations[i].index_tie_point].coordinates.y(), tie_points[colinearity_observations[i].index_tie_point].coordinates.z(),
						colinearity_observations[i].ksi, colinearity_observations[i].eta);

				std::cout << delta(0,0) << "," << delta(1,0) << "|";

				int ir = tripletListB.size();
				int ic_camera = colinearity_observations[i].index_camera * 6;
				int ic_tie_point = cameras.size() * 6 + colinearity_observations[i].index_tie_point * 3;

				tripletListA.emplace_back(ir     , ic_camera    , -jacobian(0,0));
				tripletListA.emplace_back(ir     , ic_camera + 1, -jacobian(0,1));
				tripletListA.emplace_back(ir     , ic_camera + 2, -jacobian(0,2));
				tripletListA.emplace_back(ir     , ic_camera + 3, -jacobian(0,3));
				tripletListA.emplace_back(ir     , ic_camera + 4, -jacobian(0,4));
				tripletListA.emplace_back(ir     , ic_camera + 5, -jacobian(0,5));

				tripletListA.emplace_back(ir     , ic_tie_point,     -jacobian(0,6));
				tripletListA.emplace_back(ir     , ic_tie_point + 1, -jacobian(0,7));
				tripletListA.emplace_back(ir     , ic_tie_point + 2, -jacobian(0,8));

				tripletListA.emplace_back(ir + 1 , ic_camera    , -jacobian(1,0));
				tripletListA.emplace_back(ir + 1 , ic_camera + 1, -jacobian(1,1));
				tripletListA.emplace_back(ir + 1 , ic_camera + 2, -jacobian(1,2));
				tripletListA.emplace_back(ir + 1 , ic_camera + 3, -jacobian(1,3));
				tripletListA.emplace_back(ir + 1 , ic_camera + 4, -jacobian(1,4));
				tripletListA.emplace_back(ir + 1 , ic_camera + 5, -jacobian(1,5));

				tripletListA.emplace_back(ir + 1 , ic_tie_point,     -jacobian(1,6));
				tripletListA.emplace_back(ir + 1 , ic_tie_point + 1, -jacobian(1,7));
				tripletListA.emplace_back(ir + 1 , ic_tie_point + 2, -jacobian(1,8));

				tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
				tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));

				tripletListB.emplace_back(ir    , 0,  delta(0,0));
				tripletListB.emplace_back(ir + 1, 0,  delta(1,0));

			}

			for(size_t i = 0; i < tie_points.size(); i++){
				if(!tie_points[i].is_control_point){
					continue;
				}

				int ir = tripletListB.size();
				int ic = cameras.size() * 6 + i * 3;

				double residual;
				residual_constraint_fixed_optimization_parameter(residual, tie_points[i].coordinates.x(), tie_points[i].coordinates.x());
				Eigen::Matrix<double, 1, 1> jacobian;
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, tie_points[i].coordinates.x(), tie_points[i].coordinates.x());

				tripletListA.emplace_back(ir, ic  , -jacobian(0,0));
				tripletListB.emplace_back(ir, 0, residual);
				tripletListP.emplace_back(ir, ir, 1000000000);

				residual_constraint_fixed_optimization_parameter(residual, tie_points[i].coordinates.y(), tie_points[i].coordinates.y());
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, tie_points[i].coordinates.y(), tie_points[i].coordinates.y());

				tripletListA.emplace_back(ir + 1, ic + 1 , -jacobian(0,0));
				tripletListB.emplace_back(ir + 1, 0, residual);
				tripletListP.emplace_back(ir + 1, ir + 1, 1000000000);

				residual_constraint_fixed_optimization_parameter(residual, tie_points[i].coordinates.z(), tie_points[i].coordinates.z());
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, tie_points[i].coordinates.z(), tie_points[i].coordinates.z());

				tripletListA.emplace_back(ir + 2, ic + 2, -jacobian(0,0));
				tripletListB.emplace_back(ir + 2, 0, residual);
				tripletListP.emplace_back(ir + 2, ir + 2, 1000000000);
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 6 + tie_points.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 6 + tie_points.size() * 3, cameras.size() * 6 + tie_points.size() * 3);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 6 + tie_points.size() * 3, 1);

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

			if(h_x.size() == cameras.size() * 6 + tie_points.size() * 3){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < cameras.size(); i++){
					TaitBryanPose pose = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.om += h_x[counter++];
					pose.fi += h_x[counter++];
					pose.ka += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_tait_bryan(pose);
				}

				for(size_t i = 0; i < tie_points.size(); i++){
					tie_points[i].coordinates.x() += h_x[counter++];
					tie_points[i].coordinates.y() += h_x[counter++];
					tie_points[i].coordinates.z() += h_x[counter++];
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'r':{
			for(size_t i = 0; i < cameras.size(); i++){
				TaitBryanPose posetb = pose_tait_bryan_from_affine_matrix(cameras[i].pose);
				posetb.px += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00000001;
				posetb.py += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00000001;
				posetb.pz += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00000001;
				posetb.om += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00000001;
				posetb.fi += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00000001;
				posetb.ka += (float(rand()%1000000)/1000000.0 - 0.5) * 2.0 * 0.00000001;
				cameras[i].pose = affine_matrix_from_pose_tait_bryan(posetb);
			}

			std::vector<Eigen::Triplet<double>> tripletListA;
			std::vector<Eigen::Triplet<double>> tripletListP;
			std::vector<Eigen::Triplet<double>> tripletListB;

			for(size_t i = 0; i < colinearity_observations.size(); i++){

				RodriguesPose pose = pose_rodrigues_from_affine_matrix(cameras[colinearity_observations[i].index_camera].pose);

				Eigen::Matrix<double, 2, 1> delta;
				observation_equation_metric_camera_colinearity_rodrigues_wc(delta, cam_params.ksi_0, cam_params.eta_0, cam_params.c, pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
						tie_points[colinearity_observations[i].index_tie_point].coordinates.x(), tie_points[colinearity_observations[i].index_tie_point].coordinates.y(), tie_points[colinearity_observations[i].index_tie_point].coordinates.z(),
						colinearity_observations[i].ksi, colinearity_observations[i].eta);

				Eigen::Matrix<double, 2, 9, Eigen::RowMajor> jacobian;
				observation_equation_metric_camera_colinearity_rodrigues_wc_jacobian(jacobian, cam_params.ksi_0, cam_params.eta_0, cam_params.c, pose.px, pose.py, pose.pz, pose.sx, pose.sy, pose.sz,
						tie_points[colinearity_observations[i].index_tie_point].coordinates.x(), tie_points[colinearity_observations[i].index_tie_point].coordinates.y(), tie_points[colinearity_observations[i].index_tie_point].coordinates.z(),
						colinearity_observations[i].ksi, colinearity_observations[i].eta);

				std::cout << delta(0,0) << "," << delta(1,0) << "|";

				int ir = tripletListB.size();
				int ic_camera = colinearity_observations[i].index_camera * 6;
				int ic_tie_point = cameras.size() * 6 + colinearity_observations[i].index_tie_point * 3;

				tripletListA.emplace_back(ir     , ic_camera    , -jacobian(0,0));
				tripletListA.emplace_back(ir     , ic_camera + 1, -jacobian(0,1));
				tripletListA.emplace_back(ir     , ic_camera + 2, -jacobian(0,2));
				tripletListA.emplace_back(ir     , ic_camera + 3, -jacobian(0,3));
				tripletListA.emplace_back(ir     , ic_camera + 4, -jacobian(0,4));
				tripletListA.emplace_back(ir     , ic_camera + 5, -jacobian(0,5));

				tripletListA.emplace_back(ir     , ic_tie_point,     -jacobian(0,6));
				tripletListA.emplace_back(ir     , ic_tie_point + 1, -jacobian(0,7));
				tripletListA.emplace_back(ir     , ic_tie_point + 2, -jacobian(0,8));

				tripletListA.emplace_back(ir + 1 , ic_camera    , -jacobian(1,0));
				tripletListA.emplace_back(ir + 1 , ic_camera + 1, -jacobian(1,1));
				tripletListA.emplace_back(ir + 1 , ic_camera + 2, -jacobian(1,2));
				tripletListA.emplace_back(ir + 1 , ic_camera + 3, -jacobian(1,3));
				tripletListA.emplace_back(ir + 1 , ic_camera + 4, -jacobian(1,4));
				tripletListA.emplace_back(ir + 1 , ic_camera + 5, -jacobian(1,5));

				tripletListA.emplace_back(ir + 1 , ic_tie_point,     -jacobian(1,6));
				tripletListA.emplace_back(ir + 1 , ic_tie_point + 1, -jacobian(1,7));
				tripletListA.emplace_back(ir + 1 , ic_tie_point + 2, -jacobian(1,8));

				tripletListP.emplace_back(ir    , ir    ,  cauchy(delta(0,0), 1));
				tripletListP.emplace_back(ir + 1, ir + 1,  cauchy(delta(1,0), 1));

				tripletListB.emplace_back(ir    , 0,  delta(0,0));
				tripletListB.emplace_back(ir + 1, 0,  delta(1,0));

			}

			for(size_t i = 0; i < tie_points.size(); i++){
				if(!tie_points[i].is_control_point){
					continue;
				}

				int ir = tripletListB.size();
				int ic = cameras.size() * 6 + i * 3;

				double residual;
				residual_constraint_fixed_optimization_parameter(residual, tie_points[i].coordinates.x(), tie_points[i].coordinates.x());
				Eigen::Matrix<double, 1, 1> jacobian;
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, tie_points[i].coordinates.x(), tie_points[i].coordinates.x());

				tripletListA.emplace_back(ir, ic  , -jacobian(0,0));
				tripletListB.emplace_back(ir, 0, residual);
				tripletListP.emplace_back(ir, ir, 1000000000);

				residual_constraint_fixed_optimization_parameter(residual, tie_points[i].coordinates.y(), tie_points[i].coordinates.y());
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, tie_points[i].coordinates.y(), tie_points[i].coordinates.y());

				tripletListA.emplace_back(ir + 1, ic + 1 , -jacobian(0,0));
				tripletListB.emplace_back(ir + 1, 0, residual);
				tripletListP.emplace_back(ir + 1, ir + 1, 1000000000);

				residual_constraint_fixed_optimization_parameter(residual, tie_points[i].coordinates.z(), tie_points[i].coordinates.z());
				residual_constraint_fixed_optimization_parameter_jacobian(jacobian, tie_points[i].coordinates.z(), tie_points[i].coordinates.z());

				tripletListA.emplace_back(ir + 2, ic + 2, -jacobian(0,0));
				tripletListB.emplace_back(ir + 2, 0, residual);
				tripletListP.emplace_back(ir + 2, ir + 2, 1000000000);
			}

			Eigen::SparseMatrix<double> matA(tripletListB.size(), cameras.size() * 6 + tie_points.size() * 3);
			Eigen::SparseMatrix<double> matP(tripletListB.size(), tripletListB.size());
			Eigen::SparseMatrix<double> matB(tripletListB.size(), 1);

			matA.setFromTriplets(tripletListA.begin(), tripletListA.end());
			matP.setFromTriplets(tripletListP.begin(), tripletListP.end());
			matB.setFromTriplets(tripletListB.begin(), tripletListB.end());

			Eigen::SparseMatrix<double> AtPA(cameras.size() * 6 + tie_points.size() * 3, cameras.size() * 6 + tie_points.size() * 3);
			Eigen::SparseMatrix<double> AtPB(cameras.size() * 6 + tie_points.size() * 3, 1);

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

			if(h_x.size() == cameras.size() * 6 + tie_points.size() * 3){
				for(size_t i = 0 ; i < h_x.size(); i++){
					std::cout << h_x[i] << std::endl;
				}
				std::cout << "AtPA=AtPB SOLVED" << std::endl;
				std::cout << "update" << std::endl;

				int counter = 0;

				for(size_t i = 0; i < cameras.size(); i++){
					RodriguesPose pose = pose_rodrigues_from_affine_matrix(cameras[i].pose);
					pose.px += h_x[counter++];
					pose.py += h_x[counter++];
					pose.pz += h_x[counter++];
					pose.sx += h_x[counter++];
					pose.sy += h_x[counter++];
					pose.sz += h_x[counter++];

					cameras[i].pose = affine_matrix_from_pose_rodrigues(pose);
				}

				for(size_t i = 0; i < tie_points.size(); i++){
					tie_points[i].coordinates.x() += h_x[counter++];
					tie_points[i].coordinates.y() += h_x[counter++];
					tie_points[i].coordinates.z() += h_x[counter++];
				}
			}else{
				std::cout << "AtPA=AtPB FAILED" << std::endl;
			}
			break;
		}
		case 'n':{
			for(size_t i = 0; i < tie_points.size(); i++){
				if(tie_points[i].is_control_point){
					continue;
				}

				tie_points[i].coordinates.x() += (float(rand()%1000000)/1000000.0 - 0.5) * 10.0;
				tie_points[i].coordinates.y() += (float(rand()%1000000)/1000000.0 - 0.5) * 10.0;
				tie_points[i].coordinates.z() += (float(rand()%1000000)/1000000.0 - 0.5) * 10.0;
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
		translate_z += dy * 0.5f;
	} else if (mouse_buttons & 3) {
		translate_x += dx * 0.5f;
		translate_y -= dy * 0.5f;
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
	std::cout << "t: optimize Tait-Bryan" << std::endl;
	std::cout << "r: optimize Rodrigues" << std::endl;
	std::cout << "rotate_x: " << rotate_x << " rotate_y: " << rotate_y << " translate_x " << translate_x << " translate_y " << translate_y << " translate_z " << translate_z << std::endl;

}







