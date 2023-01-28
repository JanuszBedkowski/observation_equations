#ifndef _STRUCTURES_H_
#define _STRUCTURES_H_

struct PerspectiveCameraParams{
	double cx;
	double cy;
	double fx;
	double fy;
};

struct MetricCameraParams{
	double c;
	double ksi_0;
	double eta_0;
};

struct TaitBryanPose
{
	double px;
	double py;
	double pz;
	double om;
	double fi;
	double ka;
};

struct RodriguesPose
{
	double px;
	double py;
	double pz;
	double sx;
	double sy;
	double sz;
};

struct QuaternionPose
{
	double px;
	double py;
	double pz;
	double q0;
	double q1;
	double q2;
	double q3;
};

struct Point3D
{
	double x;
	double y;
	double z;
	int index_pose;
};

#endif
