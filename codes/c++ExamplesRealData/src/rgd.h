#ifndef _RGD_H_
#define _RGD_H_

#include <Eigen/Eigen>

#include "structures.h"

struct GridParameters{
    double bounding_box_min_X;
    double bounding_box_min_Y;
    double bounding_box_min_Z;
    double bounding_box_max_X;
    double bounding_box_max_Y;
    double bounding_box_max_Z;
    double bounding_box_extension;
    int number_of_buckets_X;
    int number_of_buckets_Y;
    int number_of_buckets_Z;
    long long unsigned int number_of_buckets;
    double resolution_X;
    double resolution_Y;
    double resolution_Z;
};

struct PointBucketIndexPair{
    int index_of_point;
    long long unsigned int index_of_bucket;
    int index_pose;
};

struct Bucket{
	long long unsigned int index_begin;
	long long unsigned int index_end;
	long long unsigned int number_of_points;
	Eigen::Vector3d mean;
	Eigen::Matrix3d cov;
};

struct Job{
	long long unsigned int index_begin_inclusive;
	long long unsigned int index_end_exclusive;
};

void grid_calculate_params(const std::vector<Point3D>& point_cloud_global, GridParameters &in_out_params);
void build_rgd(std::vector<Point3D> &points, std::vector<PointBucketIndexPair>& index_pair, std::vector<Bucket>& buckets, GridParameters& rgd_params, int num_threads = 8);
std::vector<Job> get_jobs(long long unsigned int size, int num_threads = 8);
#endif
