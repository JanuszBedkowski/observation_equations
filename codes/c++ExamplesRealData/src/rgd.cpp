#include <iostream>
#include <thread>
#include "rgd.h"

void grid_calculate_params(const std::vector<Point3D>& point_cloud_global, GridParameters &in_out_params)
{
	double min_x = std::numeric_limits<double>::max();
	double max_x = std::numeric_limits<double>::min();

	double min_y = std::numeric_limits<double>::max();
	double max_y = std::numeric_limits<double>::min();

	double min_z = std::numeric_limits<double>::max();
	double max_z = std::numeric_limits<double>::min();

	for(size_t i = 0 ; i < point_cloud_global.size(); i++) {
		if(point_cloud_global[i].x < min_x) min_x = point_cloud_global[i].x;
		if(point_cloud_global[i].x > max_x) max_x = point_cloud_global[i].x;

		if(point_cloud_global[i].y < min_y) min_y = point_cloud_global[i].y;
		if(point_cloud_global[i].y > max_y) max_y = point_cloud_global[i].y;

		if(point_cloud_global[i].z < min_z) min_z = point_cloud_global[i].z;
		if(point_cloud_global[i].z > max_z) max_z = point_cloud_global[i].z;
	}
	long long unsigned int number_of_buckets_X=((max_x - min_x)/in_out_params.resolution_X)+1;
	long long unsigned int number_of_buckets_Y=((max_y - min_y)/in_out_params.resolution_Y)+1;
	long long unsigned int number_of_buckets_Z=((max_z - min_z)/in_out_params.resolution_Z)+1;


	in_out_params.number_of_buckets_X = number_of_buckets_X;
	in_out_params.number_of_buckets_Y = number_of_buckets_Y;
	in_out_params.number_of_buckets_Z = number_of_buckets_Z;
	in_out_params.number_of_buckets   = static_cast<long long unsigned int>(number_of_buckets_X) *
		static_cast<long long unsigned int>(number_of_buckets_Y) * static_cast<long long unsigned int>(number_of_buckets_Z);

	in_out_params.bounding_box_max_X = max_x;
	in_out_params.bounding_box_min_X = min_x;
	in_out_params.bounding_box_max_Y = max_y;
	in_out_params.bounding_box_min_Y = min_y;
	in_out_params.bounding_box_max_Z = max_z;
	in_out_params.bounding_box_min_Z = min_z;
}

std::vector<Job> get_jobs(long long unsigned int size, int num_threads){

	int hc = size / num_threads;
	if(hc < 1)hc = 1;

	std::vector<Job> jobs;
	for(long long unsigned int i = 0 ; i < size; i += hc){
		long long unsigned int sequence_length = hc;
		if(i + hc >= size){
			sequence_length = size - i;
		}
		if(sequence_length == 0)break;

		Job j;
		j.index_begin_inclusive = i;
		j.index_end_exclusive = i + sequence_length;
		jobs.push_back(j);
	}

	//std::cout << jobs.size()<< " jobs; chunks: ";
	//for(size_t i = 0; i < jobs.size(); i++){
	//	std::cout << "("<<jobs[i].index_begin_inclusive << ", " << jobs[i].index_end_exclusive <<") ";
	//}
	//std::cout << "\n";
	return jobs;
}

void reindex_job(int i, Job* job, std::vector<Point3D>* points, std::vector<PointBucketIndexPair>* pairs, GridParameters rgd_params){
	for(size_t ii = job->index_begin_inclusive; ii < job->index_end_exclusive; ii++){

		Point3D& p = (*points)[ii];

		(*pairs)[ii].index_of_point = ii;
		(*pairs)[ii].index_of_bucket = 0;
		(*pairs)[ii].index_pose = p.index_pose;


		long long unsigned int ix = (p.x - rgd_params.bounding_box_min_X) / rgd_params.resolution_X;
		long long unsigned int iy = (p.y - rgd_params.bounding_box_min_Y) / rgd_params.resolution_Y;
		long long unsigned int iz = (p.z - rgd_params.bounding_box_min_Z) / rgd_params.resolution_Z;

		(*pairs)[ii].index_of_bucket = ix* static_cast<long long unsigned int>(rgd_params.number_of_buckets_Y) *
				static_cast<long long unsigned int>(rgd_params.number_of_buckets_Z) + iy * static_cast<long long unsigned int>( rgd_params.number_of_buckets_Z) + iz;
	}
}


void reindex(std::vector<Point3D> &points, std::vector<PointBucketIndexPair>& index_pair, GridParameters& rgd_params, int num_threads)
{
	index_pair.resize(points.size());

	std::vector<Job> jobs = get_jobs(index_pair.size(), num_threads);

	std::vector<std::thread> threads;

	for(size_t i = 0; i < jobs.size(); i++){
		threads.push_back(std::thread( reindex_job, i, &jobs[i], &points, &index_pair, rgd_params));
	}

	for(size_t j = 0; j < threads.size(); j++){
		threads[j].join();
	}
	threads.clear();

	std::sort(index_pair.begin(), index_pair.end(),  [] (const PointBucketIndexPair & a, const PointBucketIndexPair & b) { return (  (a.index_of_bucket == b.index_of_bucket) ? (a.index_pose < b.index_pose) : (a.index_of_bucket < b.index_of_bucket)    ) ;}  );
}

void build_rgd_init_job(int i, Job* job, std::vector<Bucket>* buckets){

	for(size_t ii = job->index_begin_inclusive; ii < job->index_end_exclusive; ii++){
		(*buckets)[ii].index_begin = -1;
		(*buckets)[ii].index_end = -1;
		(*buckets)[ii].number_of_points = 0;
	}
}

void build_rgd_job(int i, Job* job, std::vector<PointBucketIndexPair>* index_pair, std::vector<Bucket>* buckets){
	for(size_t ii = job->index_begin_inclusive; ii < job->index_end_exclusive; ii++){
		int ind = ii;
		if(ind == 0)
		{
			long long unsigned int index_of_bucket = (*index_pair)[ind].index_of_bucket;
			long long unsigned int index_of_bucket_1 = (*index_pair)[ind+1].index_of_bucket;

			(*buckets)[index_of_bucket].index_begin=ind;
			if(index_of_bucket != index_of_bucket_1)
			{
				(*buckets)[index_of_bucket].index_end=ind+1;
				(*buckets)[index_of_bucket_1].index_end=ind+1;
			}
		}else if(ind == (*buckets).size()-1)
		{
			if((*index_pair)[ind].index_of_bucket < (*buckets).size()){
				(*buckets)[(*index_pair)[ind].index_of_bucket].index_end=ind+1;
			}
		}else if (ind+1 < (*index_pair).size())
		{
			long long unsigned int index_of_bucket = (*index_pair)[ind].index_of_bucket;
			long long unsigned int index_of_bucket_1 = (*index_pair)[ind+1].index_of_bucket;

			if(index_of_bucket != index_of_bucket_1)
			{
				(*buckets)[index_of_bucket].index_end=ind+1;
				(*buckets)[index_of_bucket_1].index_begin=ind+1;
			}
		}
	}
}

void build_rgd_final_job(int i, Job* job, std::vector<Bucket>* buckets){
	for(size_t ii = job->index_begin_inclusive; ii < job->index_end_exclusive; ii++){
		long long unsigned int index_begin = (*buckets)[ii].index_begin;
		long long unsigned int index_end = (*buckets)[ii].index_end;
		if(index_begin != -1 && index_end !=-1)
		{
			(*buckets)[ii].number_of_points = index_end - index_begin;
		}
	}
}

void build_rgd(std::vector<Point3D> &points, std::vector<PointBucketIndexPair>& index_pair, std::vector<Bucket>& buckets, GridParameters& rgd_params, int num_threads)
{
	if(num_threads < 1)num_threads = 1;

	index_pair.resize(points.size());
	reindex(points, index_pair, rgd_params, num_threads);

	buckets.resize(rgd_params.number_of_buckets);

	std::vector<Job> jobs = get_jobs(buckets.size(), num_threads);
	std::vector<std::thread> threads;

	for(size_t i = 0; i < jobs.size(); i++){
		threads.push_back(std::thread( build_rgd_init_job, i, &jobs[i], &buckets));
	}

	for(size_t j = 0; j < threads.size(); j++){
		threads[j].join();
	}
	threads.clear();

	jobs = get_jobs(points.size(), num_threads);

	for(size_t i = 0; i < jobs.size(); i++){
		threads.push_back(std::thread( build_rgd_job, i, &jobs[i], &index_pair, &buckets));
	}
	for(size_t j = 0; j < threads.size(); j++){
		threads[j].join();
	}
	threads.clear();


	jobs = get_jobs(buckets.size(), num_threads);

	for(size_t i = 0; i < jobs.size(); i++){
		threads.push_back(std::thread( build_rgd_final_job, i, &jobs[i], &buckets));
	}

	for(size_t j = 0; j < threads.size(); j++){
		threads[j].join();
	}
	threads.clear();
}





