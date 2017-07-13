#include "mi_cost.h"
#include <cuda.h>
#include <helper_cuda.h>
#include <iostream>
#include <stdio.h>
#include <cfloat>
#include <boost/timer/timer.hpp>

#include <thrust/extrema.h>
#include <thrust/pair.h>
#include <thrust/device_ptr.h>

MI_ScanPoint *d_scanA, *d_scanB, *d_transformedScanB;
int *d_marginalHistA, *d_marginalHistB, *d_jointHistAB;

float *d_voxelDataA, *d_voxelDataB;
int *d_voxelPoints; float *d_voxelSumZ;

// CUDA compact 
bool *d_voxelPredicateArray;

mi_transform_t d_transform;
int size_a, size_b, voxel_size;

// resolution should be greater than 1
#define MAX_VOXELS_ALONG_DIM 500
#define MAX_VOXELS_ALONG_X 250
#define MAX_VOXELS_ALONG_Y 250
#define MAX_VOXELS_ALONG_Z 100
#define MAX_POINTS_IN_VOXEL 2000 // empirical value
#define MAX_VARIANCE 1
#define VARIANCE_BIN 5e-3
#define VAR_RANGE 200

int DEBUG = false;
const int resolution = 1;

struct xyz_minmax {
	float minX, minY, minZ, maxX, maxY, maxZ;
} minmaxA;

struct transform_to_entropy {
	int totalVoxels;

	transform_to_entropy(int _totalVoxels) : totalVoxels(_totalVoxels) {}

	__device__
	double operator()(const int voxels) {
		double p = (double) voxels / totalVoxels;
		
		if (p == 0)
			return 0;

		p = p * log(p);
		return p;

	}

};

struct expand_to_xyzminmax {
  
  __device__  
  xyz_minmax operator()(const MI_ScanPoint p){
    xyz_minmax result;
    result.minX = p.x;
    result.maxX = p.x;
    result.minY = p.y;
    result.maxY = p.y;
    result.minZ = p.z;
    result.maxZ = p.z;
    return result;
  }

};

struct minmax3_functor {
  
  __device__
  xyz_minmax operator()(const xyz_minmax a, const xyz_minmax b) {
    xyz_minmax result;
    result.minX = (a.minX < b.minX) ? a.minX:b.minX;
    result.maxX = (a.maxX > b.maxX) ? a.maxX:b.maxX;
    result.minY = (a.minY < b.minY) ? a.minY:b.minY;
    result.maxY = (a.maxY > b.maxY) ? a.maxY:b.maxY;
    result.minZ = (a.minZ < b.minZ) ? a.minZ:b.minZ;
    result.maxZ = (a.maxZ > b.maxZ) ? a.maxZ:b.maxZ;
    return result;
  }
};

struct transform_point {

	double r0, r1, r2, r3, r4, r5, r6, r7, r8;
	double t0, t1, t2;	

	transform_point(double _r0, double _r1, double _r2, 
					double _r3, double _r4, double _r5, 
					double _r6, double _r7, double _r8, 
					double _t0, double _t1, double _t2) :
					r0(_r0), r1(_r1), r2(_r2),
					r3(_r3), r4(_r4), r5(_r5),
					r6(_r6), r7(_r7), r8(_r8),
					t0(_t0), t1(_t1), t2(_t2) {}
 
 	__host__ __device__
  MI_ScanPoint operator()(MI_ScanPoint basePoint)
  {
	
	MI_ScanPoint transformed_p;	
	transformed_p.x = r0 * basePoint.x + r1 * basePoint.y + r2 * basePoint.z + t0;
	transformed_p.y = r3 * basePoint.x + r4 * basePoint.y + r5 * basePoint.z + t1;
	transformed_p.z = r6 * basePoint.x + r7 * basePoint.y + r8 * basePoint.z + t2;

	transformed_p.refc = basePoint.refc;

	return transformed_p;
  }
};

struct transform_predicate_to_integer {	
 	
 	__host__ __device__
  int operator()(bool predicate)
  {
	
	if (predicate)
		return 1;
	else
		return 0;
	
  }
};

__host__ __device__ int
get_voxel_index (float px, float py, float pz, int resolution) {

	int key = 0;

	float x = px;
	float y = py;
	float z = pz;

	int keyx = x / resolution;
	int keyy = y / resolution;
	int keyz = z / resolution;

	if (x < 0) 
		keyx -= 1;		

	if (y < 0) 
		keyy -= 1;
	
	if (z < 0) 
		keyz -= 1;		

	keyx += MAX_VOXELS_ALONG_X / 2;
	keyy += MAX_VOXELS_ALONG_Y / 2;
	keyz += MAX_VOXELS_ALONG_Z / 2;
	
	
	key = 	keyx * MAX_VOXELS_ALONG_Y * MAX_VOXELS_ALONG_Z +
			keyy * MAX_VOXELS_ALONG_Z + 
			keyz;			

	return key;
}

__host__ __device__ MI_ScanPoint
get_point_for_key (int key) {

	MI_ScanPoint p;

	int keyx = key / (MAX_VOXELS_ALONG_Y * MAX_VOXELS_ALONG_Z);

	int remain1 = key % (MAX_VOXELS_ALONG_Y * MAX_VOXELS_ALONG_Z);
	int keyy = remain1 / MAX_VOXELS_ALONG_Z;
	int remain2 = remain1 % MAX_VOXELS_ALONG_Z;

	int keyz = remain2;

	keyx -= MAX_VOXELS_ALONG_X / 2;
	keyy -= MAX_VOXELS_ALONG_Y / 2;
	keyz -= MAX_VOXELS_ALONG_Z / 2;

	p.x = keyx;
	p.y = keyy;
	p.z = keyz;
	return p;
}

__global__ void
voxelize_scan (MI_ScanPoint *d_scan, float *d_voxelSumZ, int *d_voxelPoints, int size, int resolution) {

	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < size) {

		MI_ScanPoint scanPoint = d_scan[tid];
		int voxelIndex = get_voxel_index(scanPoint.x, scanPoint.y, scanPoint.z, resolution);
		
		atomicAdd(d_voxelPoints + voxelIndex, 1);		

		float z = scanPoint.z;
		if (z < 0)
			z *= -1;

		atomicAdd(d_voxelSumZ + voxelIndex, z);
				
	}	
}

__global__ void 
compute_voxel_variance (MI_ScanPoint *d_scan, float *d_voxelSumZ, int *d_voxelPoints, float *d_voxelData, int size) {

	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < size) {

		MI_ScanPoint scanPoint = d_scan[tid];
		int voxelIndex = get_voxel_index(scanPoint.x, scanPoint.y, scanPoint.z, resolution);

		int totalPointsInVoxel = d_voxelPoints[voxelIndex];
		float mean = d_voxelSumZ[voxelIndex] / totalPointsInVoxel;			

		float z = scanPoint.z;
		if (z < 0)
			z *= -1;

		float meanDistance = z - mean;

		float xVar = (meanDistance * meanDistance) / totalPointsInVoxel;
		atomicAdd (d_voxelData + voxelIndex, xVar);		
	}

}

__global__ void
transform_voxelize_scan (MI_ScanPoint* d_scan, float* d_voxelSumZ, int *d_voxelPoints, MI_ScanPoint* d_transformed_scan, 
						int size, int resolution, double r0, 
						double r1, double r2, double r3, 
						double r4, double r5, double r6, 
						double r7, double r8, double t0, double t1, double t2) {

	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < size) {

		MI_ScanPoint basePoint = d_scan[tid];

		MI_ScanPoint scanPoint;	
		scanPoint.x = r0 * basePoint.x + r1 * basePoint.y + r2 * basePoint.z + t0;
		scanPoint.y = r3 * basePoint.x + r4 * basePoint.y + r5 * basePoint.z + t1;
		scanPoint.z = r6 * basePoint.x + r7 * basePoint.y + r8 * basePoint.z + t2;
		scanPoint.refc = basePoint.refc;

		d_transformed_scan[tid] = scanPoint;
		int voxelIndex = get_voxel_index(scanPoint.x, scanPoint.y, scanPoint.z, resolution);
		
		atomicAdd(d_voxelPoints + voxelIndex, 1);		

		float z = scanPoint.z;
		if (z < 0)
			z *= -1;

		atomicAdd(d_voxelSumZ + voxelIndex, z);
	}	
}

__global__ void
create_histogram (float* d_voxelDataA, float* d_voxelDataB, bool* predicate, 
				  int* d_marginalHistA, int* d_marginalHistB, int* d_jointHistAB, 
				  int size) {

	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < size) {
		
		if (predicate[tid]) { 
			
			// calculate histogram for only the overlappint region

			// tid - voxelIdentifier
			// 2 global reads
			float dataA = d_voxelDataA[tid];
			float dataB = d_voxelDataB[tid];

			int indexA = (dataA + VARIANCE_BIN / 2) / VARIANCE_BIN;
			int indexB = (dataB + VARIANCE_BIN / 2) / VARIANCE_BIN;

			int jointIndex = indexA * VAR_RANGE + indexB;
			// 2 global ATOMIC adds
			// Use shared memory here
			atomicAdd(d_marginalHistA+indexA, 1);
			atomicAdd(d_marginalHistB+indexB, 1);		
			atomicAdd(d_jointHistAB+jointIndex, 1);		

		}
	}	
}

__global__ void
setupPredicateForOverlappingRegion(bool *d_voxelPredicateArray, int size, int resolution,
								   float minX, float minY, float minZ, float maxX, float maxY, float maxZ) {

	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < size) {

		// check if tid in overlapping region
		MI_ScanPoint p = get_point_for_key(tid);

		if (p.x < minX || p.x >= maxX ||
			p.y < minY || p.y >= maxY ||
			p.z < minZ || p.z >= maxZ ) {

			// predicate false
			d_voxelPredicateArray[tid] = false;
			
		} else {
			// predicate true
			d_voxelPredicateArray[tid] = true;		
		}

	}

}

void
initializeDeviceData (std::vector<MI_ScanPoint> h_scanA, std::vector<MI_ScanPoint> h_scanB) {	
	
	size_a = h_scanA.size();
	size_b = h_scanB.size();
	
	voxel_size = MAX_VOXELS_ALONG_X*MAX_VOXELS_ALONG_Y*MAX_VOXELS_ALONG_Z;

	std::cout << "Initializing Data" << std::endl;
	std::cout << "SizeA: " << size_a << std::endl;
	std::cout << "SizeB: " << size_b << std::endl;

	boost::timer::cpu_timer timer;
	boost::timer::cpu_times elapsed;
	
	checkCudaErrors(cudaMalloc(&d_scanA, sizeof(struct MI_ScanPoint) * size_a));
	checkCudaErrors(cudaMalloc(&d_scanB, sizeof(struct MI_ScanPoint) * size_b));
	checkCudaErrors(cudaMalloc(&d_transformedScanB, sizeof(MI_ScanPoint) * size_b)); 

	if (DEBUG) {
		elapsed = timer.elapsed();
		std::cout << "CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;	
		std::cout << "Allocating voxel space " << std::endl;
	}

	// needed to temorarily store reflectivity values
	checkCudaErrors(cudaMalloc(&d_voxelPoints, sizeof(int) * voxel_size));	// to store total points	
	checkCudaErrors(cudaMalloc(&d_voxelSumZ, sizeof(float) * voxel_size));	//store sum
	checkCudaErrors(cudaMalloc(&d_voxelDataA, sizeof(float) * voxel_size));	//store variances
	checkCudaErrors(cudaMalloc(&d_voxelDataB, sizeof(float) * voxel_size));	

	checkCudaErrors(cudaMalloc(&d_voxelPredicateArray, sizeof(bool) * voxel_size));		


	checkCudaErrors(cudaMalloc(&d_marginalHistA, sizeof(int) * VAR_RANGE)); 
	checkCudaErrors(cudaMalloc(&d_marginalHistB, sizeof(int) * VAR_RANGE)); 
	checkCudaErrors(cudaMalloc(&d_jointHistAB, sizeof(int) * VAR_RANGE * VAR_RANGE));	

	elapsed = timer.elapsed();
	std::cout << "Device Allocation completed in time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;		
	
	// loading scans in device	
	checkCudaErrors(cudaMemcpy(d_scanA, &(h_scanA[0]), sizeof(MI_ScanPoint) * size_a, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_scanB, &(h_scanB[0]), sizeof(MI_ScanPoint) * size_b, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaGetLastError());
	// checkCudaErrors(cudaDeviceSynchronize());
	
	elapsed = timer.elapsed();
	std::cout << "Data initialized at time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;	
}

void
preprocessScanA () {
	
	const int blockSize = 128;
	int gridSize = (size_a + blockSize) / blockSize;	

	// initialize voxel data to zero
	checkCudaErrors(cudaMemset(d_voxelDataA, 0, sizeof(float) * voxel_size));

	checkCudaErrors(cudaMemset(d_voxelSumZ, 0, sizeof(float) * voxel_size));
	
	// initialize temp voxel data to zero
	checkCudaErrors(cudaMemset(d_voxelPoints, 0, sizeof(int) * voxel_size));	

	// voxelizing scan A
	// This stores the count of the number of A points in voxel in d_voxelDataA 
	voxelize_scan<<<gridSize, blockSize>>> (d_scanA, d_voxelSumZ, d_voxelPoints, size_a, resolution);
	compute_voxel_variance<<<gridSize, blockSize>>> (d_scanA, d_voxelSumZ, d_voxelPoints, d_voxelDataA, size_a);
	
	gridSize = (voxel_size + blockSize) / blockSize;
	
	// transform_reduce to compute minmax3 for scan A
	thrust::device_ptr<MI_ScanPoint> d_thrustScanA (d_scanA);	

	xyz_minmax limit_init;
	limit_init.minX = FLT_MAX;
	limit_init.maxX = FLT_MIN;
	limit_init.minY = FLT_MAX;
	limit_init.maxY = FLT_MIN;
	limit_init.minZ = FLT_MAX;
	limit_init.maxZ = FLT_MIN;	

	minmaxA = thrust::transform_reduce(d_thrustScanA, d_thrustScanA + size_a, expand_to_xyzminmax(), limit_init, minmax3_functor());
}

double
calculateMIForPose (mi_transform_t t) {
	
	if (DEBUG) 
		printTransform(t);

	const int blockSize = 128;
	int gridSize;

	boost::timer::cpu_timer timer;
	boost::timer::cpu_times elapsed;

	// Initialize scan B data
	checkCudaErrors(cudaMemset(d_marginalHistA, 0, sizeof(int) * VAR_RANGE));
	checkCudaErrors(cudaMemset(d_marginalHistB, 0, sizeof(int) * VAR_RANGE));
	checkCudaErrors(cudaMemset(d_jointHistAB, 0, sizeof(int) * VAR_RANGE * VAR_RANGE));

	// initialize temp voxel data to zero
	checkCudaErrors(cudaMemset(d_voxelDataB, 0, sizeof(float) * voxel_size));	
	checkCudaErrors(cudaMemset(d_voxelPoints, 0, sizeof(int) * voxel_size));
	checkCudaErrors(cudaMemset(d_voxelSumZ, 0, sizeof(float) * voxel_size));

	// no need for this
	// checkCudaErrors(cudaMemset(d_voxelPredicateArray, false, sizeof(int) * MAX_POINTS_IN_VOXEL * MAX_POINTS_IN_VOXEL)); 
	if (DEBUG) {
		elapsed = timer.elapsed();	
		std::cout << "CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;
		
		std::cout << "KernelLaunch - Get voxel hash for transformed scan B" << std::endl;
	}

	gridSize = (size_b + blockSize) / blockSize;
	transform_voxelize_scan<<<gridSize, blockSize>>> (d_scanB, d_voxelSumZ, d_voxelPoints, d_transformedScanB, size_b, resolution,	
																													t[0], t[1], t[2], 
																													t[4], t[5], t[6], 
																													t[8], t[9], t[10], 
																													t[3], t[7], t[11]);
	compute_voxel_variance<<<gridSize, blockSize>>> (d_transformedScanB, d_voxelSumZ, d_voxelPoints, d_voxelDataB, size_b);

	// gridSize = (voxel_size + blockSize) / blockSize;
	// //set_voxel_data<<<gridSize, blockSize>>> (d_tempVoxelData, d_voxelDataB, voxel_size);

	if (DEBUG) {
		elapsed = timer.elapsed();
		std::cout << "CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;

		std::cout << "KernelLaunch - Transform reduce for minmax3 for transformed Scan B " << std::endl;
	}

	// thrust find minmax3
	// init xyz_minmax
	thrust::device_ptr<MI_ScanPoint> d_thrustTransformedScanB (d_transformedScanB);	
	thrust::device_ptr<int> d_thrustHistA (d_marginalHistA);
	thrust::device_ptr<int> d_thrustHistB (d_marginalHistB);
	thrust::device_ptr<int> d_thrustHistAB (d_jointHistAB);
	thrust::device_ptr<bool> d_thrustPredicateVoxel (d_voxelPredicateArray);

	xyz_minmax limit_init;
	limit_init.minX = FLT_MAX;
	limit_init.maxX = FLT_MIN;
	limit_init.minY = FLT_MAX;
	limit_init.maxY = FLT_MIN;
	limit_init.minZ = FLT_MAX;
	limit_init.maxZ = FLT_MIN;	
	
	xyz_minmax minmaxB = thrust::transform_reduce(d_thrustTransformedScanB, d_thrustTransformedScanB + size_b, expand_to_xyzminmax(), limit_init, minmax3_functor());

	if (DEBUG) {
		elapsed = timer.elapsed();
		std::cout << "CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;	
		std::cout << "CPU calculation overlapping region from minmax3 A and B" << std::endl;
	}

	float minX, minY, minZ, maxX, maxY, maxZ;
	minX = std::max(minmaxA.minX, minmaxB.minX);
	maxX = std::min(minmaxA.maxX, minmaxB.maxX);
	minY = std::max(minmaxA.minY, minmaxB.minY);
	maxY = std::min(minmaxA.maxY, minmaxB.maxY);
	minZ = std::max(minmaxA.minZ, minmaxB.minZ);
	maxZ = std::min(minmaxA.maxZ, minmaxB.maxZ);

	if (DEBUG) {

		std::cout << "minmaxA" << std::endl;
		std::cout << "MinX" << minmaxA.minX << " MaxX: " << minmaxA.maxX << std::endl;
		std::cout << "MinY" << minmaxA.minY << " MaxY: " << minmaxA.maxY << std::endl;
		std::cout << "MinZ" << minmaxA.minZ << " MaxZ: " << minmaxA.maxZ << std::endl;		
		std::cout << "minmaxB" << std::endl;
		std::cout << "MinX" << minmaxB.minX << " MaxX: " << minmaxB.maxX << std::endl;
		std::cout << "MinY" << minmaxB.minY << " MaxY: " << minmaxB.maxY << std::endl;
		std::cout << "MinZ" << minmaxB.minZ << " MaxZ: " << minmaxB.maxZ << std::endl;		
		std::cout << "Overlap" << std::endl;
		std::cout << "MinX" << minX << " MaxX: " << maxX << std::endl;
		std::cout << "MinY" << minY << " MaxY: " << maxY << std::endl;
		std::cout << "MinZ" << minZ << " MaxZ: " << maxZ << std::endl;
	}

	minX = floor(minX);
	minY = floor(minY);
	minZ = floor(minZ);

	maxX = ceil(maxX);
	maxY = ceil(maxY);
	maxZ = ceil(maxZ);	

	if (DEBUG) {
		elapsed = timer.elapsed();
		std::cout << "CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;

		std::cout << "Kernel Launch - SetupPredicate Array for voxels"	<< std::endl;	
	}

	// Not using compact algorithm for cuda for now
	// as the histogram calculation is not an expensive operation
	// Will launch the kernel for all voxels and only consdier the ones with true predicate for now
	gridSize = (voxel_size + blockSize) / blockSize;
	setupPredicateForOverlappingRegion<<<gridSize, blockSize>>> (d_voxelPredicateArray, voxel_size, 1, minX, minY, minZ, maxX, maxY, maxZ);
	
	if (DEBUG) {
		elapsed = timer.elapsed();
		std::cout << "CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;
		std::cout << "Creating histograms from predicate and voxelMapping"	<< std::endl;	
	}	
	
	//std::cout << "Calulate overlapping region size using predicate transform_reduce" << std::endl;
	//int initSum = 0;
	//int numberOfVoxels = thrust::transform_reduce(d_thrustPredicateVoxel, d_thrustPredicateVoxel+voxel_size , transform_predicate_to_integer(), initSum, thrust::plus<int>());
	int numberOfVoxels = (maxX - minX) * (maxY - minY) * (maxZ - minZ);
	//elapsed = timer.elapsed();
	//std::cout << "CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;	
	
	if (DEBUG) {

		// float *h_voxelA = (float*) malloc (sizeof(float) * voxel_size);		
		// float *h_voxelB = (float*) malloc (sizeof(float) * voxel_size);
		// bool *h_predicate = (bool*) malloc (sizeof(bool) * voxel_size);

		// checkCudaErrors(cudaMemcpy(h_voxelA, d_voxelDataA, sizeof(float) * voxel_size, cudaMemcpyDeviceToHost));		
		// checkCudaErrors(cudaMemcpy(h_voxelB, d_voxelDataB, sizeof(float) * voxel_size, cudaMemcpyDeviceToHost));		
		// checkCudaErrors(cudaMemcpy(h_predicate, d_voxelPredicateArray, sizeof(bool) * voxel_size, cudaMemcpyDeviceToHost));

		// saveFloatArray ("DEBUG/voxel_A.txt", h_voxelA, voxel_size);		
		// saveFloatArray ("DEBUG/voxel_B.txt", h_voxelB, voxel_size);	
		// savePredicate ("DEBUG/predicate.txt", h_predicate, voxel_size);
	}
	
	gridSize = (voxel_size + blockSize) / blockSize;
	create_histogram<<<gridSize, blockSize>>> (d_voxelDataA, d_voxelDataB, d_voxelPredicateArray, d_marginalHistA, d_marginalHistB, d_jointHistAB, voxel_size);

	if (DEBUG) {
		elapsed = timer.elapsed();
		std::cout << "CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;	
		std::cout << "Kernel launch - Transform to entropy and reduce " << std::endl;
	}

	double init_marginal = 0.;
	// Convert this into single kernel
	// Entropy of A on the other hand will change as the overlapping region changes
	double marginalEntropyA = thrust::transform_reduce(d_thrustHistA, d_thrustHistA+VAR_RANGE , transform_to_entropy(numberOfVoxels), init_marginal, thrust::plus<double>());
	double marginalEntropyB = thrust::transform_reduce(d_thrustHistB, d_thrustHistB+VAR_RANGE , transform_to_entropy(numberOfVoxels), init_marginal, thrust::plus<double>());
	double jointEntropyAB   = thrust::transform_reduce(d_thrustHistAB, d_thrustHistAB+(VAR_RANGE*VAR_RANGE) , transform_to_entropy(numberOfVoxels), init_marginal, thrust::plus<double>());

	checkCudaErrors(cudaGetLastError());	

	double mi = -(marginalEntropyA + marginalEntropyB - jointEntropyAB);

	if (DEBUG) {
		elapsed = timer.elapsed();
		std::cout << "CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;	

		std::cout << "Total voxels: " << numberOfVoxels << std::endl;
		std::cout << "Entropy A: " << marginalEntropyA << std::endl;
		std::cout << "Entropy B: " << marginalEntropyB << std::endl;
		std::cout << "Joint Entropy AB: " << jointEntropyAB << std::endl;
		std::cout << "MI: " << mi << std::endl;
	}

	
	return mi;	
}

void
freeDeviceData() {
	cudaFree (d_scanA);
	cudaFree (d_transformedScanB);
	cudaFree (d_scanB);
	cudaFree (d_voxelDataA);
	cudaFree (d_voxelDataB);
	cudaFree (d_marginalHistA);
	cudaFree (d_marginalHistB);
	cudaFree (d_jointHistAB);
	cudaFree (d_voxelPredicateArray);
	cudaFree (d_voxelSumZ);
	cudaFree (d_voxelPoints);
}

void
setDebug(bool _debug) {
	DEBUG = _debug;
}
