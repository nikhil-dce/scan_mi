#include "mi_util.h"
#include "mi_cost.h"

#include <string>
#include <queue>
#include <climits>
#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>
#include <boost/timer/timer.hpp>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>

// static std::string DATA_DIR = "/media/data_raid/nikhil/EXPERIMENT_DATA/";
static std::string DATA_DIR = "/home/nikhil/Desktop/";
static std::string transFilename;
static bool onlyCost = false;
static bool costPlot = false;

int 
registerPointCloudsUsingPatternSearch(mi_transform_t baseT, mi_transform_t result);

int initOptions(int argc, char* argv[]) {

	namespace po = boost::program_options;
	po::options_description desc ("Allowed Options");

	desc.add_options()
    		("help,h", "Usage <Scan 1 Path> <Scan 2 Path> <Transform File>")
    		("costPlot", po::value<bool>(&costPlot), "Plot the cost after naive search")
    		("onlyCost", po::value<bool>(&onlyCost), "Print only the cost");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << "Cuda implementation of Scan MI Registration" << std::endl << desc << std::endl;
		return 1;
	}

	return 0;

}

int main (int argc, char *argv[]) {

	if (initOptions(argc, argv))
		return 0;

	if (argc < 4) {
		std::cerr << "Invalid Arguments " << std::endl;
		return 1;
	}

	int s1, s2;

	s1 = atoi(argv[1]);
	s2 = atoi(argv[2]);
	transFilename = argv[3];

	std::stringstream ss;
	std::vector<MI_ScanPoint> scanA, scanB;
	mi_transform_t trans, transInverse;

	trans = (mi_transform_t) calloc (12, sizeof(double));
	transInverse = (mi_transform_t) calloc (12, sizeof(double));

	 ss << DATA_DIR << "XYZRGBA/" << boost::format("xyzrgba_%04d.txt")%s1;
//	ss << DATA_DIR << "XYZA/" << boost::format("xyza_%04d")%s1;
	loadScan(ss.str(), scanA);

	ss.str(std::string());
	 ss << DATA_DIR << "XYZRGBA/" << boost::format("xyzrgba_%04d.txt")%s2;
//	ss << DATA_DIR << "XYZA/" << boost::format("xyza_%04d")%s2;
	loadScan(ss.str(), scanB);

	ss.str(std::string());
	loadTransform(transFilename, trans);	

	// transform inverse for base trans
	transform_inverse(trans, transInverse);

	boost::timer::cpu_timer timer;
	boost::timer::cpu_times elapsed;

	// allocate memory in gpu
	// load the two scans		
	initializeDeviceData(scanA, scanB);
	
	// preprocess voxel data for scan A
	preprocessScanA();

	elapsed = timer.elapsed();
	std::cout << "Preprocess scan A - CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;

	if (onlyCost) {
		
		double x, y, z, roll, pitch, yaw;
		transform_get_translation_from_affine(trans, &x, &y, &z);
		transform_get_rotation_xyz_from_affine(trans, &roll, &pitch, &yaw);

		std::cout << "X: " << x << std::endl;
		std::cout << "Y: " << y << std::endl;
		std::cout << "Z: " << z << std::endl;
		std::cout << "Roll: " << roll << std::endl;
		std::cout << "Pitch: " << pitch << std::endl;
		std::cout << "Yaw: " << yaw << std::endl;

		// Calculate MI and quit
		setDebug(true);
		double mi = calculateMIForPose (transInverse);			
		std::cout << "Complete search: " << std::endl;
		std::cout << "MI: " << mi << std::endl;
		
		// z = 0;
		// roll = 0;
		// pitch = 0;
		// mi_transform_create (trans, x, y, z, roll, pitch, yaw);	
		// transform_inverse(trans, transInverse);
		// mi = calculateMIForPose (transInverse);		

		// std::cout << "Mini search: " << std::endl;
		// std::cout << "MI: " << mi << std::endl;


		std::cout << "Total time - CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;

		printTransform(transInverse);
		freeDeviceData();

		return 0;

	} else if (costPlot) {

		std::ofstream fout_cost("result_heading_cost_plot");

		mi_transform_t queryTrans = (mi_transform_t) calloc (12, sizeof(double));
		mi_transform_t queryInverse = (mi_transform_t) calloc (12, sizeof(double));
		
		std::cout << "Ground Truth" << std::endl;

		double x, y;
		double tx, ty, tz, rx, ry, rz;
		
		transform_get_rotation_xyz_from_affine(trans, &rx, &ry, &rz);
		transform_get_translation_from_affine(trans, &tx, &ty, &tz);

		std::cout << "GT Roll: " << rx << std::endl;
		std::cout << "GT Pitch: " << ry << std::endl;
		std::cout << "GT Heading: " << rz << std::endl;
		// transform_get_rotation_xyz_from_affine(transInverse, &rx, &ry, &rz);


		// std::cout << "X: " << tx << std::endl;
		// std::cout << "Y: " << ty << std::endl;
		// std::cout << "Z: " << tz << std::endl;
		// std::cout << "Roll: " << rx << std::endl;
		// std::cout << "Pitch: " << ry << std::endl;
		// std::cout << "Yaw: " << rz << std::endl;

		// z = tz;
		// roll = rx;
		// pitch = ry;
		// yaw = rz;

		std::cout << "Inverse \n";
		printTransform(transInverse);
		std::cout << "Trans \n";
		printTransform(trans);

		// for (float heading = rz-0.6; heading <= rz+0.6; heading+=0.02) {

		// 	mi_transform_create(queryTrans, tx, ty, tz, rx, ry, heading);
			
		// 	if (fabs(heading - 0.479107) < 0.02) {
		// 		printTransform(queryTrans);
		// 		double a,b,c;
		// 		transform_get_rotation_xyz_from_affine(queryTrans, &a, &b, &c);
		// 		std::cout << a << ' ' << b << ' ' << c << std::endl;
		// 	}
				
		// 	transform_inverse(queryTrans, queryInverse);
		// 	double mi = calculateMIForPose (queryInverse);
		// 	fout_cost << heading << ' ' << mi << '\n';			
		// }

		for (float x = tx-8; x <= tx+8; x+=0.2) {
			for (float y = ty-8; y <= ty+8; y+=0.2) {

				// take gt as it is

				transform_copy(queryTrans, trans);
				queryTrans[3] = x;
				queryTrans[7] = y;

				// mi_transform_create (queryTrans, x, y, z, roll, pitch, yaw);	
				transform_inverse(queryTrans, queryInverse);
				double mi = calculateMIForPose (queryInverse);			

				// transform_inverse(queryTrans, queryInverse);
				// double a, b, c;
				// transform_get_translation_from_affine(queryInverse, &a, &b, &c);
				//transform_get_rotation_xyz_from_affine(queryInverse, &rx, &ry, &rz);
				fout_cost << x << ' ' << y << ' ' << mi << '\n';

			}
		}

		fout_cost.close();
		freeDeviceData();
		return 0;
	}

	
	// plain old pattern search for 6 dimensions
	mi_transform_t result, resultInverse;
	result = (mi_transform_t) calloc (12, sizeof(double));
	resultInverse = (mi_transform_t) calloc (12, sizeof(double));

	int iterations = registerPointCloudsUsingPatternSearch(transInverse, result);	
	transform_inverse (result, resultInverse);	

	elapsed = timer.elapsed();
	std::cout << "Registration completed with " << iterations << " iterations - CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;

	freeDeviceData();

	ss.str(std::string());
	ss << DATA_DIR << "trans_result_varz_cuda/" << "trans_result_" << s1 << "_" << s2;
	double elapsedTime = elapsed.wall / 1e9;		
	saveTransform (ss.str(), resultInverse, elapsedTime);

	delete[] result;
	delete[] transInverse;
	delete[] trans;
	return 0;
}

double f (const gsl_vector *pose, void* params) {

	// Initialize All Data
	double x, y, z, roll, pitch ,yaw;
	x = gsl_vector_get(pose, 0);
	y = gsl_vector_get(pose, 1);
	z = gsl_vector_get(pose, 2);
	roll = gsl_vector_get(pose, 3);
	pitch = gsl_vector_get(pose, 4);
	yaw = gsl_vector_get(pose, 5);

	mi_transform_t transform;
	transform = (mi_transform_t) calloc(12, sizeof(double));
	mi_transform_create(transform, x, y, z, roll, pitch, yaw);

	double mi  = calculateMIForPose (transform);	
	double cost = -mi;

	return cost;
}

int 
registerPointCloudsUsingPatternSearch(mi_transform_t baseT, mi_transform_t result) {

	double x, y, z, roll, pitch, yaw;

	transform_get_translation_from_affine(baseT, &x, &y, &z);
	transform_get_rotation_xyz_from_affine(baseT, &roll, &pitch, &yaw);

	std::cout << "Base Pose: " << std::endl;

	std::cout << "X: " << x << std::endl;
	std::cout << "Y: " << y << std::endl;
	std::cout << "Z: " << z << std::endl;
	std::cout << "Roll: " << roll << std::endl;
	std::cout << "Pitch: " << pitch << std::endl;
	std::cout << "Yaw: " << yaw << std::endl;

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;

	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *stepSize, *basePose;

	//gsl_multimin_function func;
	gsl_multimin_function func;

	// Base Pose 
	basePose = gsl_vector_alloc (6);
	gsl_vector_set (basePose, 0, x);
	gsl_vector_set (basePose, 1, y);
	gsl_vector_set (basePose, 2, z);
	gsl_vector_set (basePose, 3, roll);
	gsl_vector_set (basePose, 4, pitch);
	gsl_vector_set (basePose, 5, yaw);

    stepSize = gsl_vector_alloc (6);

    gsl_vector_set (stepSize, 0, 8);
    gsl_vector_set (stepSize, 1, 8);
    gsl_vector_set (stepSize, 2, 1);
    gsl_vector_set (stepSize, 3, 0.1);
    gsl_vector_set (stepSize, 4, 0.1);
    gsl_vector_set (stepSize, 5, 0.8);

	// Initialize method and iterate 
	func.n = 6;
	func.f = f;
	func.params = 0;

	s = gsl_multimin_fminimizer_alloc (T, 6);
	gsl_multimin_fminimizer_set (s, &func, basePose, stepSize);

	size_t iter = 0;
	int status;
	double size;

	do  {

		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		// std::cout << gsl_strerror(status) << std::endl;

		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-3);

		if (status == GSL_SUCCESS)
		{
			std::cout << "Converged to minimum\n";			
		}

	}   while (status == GSL_CONTINUE && iter < 500);

	gsl_vector_free (basePose);

	x = gsl_vector_get (s->x, 0);
	y = gsl_vector_get (s->x, 1);
	z = gsl_vector_get (s->x, 2);
	roll = gsl_vector_get (s->x, 3);
	pitch = gsl_vector_get (s->x, 4);
	yaw = gsl_vector_get (s->x, 5);

	
	mi_transform_create (result, x, y, z, roll, pitch, yaw);	
	return iter;
}
