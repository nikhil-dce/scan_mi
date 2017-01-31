#include "mi_util.h"
#include "mi_optimize.h"
#include "mi_plotter.h"

#include <string>
#include <queue>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>
#include <boost/timer/timer.hpp>

#include <pcl/common/transforms.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/octree/octree_pointcloud.h>
#include <pcl/visualization/pcl_plotter.h>
#include <pcl/io/pcd_io.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>

static std::string DATA_DIR = "../../scan_registration_data/";
static std::string transFilename;
static bool plotJointHistogram = false;
static bool plotX = false;
static bool plotY = false;
static bool plotZ = false;
static bool plotRoll = false;
static bool plotPitch = false;
static bool plotYaw = false;
static float RESOLUTION = 1.;
static bool onlyCost = false;
static bool plotStats = false;
static bool plotError = false;

void
registerPointClouds (Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, Eigen::Affine3d& result);

void 
naivePoseSearch (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, Eigen::Affine3d& transform);

void
registerPointCloudsUsingPatternSearch (Eigen::Affine3d& baseT, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, Eigen::Affine3d& result);

void
branchAndBoundSearch (
		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A,
		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B,
		Eigen::Affine3d& transform);


int initOptions(int argc, char* argv[]) {

	namespace po = boost::program_options;
	po::options_description desc ("Allowed Options");

	desc.add_options()
    		("help,h", "Usage <Scan 1 Path> <Scan 2 Path> <Transform File>")
    		("plotStats", po::value<bool>(&plotStats), "Plot Statistical Data")
    		("plotError", po::value<bool>(&plotError), "Plot Error")
    		("onlyCost", po::value<bool>(&onlyCost), "Print only the cost")

    		("plotX", po::value<bool>(&plotX), "PlotX Score Data")
    		("plotY", po::value<bool>(&plotY), "PlotY Score Data")
    		("plotZ", po::value<bool>(&plotZ), "PlotZ Score Data")
    		("plotRoll", po::value<bool>(&plotRoll), "PlotRoll Score Data")
    		("plotPitch", po::value<bool>(&plotPitch), "PlotPitch Score Data")
    		("plotYaw", po::value<bool>(&plotYaw), "PlotYaw Score Data")

    		("plotRangeTrans", po::value<float>(&plotRangeTrans), "Translation Range to plot")
    		("plotStepTrans", po::value<float>(&plotStepTrans), "Translation steps while plotting")

    		("plotRangeAngle", po::value<float>(&plotRangeAngle), "Rotation Range to plot")
    		("plotStepAngle", po::value<float>(&plotStepAngle), "Rotation steps while plotting")

    		("joint_histogram", po::value<bool>(&plotJointHistogram), "Plot Joint Histogram");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << "Scan MI Registration" << std::endl << desc << std::endl;
		return 1;
	}

	return 0;

}

int main (int argc, char *argv[]) {

	if (initOptions(argc, argv))
		return 0;

	if (plotStats) {
		statisticalPlot();
		return 0;
	} else if (plotError) {
		plotErrorGraph();
		return 0;
	}

	if (argc < 4) {
		std::cerr << "Invalid Arguments " << std::endl;
		return 1;
	}

	int s1, s2;

	s1 = atoi(argv[1]);
	s2 = atoi(argv[2]);
	transFilename = argv[3];

	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanA = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());
	Eigen::Affine3d trans = Eigen::Affine3d::Identity();

	std::stringstream ss;
	ss << DATA_DIR << "scans_pcd/" << boost::format("scan_%04d.pcd")%s1;
	if (pcl::io::loadPCDFile<pcl::PointXYZRGBA> (ss.str(), *scanA)) {
		std::cout << "Error loading cloud file: " << ss.str() << std::endl;
		return (1);
	}

	ss.str(std::string());
	ss << DATA_DIR << "scans_pcd/" << boost::format("scan_%04d.pcd")%s2;
	if (pcl::io::loadPCDFile<pcl::PointXYZRGBA> (ss.str(), *scanB)) {
		std::cout << "Error loading cloud file: " << ss.str() << std::endl;
		return (1);
	}

	ss.str(std::string());
	loadTransform(transFilename, trans);

	boost::timer::cpu_timer timer;


	/*
  // Numerical Gradient Test
  optiData* data = new optiData();
  data->scanA = scanA;
  data->scanB = scanB;

  gsl_vector* testPose = gsl_vector_alloc(6);
  gsl_vector_set_all (testPose, 0);
  gsl_vector_set (testPose, 1, 5.);

  f(testPose, data);
  return 0;

  std :: cout << "dx: " << gsl_vector_get (g, 0) << std::endl;
  std :: cout << "dy: " << gsl_vector_get (g, 1) << std::endl;
  std :: cout << "dz: " << gsl_vector_get (g, 2) << std::endl;
  std :: cout << "droll: " << gsl_vector_get (g, 3) << std::endl;
  std :: cout << "dpitch: " << gsl_vector_get (g, 4) << std::endl; 
  std :: cout << "dyaw: " << gsl_vector_get (g, 5) << std::endl;

  boost::timer::cpu_times elapsed = timer.elapsed();
  std::cout << "Time take for derivative: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;
  onlyCost = true;




  // Naive Search
  Eigen::Affine3d test = trans.inverse();
    branchAndBoundSearch (scanA, scanB, test);
  //naivePoseSearch ( scanA, scanB, test);
  boost::timer::cpu_times elapsed = timer.elapsed();
  std::cout << "Time taken for naive search: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;
  return 0;
	 */

	/*
    // Check mapping function
    pcl::PointCloud<pcl::PointXYZRGBL>::Ptr c = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBL> > (new pcl::PointCloud<pcl::PointXYZRGBL> ());
    Eigen::Affine3d test = trans.inverse();

    double mi = calculateMIFromMap(scanA, scanB, 1);
    std::cout << "MI: " << mi << std::endl;

    boost::timer::cpu_times elapsed = timer.elapsed();
    std::cout << "Time taken for mapping function : " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;

    onlyCost = true;
	 */

	ss.str(std::string());
	// data loaded

	Eigen::Affine3d transformation = trans.inverse();

	if (onlyCost) {

		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());
		pcl::transformPointCloud(*scanB, *transformedB, transformation);
		double mi = calculateMIFromMap(scanA, transformedB, 1);

		std::cout << "MI: " << mi << std::endl;

		boost::timer::cpu_times elapsed = timer.elapsed();
		std::cout << "MI Calculation CPU Time: " << (elapsed.user + elapsed.system) / 1e9 << " seconds" << " Actual Time: " << elapsed.wall / 1e9 << " seconds" << std::endl;

		return 0;
	}

	if (plotX || plotY || plotZ || plotRoll || plotPitch || plotYaw) {

		pcl::visualization::PCLPlotter *plotter;
		plotter = new pcl::visualization::PCLPlotter("Mutual Information");
		plotter->setYTitle("Mutual Information");
		plotter->setTitle(transFilename.c_str());

		if (plotX) {
			plotForX(transformation, scanA, scanB, plotter);
		}

		if (plotY) {
			plotForY(transformation, scanA, scanB, plotter);
		}

		if (plotZ) {
			plotForZ(transformation, scanA, scanB, plotter);
		}

		if (plotRoll) {
			plotForRoll(transformation, scanA, scanB, plotter);
		}

		if (plotPitch) {
			plotForPitch(transformation, scanA, scanB, plotter);
		}

		if (plotYaw) {
			plotForYaw(transformation, scanA, scanB, plotter);
		}

		plotter->plot();

	} else {

		Eigen::Affine3d result;
		registerPointCloudsUsingPatternSearch(transformation, scanA, scanB, result);
		//registerPointClouds(transformation, scanA, scanB, result);

		boost::timer::cpu_times elapsed = timer.elapsed();

		double cpuTime = (elapsed.user + elapsed.system) / 1e9;
		double elapsedTime = elapsed.wall / 1e9;

		std::cout << "Total CPU Time: " <<  cpuTime << " seconds" << " Actual Time: " << elapsedTime << " seconds" << std::endl;

		// save transformation

		ss.str(std::string());
		ss << DATA_DIR << "trans_scan_mi/" << boost::format("trans_result_%d_%d")%s1%s2;

		std::string transResultFile = ss.str();
		std::ofstream fout(transResultFile.c_str());
		fout << result.inverse().matrix() << std::endl << elapsedTime << std::endl;
		fout.close();

	}

	return 0;
}

void
registerPointClouds(Eigen::Affine3d& baseT, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, Eigen::Affine3d& result) {

	// cost -->

	optiData* data = new optiData();
	data->scanA = A;
	data->scanB = B;

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

	//const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_steepest_descent;//gsl_multimin_fdfminimizer_vector_bfgs2;
	gsl_multimin_fdfminimizer *s = NULL;
	gsl_vector *stepSize, *basePose;

	//gsl_multimin_function func;
	gsl_multimin_function_fdf func;

	/* Base Pose */
	basePose = gsl_vector_alloc (6);
	gsl_vector_set (basePose, 0, x);
	gsl_vector_set (basePose, 1, y);
	gsl_vector_set (basePose, 2, z);
	gsl_vector_set (basePose, 3, roll);
	gsl_vector_set (basePose, 4, pitch);
	gsl_vector_set (basePose, 5, yaw);

	/* Set initial step sizes to 1 */
	stepSize = gsl_vector_alloc (6);
	gsl_vector_set (stepSize, 0, 5);
	gsl_vector_set (stepSize, 1, 5);
	gsl_vector_set (stepSize, 2, 5);
	gsl_vector_set (stepSize, 3, 0.5);
	gsl_vector_set (stepSize, 4, 0.5);
	gsl_vector_set (stepSize, 5, 0.5);

	/* Initialize method and iterate */
	func.n = 6;
	func.f = f;
	func.df = df;
	func.fdf = fdf;
	func.params = (void*) data;

	//s = gsl_multimin_fminimizer_alloc (T, 6);
	//gsl_multimin_fminimizer_set (s, &func, basePose, stepSize);

	s = gsl_multimin_fdfminimizer_alloc (T, 6);
	gsl_multimin_fdfminimizer_set (s, &func, basePose, 0.01, 1e-4);

	size_t iter = 0;
	int status;
	double size;

	do  {

		iter++;
		//status = gsl_multimin_fminimizer_iterate(s);

		status = gsl_multimin_fdfminimizer_iterate (s);

		std::cout << gsl_strerror(status) << std::endl;

		if (status)
			break;

		//size = gsl_multimin_fminimizer_size (s);
		//status = gsl_multimin_test_size (size, 1e-3);

		status = gsl_multimin_test_gradient (s->gradient, 1e-3);

		if (status == GSL_SUCCESS)
		{
			std::cout << "Converged to minimum at\n";
		}

		std::cout << "Iteration: " << iter << std::endl;
		std::cout << "Cost: " << s->f << std::endl;

		std::cout << "X: " << gsl_vector_get (s->x, 0) << std::endl;
		std::cout << "Y: " << gsl_vector_get (s->x, 1) << std::endl;
		std::cout << "Z: " << gsl_vector_get (s->x, 2) << std::endl;
		std::cout << "Roll: " << gsl_vector_get (s->x, 3) << std::endl;
		std::cout << "Pitch: " << gsl_vector_get (s->x, 4) << std::endl;
		std::cout << "Yaw: " << gsl_vector_get (s->x, 5) << std::endl;

		//std::cout << "Size: " << size << std::endl;

	}   while (status == GSL_CONTINUE && iter < 500);

	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (basePose);

	x = gsl_vector_get (s->x, 0);
	y = gsl_vector_get (s->x, 1);
	z = gsl_vector_get (s->x, 2);
	roll = gsl_vector_get (s->x, 3);
	pitch = gsl_vector_get (s->x, 4);
	yaw = gsl_vector_get (s->x, 5);

	result = Eigen::Affine3d::Identity();
	result.translation() << x,y,z;
	result.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
	result.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
	result.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

}

void
registerPointCloudsUsingPatternSearch (Eigen::Affine3d& baseT, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, Eigen::Affine3d& result) {

	// cost -->

	optiData* data = new optiData();
	data->scanA = A;
	data->scanB = B;

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

	/* Base Pose */
	basePose = gsl_vector_alloc (6);
	gsl_vector_set (basePose, 0, x);
	gsl_vector_set (basePose, 1, y);
	gsl_vector_set (basePose, 2, z);
	gsl_vector_set (basePose, 3, roll);
	gsl_vector_set (basePose, 4, pitch);
	gsl_vector_set (basePose, 5, yaw);

	/* Set initial step sizes to 1 */
	stepSize = gsl_vector_alloc (6);
	gsl_vector_set (stepSize, 0, 5);
	gsl_vector_set (stepSize, 1, 5);
	gsl_vector_set (stepSize, 2, 2);
	gsl_vector_set (stepSize, 3, 0.5);
	gsl_vector_set (stepSize, 4, 0.5);
	gsl_vector_set (stepSize, 5, 0.5);

	/* Initialize method and iterate */
	func.n = 6;
	func.f = f;
	func.params = (void*) data;

	s = gsl_multimin_fminimizer_alloc (T, 6);
	gsl_multimin_fminimizer_set (s, &func, basePose, stepSize);

	size_t iter = 0;
	int status;
	double size;

	do  {

		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		std::cout << gsl_strerror(status) << std::endl;

		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-3);

		if (status == GSL_SUCCESS)
		{
			std::cout << "Converged to minimum at\n";
		}

		std::cout << "Iteration: " << iter << std::endl;
		std::cout << "Cost: " << s->fval << std::endl;

		std::cout << "X: " << gsl_vector_get (s->x, 0) << std::endl;
		std::cout << "Y: " << gsl_vector_get (s->x, 1) << std::endl;
		std::cout << "Z: " << gsl_vector_get (s->x, 2) << std::endl;
		std::cout << "Roll: " << gsl_vector_get (s->x, 3) << std::endl;
		std::cout << "Pitch: " << gsl_vector_get (s->x, 4) << std::endl;
		std::cout << "Yaw: " << gsl_vector_get (s->x, 5) << std::endl;

	}   while (status == GSL_CONTINUE && iter < 500);

	gsl_vector_free (basePose);

	x = gsl_vector_get (s->x, 0);
	y = gsl_vector_get (s->x, 1);
	z = gsl_vector_get (s->x, 2);
	roll = gsl_vector_get (s->x, 3);
	pitch = gsl_vector_get (s->x, 4);
	yaw = gsl_vector_get (s->x, 5);

	result = Eigen::Affine3d::Identity();
	result.translation() << x,y,z;
	result.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
	result.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
	result.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

}

void naivePoseSearch (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, Eigen::Affine3d& transform) {

	double x,y,z,roll,pitch,yaw;

	transform_get_translation_from_affine(transform, &x, &y, &z);
	transform_get_rotation_xyz_from_affine(transform, &roll, &pitch, &yaw);

	double transRange = 5;
	double angleRange = 1.57;

	double transStep = 1.; // 50 cm
	double angleStep = 0.2; // appx - 5.7 degree

	double headingVar = yaw - angleRange;

	double maxMi = -1;

	double Xo, Yo, Yawo;

	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedScanB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());

	while ( headingVar < (yaw+angleRange) ) {

		double xVar = x - transRange;

		while (xVar < (x+transRange) ) {

			double yVar = y - transRange;

			while (yVar < (y+transRange) ) {

				Eigen::Affine3d t = Eigen::Affine3d::Identity();

				t.translation() << xVar, yVar, z;
				t.rotate(Eigen::AngleAxisd (headingVar, Eigen::Vector3d::UnitZ()));
				t.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
				t.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

				pcl::transformPointCloud(*B, *transformedScanB, t);

				double mi = calculateMI(A, transformedScanB, 1);
				std::cout << "yaw: " << headingVar << " XVar: " << xVar << " yVar: " << yVar << " MI: " << mi << std::endl;

				if (mi>maxMi) {
					Xo = xVar;
					Yo = yVar;
					Yawo = headingVar;
					maxMi = mi;
				}

				yVar += transStep;
			}
			xVar += transStep;
		}
		headingVar += angleStep;
	}

	std::cout << "Optimum MI " << std::endl;
	std::cout << "yawO: " << Yawo << " Xo: " << Xo << " yVar: " << Yo << " MI: " << maxMi << std::endl;


}

struct pose_search {

	double x;
	double y;
	double heading;
	double cost;

};

class ComparePose {
public:
	bool operator() (pose_search& t1, pose_search& t2)
	{

		return t1.cost > t2.cost;

	}
};

#define TRANS_STEP 0.5
#define ANGLE_STEP 0.1

void
branchAndBoundSearch (
		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A,
		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B,
		Eigen::Affine3d& transform) {


	double x,y,z,roll,pitch,yaw;

	transform_get_translation_from_affine(transform, &x, &y, &z);
	transform_get_rotation_xyz_from_affine(transform, &roll, &pitch, &yaw);

	double transRange = 5;
	double angleRange = 0.5;

	float resolution = 8;
	float angleStep = resolution*ANGLE_STEP;
	float transStep = resolution*TRANS_STEP;

	double headingVar = yaw - angleRange;

	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedScanB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());
	Eigen::Affine3d t = Eigen::Affine3d::Identity();
	std::priority_queue<pose_search, std::vector<pose_search>, ComparePose > miQueue;

	std::cout << "RESOLUTION: " << resolution << std::endl;
	while ( headingVar < (yaw+angleRange) ) {

		t = Eigen::Affine3d::Identity();
		t.rotate(Eigen::AngleAxisd (headingVar, Eigen::Vector3d::UnitZ()));
		t.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
		t.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

		//
		double xVar = x - transRange;

		while (xVar < (x+transRange) ) {

			double yVar = y - transRange;

			while (yVar < (y+transRange) ) {

				t.translation() << xVar, yVar, z;
				pcl::transformPointCloud(*B, *transformedScanB, t);

				double mi = calculateMIFromMap(A, transformedScanB, resolution);

				pose_search ps;
				ps.x = xVar;
				ps.y = yVar;
				ps.heading = headingVar;
				ps.cost = -mi;

				miQueue.push(ps);

				std::cout << "yaw: " << headingVar << " XVar: " << xVar << " yVar: " << yVar << " MI: " << mi << std::endl;

				yVar += transStep;
			}
			xVar += transStep;
		}
		headingVar += angleStep;
	}

	std::cout << "Optimum MI " << std::endl;
	pose_search Xo = miQueue.top();

	std::cout << "yaw: " << Xo.heading << " xVar: " << Xo.x << " yVar: " << Xo.y << " MI: " << Xo.cost << std::endl;

	//while (!miQueue.em)

}


