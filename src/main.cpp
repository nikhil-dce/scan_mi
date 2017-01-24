#include <string>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>
#include <boost/unordered_map.hpp>

#include <pcl/common/transforms.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/octree/octree_pointcloud.h>
#include <pcl/visualization/pcl_plotter.h>

static std::string DATA_DIR = "../../scan_registration_data/";
static float RESOLUTION = 1.;
static std::string transFilename;
static bool plotJointHistogram = false;
static bool plotX = false;
static bool plotY = false;
static bool plotZ = false;
static bool plotRoll = false;
static bool plotPitch = false;
static bool plotYaw = false;

static float plotRangeTrans = 0.;
static float plotStepTrans = 0.;

static float plotRangeAngle = 0;
static float plotStepAngle = 0;

void inline
transform_get_translation_from_affine(Eigen::Affine3d& t, double *x, double *y, double *z);
void inline
transform_get_rotation_xyz_from_affine(Eigen::Affine3d& t, double *x, double *y, double *z);

void
loadTransform(std::string filename, Eigen::Affine3d& transform);

double
calculateMI(pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanA, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanB);

void
plotForX(Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
plotForY(Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
plotForZ(Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
plotForRoll(Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
plotForPitch(Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
plotForYaw(Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);


int initOptions(int argc, char* argv[]) {

	namespace po = boost::program_options;
	po::options_description desc ("Allowed Options");

	desc.add_options()
							("help,h", "Usage <Scan 1 Path> <Scan 2 Path> <Transform File>")

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

int
main (int argc, char *argv[]) {

	if (argc < 4) {
		std::cerr << "Invalid Arguments " << std::endl;
		return 1;
	}

	if (initOptions(argc, argv))
		return 0;

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

	ss.str(std::string());
	// data loaded

	Eigen::Affine3d transformation = trans.inverse();

	// transformation

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

	return 0;
}

void
plotForX(Eigen::Affine3d& t, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter) {

	double x,y,z,roll,pitch,yaw;

	transform_get_rotation_xyz_from_affine(t, &roll, &pitch, &yaw);
	transform_get_translation_from_affine(t, &x, &y, &z);

	double xVar = x - plotRangeTrans;
	std::vector<std::pair<double, double> > data;

	while (xVar < (x + plotRangeTrans) ) {

		Eigen::Affine3d transform = Eigen::Affine3d::Identity();
		transform.translation() << xVar,y,z;
		transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
		transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
		transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedScanB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());
		pcl::transformPointCloud(*B, *transformedScanB, transform);

		double mi = calculateMI(A, transformedScanB);

		data.push_back(std::pair<double, double> (xVar, mi));

		xVar += plotStepTrans;
	}

	std::string labelString = (boost::format("X: %f")%x).str();
	plotter->addPlotData(data, labelString.c_str());

}

void
plotForY(Eigen::Affine3d& t, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter) {

	double x,y,z,roll,pitch,yaw;

	transform_get_rotation_xyz_from_affine(t, &roll, &pitch, &yaw);
	transform_get_translation_from_affine(t, &x, &y, &z);

	double yVar = y - plotRangeTrans;
	std::vector<std::pair<double, double> > data;

	while (yVar < (y + plotRangeTrans) ) {

		Eigen::Affine3d transform = Eigen::Affine3d::Identity();
		transform.translation() << x, yVar, z;
		transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
		transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
		transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedScanB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());
		pcl::transformPointCloud(*B, *transformedScanB, transform);

		double mi = calculateMI(A, transformedScanB);

		data.push_back(std::pair<double, double> (yVar, mi));

		yVar += plotStepTrans;
	}

	std::string labelString = (boost::format("Y: %f")%y).str();
	plotter->addPlotData(data, labelString.c_str());

}

void
plotForZ(Eigen::Affine3d& t, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter) {

	double x,y,z,roll,pitch,yaw;

	transform_get_rotation_xyz_from_affine(t, &roll, &pitch, &yaw);
	transform_get_translation_from_affine(t, &x, &y, &z);

	double zVar = z - plotRangeTrans;
	std::vector<std::pair<double, double> > data;

	while (zVar < (z + plotRangeTrans) ) {

		Eigen::Affine3d transform = Eigen::Affine3d::Identity();
		transform.translation() << x, y, zVar;
		transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
		transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
		transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedScanB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());
		pcl::transformPointCloud(*B, *transformedScanB, transform);

		double mi = calculateMI(A, transformedScanB);

		data.push_back(std::pair<double, double> (zVar, mi));

		zVar += plotStepTrans;
	}

	std::string labelString = (boost::format("Z: %f")%z).str();
	plotter->addPlotData(data, labelString.c_str());

}

void
plotForRoll(Eigen::Affine3d& t, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter) {

	double x,y,z,roll,pitch,yaw;

	transform_get_rotation_xyz_from_affine(t, &roll, &pitch, &yaw);
	transform_get_translation_from_affine(t, &x, &y, &z);

	std::vector<std::pair<double, double> > data;

	double actualRollDeg = roll * 180 / M_PI;
	double rollDegree = actualRollDeg - plotRangeAngle;

	while (rollDegree <= (actualRollDeg + plotRangeAngle)) {

		double rollVar = rollDegree * M_PI / 180;

		Eigen::Affine3d transform = Eigen::Affine3d::Identity();
		transform.translation() << x, y, z;
		transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
		transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
		transform.rotate (Eigen::AngleAxisd (rollVar, Eigen::Vector3d::UnitX()));

		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedScanB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());
		pcl::transformPointCloud(*B, *transformedScanB, transform);

		double mi = calculateMI(A, transformedScanB);

		data.push_back(std::pair<double, double> (rollDegree, mi));
		rollDegree += plotStepAngle;

	}

	std::string labelString = (boost::format("Roll: %f")%actualRollDeg).str();
	plotter->addPlotData(data, labelString.c_str());

}

void plotForPitch (Eigen::Affine3d& t, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter) {

	double x,y,z,roll,pitch,yaw;

	transform_get_rotation_xyz_from_affine(t, &roll, &pitch, &yaw);
	transform_get_translation_from_affine(t, &x, &y, &z);

	std::vector<std::pair<double, double> > data;

	double actualPitchDeg = pitch * 180 / M_PI;
	double pitchDegree = actualPitchDeg - plotRangeAngle;

	while (pitchDegree <= (actualPitchDeg + plotRangeAngle)) {

		double pitchVar = pitchDegree * M_PI / 180;

		Eigen::Affine3d transform = Eigen::Affine3d::Identity();
		transform.translation() << x, y, z;
		transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
		transform.rotate (Eigen::AngleAxisd (pitchVar, Eigen::Vector3d::UnitY()));
		transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedScanB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());
		pcl::transformPointCloud(*B, *transformedScanB, transform);

		double mi = calculateMI(A, transformedScanB);

		data.push_back(std::pair<double, double> (pitchDegree, mi));
		pitchDegree += plotStepAngle;

	}

	std::string labelString = (boost::format("Pitch: %f")%actualPitchDeg).str();
	plotter->addPlotData(data, labelString.c_str());

}

void plotForYaw (Eigen::Affine3d& t, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter) {

	double x,y,z,roll,pitch,yaw;

	transform_get_rotation_xyz_from_affine(t, &roll, &pitch, &yaw);
	transform_get_translation_from_affine(t, &x, &y, &z);

	std::vector<std::pair<double, double> > data;

	double actualYawDeg = yaw * 180 / M_PI;
	double yawDegree = actualYawDeg - plotRangeAngle;

	while (yawDegree <= (actualYawDeg + plotRangeAngle)) {

		double yawVar = yawDegree * M_PI / 180;

		Eigen::Affine3d transform = Eigen::Affine3d::Identity();
		transform.translation() << x, y, z;
		transform.rotate(Eigen::AngleAxisd (yawVar, Eigen::Vector3d::UnitZ()));
		transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
		transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedScanB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());
		pcl::transformPointCloud(*B, *transformedScanB, transform);

		double mi = calculateMI(A, transformedScanB);

		data.push_back(std::pair<double, double> (yawDegree, mi));
		yawDegree += plotStepAngle;

	}

	std::string labelString = (boost::format("Yaw: %f")%actualYawDeg).str();
	plotter->addPlotData(data, labelString.c_str());

}

void getCombinedScan (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A,
		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B,
		pcl::PointCloud<pcl::PointXYZRGBL>::Ptr combinedScan,
		pcl::octree::OctreePointCloud<pcl::PointXYZRGBL>::Ptr octree,
		int& minXo, int& minYo,
		int& minZo, int& maxXo,
		int& maxYo, int& maxZo) {

	octree->setResolution(RESOLUTION);
	double minX = std::numeric_limits<float>::max (), minY = std::numeric_limits<float>::max (), minZ = std::numeric_limits<float>::max ();
	double maxX = -std::numeric_limits<float>::max(), maxY = -std::numeric_limits<float>::max(), maxZ = -std::numeric_limits<float>::max();

	for (int i = 0; i < A->size(); i++) {
		pcl::PointXYZRGBA temp (A->points[i]);

		if (!pcl::isFinite (temp)) //Check to make sure transform didn't make point not finite
			continue;

		if (temp.x < minX)
			minX = temp.x;
		if (temp.y < minY)
			minY = temp.y;
		if (temp.z < minZ)
			minZ = temp.z;
		if (temp.x > maxX)
			maxX = temp.x;
		if (temp.y > maxY)
			maxY = temp.y;
		if (temp.z > maxZ)
			maxZ = temp.z;

		pcl::PointXYZRGBL p;
		p.x = temp.x;
		p.y = temp.y;
		p.z = temp.z;
		p.r = temp.r;
		p.g = temp.g;
		p.b = temp.b;
		p.label = 0;

		combinedScan->push_back(p);
	}

	double minX_A = minX, minY_A = minY, minZ_A = minZ;
	double maxX_A = maxX, maxY_A = maxY, maxZ_A = maxZ;

	minX = std::numeric_limits<float>::max (), minY = std::numeric_limits<float>::max (), minZ = std::numeric_limits<float>::max ();
	maxX = -std::numeric_limits<float>::max(), maxY = -std::numeric_limits<float>::max(), maxZ = -std::numeric_limits<float>::max();

	for (int i = 0; i < B->size(); i++) {

		pcl::PointXYZRGBA temp (B->points[i]);

		if (!pcl::isFinite (temp)) //Check to make sure transform didn't make point not finite
			continue;

		if (temp.x < minX)
			minX = temp.x;
		if (temp.y < minY)
			minY = temp.y;
		if (temp.z < minZ)
			minZ = temp.z;
		if (temp.x > maxX)
			maxX = temp.x;
		if (temp.y > maxY)
			maxY = temp.y;
		if (temp.z > maxZ)
			maxZ = temp.z;

		pcl::PointXYZRGBL p;
		p.x = temp.x;
		p.y = temp.y;
		p.z = temp.z;
		p.r = temp.r;
		p.g = temp.g;
		p.b = temp.b;
		p.label = 1;

		combinedScan->push_back(p);
	}

	double minX_B = minX, minY_B = minY, minZ_B = minZ;
	double maxX_B = maxX, maxY_B = maxY, maxZ_B = maxZ;

	minX = std::min(minX_A, minX_B);
	maxX = std::max(maxX_A, maxX_B);
	minY = std::min(minY_A, minY_B);
	maxY = std::max(maxY_A, maxY_B);
	minZ = std::min(minZ_A, minZ_B);
	maxZ = std::max(maxZ_A, maxZ_B);

	octree->defineBoundingBox (floor(minX), floor(minY), floor(minZ), ceil(maxX), ceil(maxY), ceil(maxZ) );
	octree->setInputCloud(combinedScan);
	octree->addPointsFromInputCloud();
	octree->getBoundingBox(minX, minY, minZ, maxX, maxY, maxZ);

//	cout << "MinX: " << minX << std::endl;
//	cout << "MinY: " << minY << std::endl;
//	cout << "MinZ: " << minZ << std::endl;
//	cout << "MaxX: " << maxX << std::endl;
//	cout << "MaxY: " << maxY << std::endl;
//	cout << "MaxZ: " << maxZ << std::endl;

	minXo = floor(std::max(minX_A, minX_B));
	minYo = floor(std::max(minY_A, minY_B));
	minZo = floor(std::max(minZ_A, minZ_B));

	maxXo = ceil(std::min(maxX_A, maxX_B));
	maxYo = ceil(std::min(maxY_A, maxY_B));
	maxZo = ceil(std::min(maxZ_A, maxZ_B));

//	cout << "MinXo: " << minXo << std::endl;
//	cout << "MinYo: " << minYo << std::endl;
//	cout << "MinZo: " << minZo << std::endl;
//	cout << "MaxXo: " << maxXo << std::endl;
//	cout << "MaxYo: " << maxYo << std::endl;
//	cout << "MaxZo: " << maxZo << std::endl;

}

double calculateMI (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanA,
		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanB) {

	long maxPointsInVoxelForEachScan = 10000;
	int minXo, minYo, minZo;
	int maxXo, maxYo, maxZo;

	std::vector<std::pair<double, double> > data;

	pcl::PointCloud<pcl::PointXYZRGBL>::Ptr combinedScan = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBL> > (new pcl::PointCloud<pcl::PointXYZRGBL> ());
	pcl::octree::OctreePointCloud<pcl::PointXYZRGBL>::Ptr octree = boost::shared_ptr <pcl::octree::OctreePointCloud<pcl::PointXYZRGBL> > (new pcl::octree::OctreePointCloud<pcl::PointXYZRGBL> (RESOLUTION));

	getCombinedScan(scanA, scanB, combinedScan, octree, minXo, minYo, minZo, maxXo, maxYo, maxZo);

	typename pcl::octree::OctreePointCloud<pcl::PointXYZRGBL>::LeafNodeIterator iterator = octree->leaf_begin();
	std::vector<int> indexVector;

	boost::unordered::unordered_map<int, int> umapA;
	boost::unordered::unordered_map<int, int> umapB;
	boost::unordered::unordered_map<long, int> umapAB;

	int overlappingLeavesInOctree = 0;
	while (iterator != octree->leaf_end()) {

		typename pcl::octree::OctreePointCloud<pcl::PointXYZRGBA>::LeafNode *leafNode = (typename pcl::octree::OctreePointCloud<pcl::PointXYZRGBA>::LeafNode*) *iterator;

		indexVector.clear();
		leafNode->getContainer().getPointIndices(indexVector);

		int pointsFromA(0), pointsFromB(0);

		pcl::PointXYZRGBL p;
		for (int i = 0; i < indexVector.size(); i++) {

			int scanIndex = indexVector[i];
			p = combinedScan->at(scanIndex);

			if (p.label == 0)
				pointsFromA++;
			else
				pointsFromB++;
		}

		int xleafMin, xLeafMax, yLeafMin, yLeafMax, zLeafMin, zLeafMax;
		xleafMin = p.x;
		yLeafMin = p.y;
		zLeafMin = p.z;

		// TODO this only work for RESOLUTION = 1
		// change logic
		if (p.x < 0)
			xleafMin--;
		if (p.y < 0)
			yLeafMin--;
		if (p.z < 0)
			zLeafMin--;

		xLeafMax = xleafMin + RESOLUTION;
		yLeafMax = yLeafMin + RESOLUTION;
		zLeafMax = zLeafMin + RESOLUTION;

		// check leaf in overlapping region

		if (	xleafMin > minXo && xLeafMax < maxXo &&
				yLeafMin > minYo && yLeafMax < maxYo &&
				zLeafMin > minZo && zLeafMax < maxZo) {

			// leaf in overlapping region
			overlappingLeavesInOctree++;

			if (umapA.find(pointsFromA) == umapA.end())
				umapA.insert(std::pair<int, int> (pointsFromA, 1));
			else
				umapA[pointsFromA]++;

			if (umapB.find(pointsFromB) == umapB.end())
				umapB.insert(std::pair<int, int> (pointsFromB, 1));
			else
				umapB[pointsFromB]++;

			long key = pointsFromA * maxPointsInVoxelForEachScan + pointsFromB;

			if (umapAB.find(key) == umapAB.end()) {
				umapAB.insert(std::pair<long, int> (key, 1));
			} else {
				umapAB[key]++;
			}

			data.push_back(std::pair<int, int> (pointsFromA, pointsFromB));

		}

		iterator++;
	}


	long total = (maxXo - minXo) * (maxYo - minYo) * (maxZo - minZo);
//	cout << "Total voxels including empty space in overlapping cuboid: " << total << std::endl;

	long leavesNotInOctree = total-overlappingLeavesInOctree;

	if (umapA.find(0) == umapA.end()) {
		umapA.insert(std::pair<int, int> (0, leavesNotInOctree));
	} else {
		umapA[0] += leavesNotInOctree;
	}

	if (umapB.find(0) == umapB.end()) {
		umapB.insert(std::pair<int, int> (0, leavesNotInOctree));
	} else {
		umapB[0] += leavesNotInOctree;
	}

	umapAB.insert(std::pair<long, int> (0, leavesNotInOctree));

	double jointEntropy(0), marginalEntropyA(0), marginalEntropyB(0);

	boost::unordered::unordered_map<int, int>::iterator itr = umapA.begin();
	while (itr != umapA.end()) {

		double probability = ( (double) itr->second) / total;
		marginalEntropyA += -probability * log(probability);

		itr++;
	}

	itr = umapB.begin();
	while (itr != umapB.end()) {

		double probability = ( (double) itr->second) / total;
		marginalEntropyB += -probability * log(probability);

		itr++;
	}

	boost::unordered::unordered_map<long, int>::iterator itrAB = umapAB.begin();
	while (itrAB != umapAB.end()) {

		double probability = ( (double) itrAB->second) / total;
		jointEntropy += -probability * log(probability);

		itrAB++;
	}

//	cout << "Marginal Entropy A: " << marginalEntropyA << std::endl;
//	cout << "Marginal Entropy B: " << marginalEntropyB << std::endl;
//	cout << "Joint Entropy: " << jointEntropy << std::endl;

	double mi = marginalEntropyA + marginalEntropyB - jointEntropy;

	if (plotJointHistogram) {
		pcl::visualization::PCLPlotter *plotter;
		plotter = new pcl::visualization::PCLPlotter("Joint Histogram Plot");
		plotter->setXTitle("A");
		plotter->setYTitle("B");

		plotter->setTitle(transFilename.c_str());
		plotter->addPlotData(data, "Joint Histogram", vtkChart::POINTS);
		plotter->plot();
	}

	return mi;
}

void
loadTransform(std::string filename, Eigen::Affine3d& transform) {

	std::cout << "Loading Transformation: " << filename << std::endl;
	transform = Eigen::Affine3d::Identity();
	std::ifstream in(filename.c_str());
	if (!in) {
		std::stringstream err;
		err << "Error loading transformation " << filename.c_str() << std::endl;
		std::cerr << err.str();
	}

	std::string line;
	for (int i = 0; i < 4; i++) {
		std::getline(in,line);

		std::istringstream sin(line);
		for (int j = 0; j < 4; j++) {
			sin >> transform (i,j);
		}
	}

	in.close();

}


void inline
transform_get_translation_from_affine(Eigen::Affine3d& t, double *x, double *y, double *z) {
	*x = t(0,3);
	*y = t(1,3);
	*z = t(2,3);
}

void inline
transform_get_rotation_xyz_from_affine(Eigen::Affine3d& t, double *x, double *y, double *z) {
	double a = t(0,0); // cycz
	double b = t(0,1); // -cysz
	double c = t(0,2); // sy
	double d = t(1,2); // -sxcy
	double e = t(2,2); // cxcy

	*y = asin(c);
	*z = atan2(-b, a);
	*x = atan2(-d, e);
}
