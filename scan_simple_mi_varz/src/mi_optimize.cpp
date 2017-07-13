#include "mi_optimize.h"

#include <boost/unordered_map.hpp>
#include <pcl/common/transforms.h>

double
numericalGradient (double p, void* params) {

	// get the data
	numDerivativeData* numDerivData = (numDerivativeData*) params;

	const gsl_vector* pose = numDerivData->pose;

	// Initialize All Data
	double x, y, z, roll, pitch ,yaw;
	x = gsl_vector_get(pose, 0);
	y = gsl_vector_get(pose, 1);
	z = gsl_vector_get(pose, 2);
	roll = gsl_vector_get(pose, 3);
	pitch = gsl_vector_get(pose, 4);
	yaw = gsl_vector_get(pose, 5);

	if (numDerivData->derivWrt == 1) {
		// wrt x
		x = p;
	} else if (numDerivData->derivWrt == 2) {
		// wrt y
		y = p;
	} else if (numDerivData->derivWrt == 3) {
		// wrt z
		z = p;
	} else if (numDerivData->derivWrt == 4) {
		// wrt roll
		roll = p;
	} else if (numDerivData->derivWrt == 5) {
		// wrt pitch
		pitch = p;
	} else if (numDerivData->derivWrt == 6) {
		// wrt yaw
		yaw = p;
	}

	optiData* data = (optiData*) numDerivData->data;

	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A = data->scanA;
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr baseB = data->scanB;
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());

	// Create Transformation
	Eigen::Affine3d transform = Eigen::Affine3d::Identity();
	transform.translation() << x,y,z;
	transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
	transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
	transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

	// Transform point cloud
	pcl::transformPointCloud(*baseB, *transformedB, transform);

	double mi = calculateMIFromMap(A, transformedB, 1);
	double cost = -mi;

	return cost;
}

void
fdf (const gsl_vector *pose, void *params,
		double *fValue, gsl_vector *g)
{

	//    *fValue = f(pose, params);
	//
	//    df(pose, params, g);

	// Initialize All Data
	double x, y, z, roll, pitch ,yaw;
	x = gsl_vector_get(pose, 0);
	y = gsl_vector_get(pose, 1);
	z = gsl_vector_get(pose, 2);
	roll = gsl_vector_get(pose, 3);
	pitch = gsl_vector_get(pose, 4);
	yaw = gsl_vector_get(pose, 5);

	optiData* data = (optiData*) params;

	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A = data->scanA;
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr baseB = data->scanB;
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());

	// Create Transformation
	Eigen::Affine3d transform = Eigen::Affine3d::Identity();
	transform.translation() << x,y,z;
	transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
	transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
	transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

	// Transform point cloud
	pcl::transformPointCloud(*baseB, *transformedB, transform);

	double mi  = calculateMIFromMap(A, transformedB, 1);
	double cost = -mi;
	*fValue = cost;


//	std::cout << "fdf "  << std::endl;
//	std::cout << "x: " << x << std::endl;
//	std::cout << "y: " << y << std::endl;
//	std::cout << "z: " << z << std::endl;
//	std::cout << "roll: " << roll << std::endl;
//	std::cout << "pitch: " << pitch << std::endl;
//	std::cout << "yaw: " << yaw << std::endl;
//	std::cout << "cost: " << cost << std::endl;

	// Numerical gradient

	double transH = 0.1;
	double angleH = 0.01745; // Appx 1 degree
	double mi_h, cost_h, g_h;

	// x+h
	transform(0,3) = x+transH;
	transform(1,3) = y;
	transform(2,3) = z;
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / transH;
	gsl_vector_set (g, 0, g_h);

	// y+h
	transform(0,3) = x;
	transform(1,3) = y+transH;
	transform(2,3) = z;
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / transH;
	gsl_vector_set (g, 1, g_h);

	// z+h
	transform(0,3) = x;
	transform(1,3) = y;
	transform(2,3) = z+transH;
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / transH;
	gsl_vector_set (g, 2, g_h);

	// roll+h
	transform = Eigen::Affine3d::Identity();
	transform.translation() << x,y,z;
	transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
	transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
	transform.rotate (Eigen::AngleAxisd (roll+angleH, Eigen::Vector3d::UnitX()));
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / angleH;
	gsl_vector_set (g, 3, g_h);

	// pitch+h
	transform = Eigen::Affine3d::Identity();
	transform.translation() << x,y,z;
	transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
	transform.rotate (Eigen::AngleAxisd (pitch+angleH, Eigen::Vector3d::UnitY()));
	transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / angleH;
	gsl_vector_set (g, 4, g_h);

	// yaw+h
	transform = Eigen::Affine3d::Identity();
	transform.translation() << x,y,z;
	transform.rotate(Eigen::AngleAxisd (yaw+angleH, Eigen::Vector3d::UnitZ()));
	transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
	transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / angleH;
	gsl_vector_set (g, 5, g_h);



}

void df (const gsl_vector* pose, void* params, gsl_vector* g) {

	// Initialize All Data
	double x, y, z, roll, pitch ,yaw;
	x = gsl_vector_get(pose, 0);
	y = gsl_vector_get(pose, 1);
	z = gsl_vector_get(pose, 2);
	roll = gsl_vector_get(pose, 3);
	pitch = gsl_vector_get(pose, 4);
	yaw = gsl_vector_get(pose, 5);

//	std::cout << "df "  << std::endl;
//	std::cout << "x: " << x << std::endl;
//	std::cout << "y: " << y << std::endl;
//	std::cout << "z: " << z << std::endl;
//	std::cout << "roll: " << roll << std::endl;
//	std::cout << "pitch: " << pitch << std::endl;
//	std::cout << "yaw: " << yaw << std::endl;


	optiData* data = (optiData*) params;

	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A = data->scanA;
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr baseB = data->scanB;
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());

	// Create Transformation
	Eigen::Affine3d transform = Eigen::Affine3d::Identity();
	transform.translation() << x,y,z;
	transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
	transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
	transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

	// Transform point cloud
	pcl::transformPointCloud(*baseB, *transformedB, transform);

	double mi  = calculateMIFromMap(A, transformedB, 1);
	double cost = -mi;

	// Numerical gradient

	double transH = 0.1;
	double angleH = 0.01745; // Appx 1 degree
	double mi_h, cost_h, g_h;

	// x+h
	transform(0,3) = x+transH;
	transform(1,3) = y;
	transform(2,3) = z;
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / transH;
	gsl_vector_set (g, 0, g_h);

	// y+h
	transform(0,3) = x;
	transform(1,3) = y+transH;
	transform(2,3) = z;

	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / transH;
	gsl_vector_set (g, 1, g_h);

	// z+h
	transform(0,3) = x;
	transform(1,3) = y;
	transform(2,3) = z+transH;
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / transH;
	gsl_vector_set (g, 2, g_h);

	// roll+h
	transform = Eigen::Affine3d::Identity();
	transform.translation() << x,y,z;
	transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
	transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
	transform.rotate (Eigen::AngleAxisd (roll+angleH, Eigen::Vector3d::UnitX()));
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / angleH;
	gsl_vector_set (g, 3, g_h);

	// pitch+h
	transform = Eigen::Affine3d::Identity();
	transform.translation() << x,y,z;
	transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
	transform.rotate (Eigen::AngleAxisd (pitch+angleH, Eigen::Vector3d::UnitY()));
	transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / angleH;
	gsl_vector_set (g, 4, g_h);

	// yaw+h
	transform = Eigen::Affine3d::Identity();
	transform.translation() << x,y,z;
	transform.rotate(Eigen::AngleAxisd (yaw+angleH, Eigen::Vector3d::UnitZ()));
	transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
	transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));
	pcl::transformPointCloud(*baseB, *transformedB, transform);
	mi_h  = calculateMIFromMap(A, transformedB, 1);
	cost_h = -mi_h;
	g_h = (cost_h - cost) / angleH;
	gsl_vector_set (g, 5, g_h);

	//    numDerivativeData* numData = new numDerivativeData();
	//    numData->data =  params;
	//    numData->pose = pose;
	//
	//    gsl_function F;
	//    double result, error;
	//
	//    F.function = &numericalGradient;
	//
	//    numData->derivWrt = 1;
	//    F.params = numData;
	//
	//    gsl_deriv_central (&F, x, transH, &result, &error);
	//    gsl_vector_set (g, 0, result);
	//
	//    numData->derivWrt = 2;
	//    F.params = numData;
	//
	//    gsl_deriv_central (&F, y, transH, &result, &error);
	//    gsl_vector_set (g, 1, result);
	//
	//    numData->derivWrt = 3;
	//    F.params = numData;
	//
	//    gsl_deriv_central (&F, z, transH, &result, &error);
	//    gsl_vector_set (g, 2, result);
	//
	//    numData->derivWrt = 4;
	//    F.params = numData;
	//
	//    gsl_deriv_central (&F, roll, angleH, &result, &error);
	//    gsl_vector_set (g, 3, result);
	//
	//    numData->derivWrt = 5;
	//    F.params = numData;
	//
	//    gsl_deriv_central (&F, pitch, angleH, &result, &error);
	//    gsl_vector_set (g, 4, result);
	//
	//    numData->derivWrt = 6;
	//    F.params = numData;
	//
	//    gsl_deriv_central (&F, yaw, angleH, &result, &error);
	//    gsl_vector_set (g, 5, result);

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

//	std::cout << "f "  << std::endl;
//	std::cout << "x: " << x << std::endl;
//	std::cout << "y: " << y << std::endl;
//	std::cout << "z: " << z << std::endl;
//	std::cout << "roll: " << roll << std::endl;
//	std::cout << "pitch: " << pitch << std::endl;
//	std::cout << "yaw: " << yaw << std::endl;

	optiData* data = (optiData*) params;

	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A = data->scanA;
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr baseB = data->scanB;
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr transformedB = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBA> > (new pcl::PointCloud<pcl::PointXYZRGBA> ());

	// Create Transformation
	Eigen::Affine3d transform = Eigen::Affine3d::Identity();
	transform.translation() << x,y,z;
	transform.rotate(Eigen::AngleAxisd (yaw, Eigen::Vector3d::UnitZ()));
	transform.rotate (Eigen::AngleAxisd (pitch, Eigen::Vector3d::UnitY()));
	transform.rotate (Eigen::AngleAxisd (roll, Eigen::Vector3d::UnitX()));

	// Transform point cloud
	pcl::transformPointCloud(*baseB, *transformedB, transform);

	double mi  = calculateMIFromMap(A, transformedB, 1);
	double cost = -mi;

//	std::cout << "cost: " << cost << std::endl;
	return cost;
}
