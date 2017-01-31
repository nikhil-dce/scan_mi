#include "mi_plotter.h"

#include <pcl/common/transforms.h>
#include <string>
#include <boost/format.hpp>

void plotErrorGraph() {

    const std::string dataDir = "../../scan_registration_data/";

	int scanNo1 = 180;
	int step = 5;


    std::string gicpFileName = "%strans_gicp/trans_gicp_%d_%d";
    std::string resultFileName = "%strans_scan_mi/trans_result_%d_%d";
    std::string gtFileName = "%sTRANS_GT/trans_%d_%d";
    std::string fileString;

	int scanNo2 = 0;

	pcl::visualization::PCLPlotter *plotterTranslation, *plotterRotation;
	plotterTranslation = new pcl::visualization::PCLPlotter("Ground Truth Translation Plot");
	plotterTranslation->setXTitle("Scan Number");
	plotterTranslation->setYTitle("Mean Squared Translation Error");

	plotterRotation = new pcl::visualization::PCLPlotter("Ground Truth Rotation Plot");
	plotterRotation->setXTitle("Scan Difference");
	plotterRotation->setYTitle("Mean Squared Rotation Error");

	Eigen::Affine3d transformGicp, transformResult, transformGT;

	double txGt(0), tyGt(0), tzGt(0), rollGt(0), pitchGt(0), yawGt(0);
	double tx(0), ty(0), tz(0), roll(0), pitch(0), yaw(0);

	std::vector<std::pair<double, double> > errorTransGicp, errorRotGicp;
	std::vector<std::pair<double, double> > errorTransResult, errorRotResult;
	std::vector<std::pair<double, double> > errorTransBase, errorRotBase;

	int scanCounter = 0;
	while (scanNo1 < 990) {

		scanNo2 = scanNo1 + step;
		fileString = (boost::format(gtFileName)%dataDir%scanNo1%scanNo2).str();
		loadTransform(fileString, transformGT);

		transform_get_rotation_xyz_from_affine(transformGT, &rollGt, &pitchGt, &yawGt);
		transform_get_translation_from_affine(transformGT, &txGt, &tyGt, &tzGt);

		std::cout << "Scan: " << scanNo1 << std::endl;
		std::cout << "Roll: " << rollGt << " Pitch: " << pitchGt << " Yaw: " << yawGt << std::endl;
		std::cout << "x: " << txGt << " y: " << tyGt << " z: " << tzGt << std::endl;

		double errorTrans = square(txGt) + square(tyGt) + square(tzGt);
		errorTrans = sqrt(errorTrans);
		errorTransBase.push_back(std::pair<double, double> (scanNo1, errorTrans));

		double errorRot = square(rollGt) + square(pitchGt) + square(yawGt);
		errorRot = sqrt(errorRot);
		errorRotBase.push_back(std::pair<double, double> (scanNo1, errorRot));

		fileString = (boost::format(resultFileName)%dataDir%scanNo1%scanNo2).str();
		loadTransform(fileString, transformResult);
		transform_get_rotation_xyz_from_affine(transformResult, &roll, &pitch, &yaw);
		transform_get_translation_from_affine(transformResult, &tx, &ty, &tz);

		errorTrans = square(txGt-tx) + square(tyGt-ty) + square(tzGt-tz);
		errorTrans = sqrt(errorTrans);
		errorTransResult.push_back(std::pair<double, double> (scanNo1, errorTrans));

		errorRot = square(rollGt-roll) + square(pitchGt-pitch) + square(yawGt-yaw);
		errorRot = sqrt(errorRot);
		errorRotResult.push_back(std::pair<double, double> (scanNo1, errorRot));

		fileString = (boost::format(gicpFileName)%dataDir%scanNo1%scanNo2).str();
		loadTransform(fileString, transformGicp);
		transform_get_rotation_xyz_from_affine(transformGicp, &roll, &pitch, &yaw);
		transform_get_translation_from_affine(transformGicp, &tx, &ty, &tz);

		errorTrans = square(txGt-tx) + square(tyGt-ty) + square(tzGt-tz);
		errorTrans = sqrt(errorTrans);
		errorTransGicp.push_back(std::pair<double, double> (scanNo1, errorTrans));

		errorRot = square(rollGt-roll) + square(pitchGt-pitch) + square(yawGt-yaw);
		errorRot = sqrt(errorRot);
		errorRotGicp.push_back(std::pair<double, double> (scanNo1, errorRot));

		//
		scanNo1 += 10;
		scanCounter++;
	}

	std::cout << "Total Scans: " << scanCounter << std::endl;

	plotterTranslation->addPlotData(errorTransBase, "Base Error");
	plotterTranslation->addPlotData(errorTransResult, "Result Error");
	plotterTranslation->addPlotData(errorTransGicp, "GICP Error");

	plotterRotation->addPlotData(errorRotBase, "Base Error");
	plotterRotation->addPlotData(errorRotResult, "Result Error");
	plotterRotation->addPlotData(errorRotGicp, "GICP Error");

	plotterRotation->plot();
	plotterTranslation->plot();

}

void
statisticalPlot() {
    
    const std::string dataDir = "../../scan_registration_data/";
    
    int scanNo1 = 180;
    int scanStepS = 5;
    int scanStepM = 10;
    int scanStepL = 15;
    int scanStepXL = 20;
    
    std::string gicpFileName = "%strans_gicp/trans_gicp_%d_%d";
    std::string resultFileName = "%strans_scan_mi_nmsimplex/trans_result_%d_%d";
    std::string gtFileName = "%sTRANS_GT/trans_%d_%d";
    std::string fileString;
    int scanNo2 = 0;
    
    pcl::visualization::PCLPlotter *plotterTranslation, *plotterRotation;
    plotterTranslation = new pcl::visualization::PCLPlotter("Ground Truth Translation Plot");
    plotterTranslation->setXTitle("Scan Difference");
    plotterTranslation->setYTitle("Mean Squared Translation Error");
    
    plotterRotation = new pcl::visualization::PCLPlotter("Ground Truth Rotation Plot");
    plotterRotation->setXTitle("Scan Difference");
    plotterRotation->setYTitle("Mean Squared Rotation Error");
    
    double transGicpErrorS(0), rotGicpErrorS(0), transResultErrorS(0), rotResultErrorS(0), transBaseErrorS(0), rotBaseErrorS(0);
    double transGicpErrorM(0), rotGicpErrorM(0), transResultErrorM(0), rotResultErrorM(0), transBaseErrorM(0), rotBaseErrorM(0);
    double transGicpErrorL(0), rotGicpErrorL(0), transResultErrorL(0), rotResultErrorL(0), transBaseErrorL(0), rotBaseErrorL(0);
    double transGicpErrorXL(0), rotGicpErrorXL(0), transResultErrorXL(0), rotResultErrorXL(0), transBaseErrorXL(0), rotBaseErrorXL(0);
    
    Eigen::Affine3d transformGicp, transformResult, transformGT;
    
    double txGt(0), tyGt(0), tzGt(0), rollGt(0), pitchGt(0), yawGt(0);
    double tx(0), ty(0), tz(0), roll(0), pitch(0), yaw(0);
    
    std::vector<std::pair<double, double> > errorTransGicp, errorRotGicp;
    std::vector<std::pair<double, double> > errorTransResult, errorRotResult;
    std::vector<std::pair<double, double> > errorTransBase, errorRotBase;
    
    int scanCounter = 0;
    while (scanNo1 < 990) {
        
        // step = 5
        scanNo2 = scanNo1 + scanStepS;
        fileString = (boost::format(gtFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformGT);
        
        transform_get_rotation_xyz_from_affine(transformGT, &rollGt, &pitchGt, &yawGt);
        transform_get_translation_from_affine(transformGT, &txGt, &tyGt, &tzGt);
        
        double errorTrans = square(txGt) + square(tyGt) + square(tzGt);
        errorTrans = sqrt(errorTrans);
        transBaseErrorS += errorTrans;
        
        double errorRot = square(rollGt) + square(pitchGt) + square(yawGt);
        errorRot = sqrt(errorRot);
        rotBaseErrorS += errorRot;
        
        fileString = (boost::format(resultFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformResult);
        transform_get_rotation_xyz_from_affine(transformResult, &roll, &pitch, &yaw);
        transform_get_translation_from_affine(transformResult, &tx, &ty, &tz);
        
        errorTrans = square(txGt-tx) + square(tyGt-ty) + square(tzGt-tz);
        errorTrans = sqrt(errorTrans);
        transResultErrorS += errorTrans;
        
        errorRot = square(rollGt-roll) + square(pitchGt-pitch) + square(yawGt-yaw);
        errorRot = sqrt(errorRot);
        rotResultErrorS += errorRot;
        
        fileString = (boost::format(gicpFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformGicp);
        transform_get_rotation_xyz_from_affine(transformGicp, &roll, &pitch, &yaw);
        transform_get_translation_from_affine(transformGicp, &tx, &ty, &tz);
        
        errorTrans = square(txGt-tx) + square(tyGt-ty) + square(tzGt-tz);
        errorTrans = sqrt(errorTrans);
        transGicpErrorS += errorTrans;
        
        errorRot = square(rollGt-roll) + square(pitchGt-pitch) + square(yawGt-yaw);
        errorRot = sqrt(errorRot);
        rotGicpErrorS += errorRot;
        
        // step = 10
        
        scanNo2 = scanNo1 + scanStepM;
        fileString = (boost::format(gtFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformGT);
        
        transform_get_rotation_xyz_from_affine(transformGT, &rollGt, &pitchGt, &yawGt);
        transform_get_translation_from_affine(transformGT, &txGt, &tyGt, &tzGt);
        
        errorTrans = square(txGt) + square(tyGt) + square(tzGt);
        errorTrans = sqrt(errorTrans);
        transBaseErrorM += errorTrans;
        
        errorRot = square(rollGt) + square(pitchGt) + square(yawGt);
        errorRot = sqrt(errorRot);
        rotBaseErrorM += errorRot;
        
        fileString = (boost::format(resultFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformResult);
        transform_get_rotation_xyz_from_affine(transformResult, &roll, &pitch, &yaw);
        transform_get_translation_from_affine(transformResult, &tx, &ty, &tz);
        
        errorTrans = square(txGt-tx) + square(tyGt-ty) + square(tzGt-tz);
        errorTrans = sqrt(errorTrans);
        transResultErrorM += errorTrans;
        
        errorRot = square(rollGt-roll) + square(pitchGt-pitch) + square(yawGt-yaw);
        errorRot = sqrt(errorRot);
        rotResultErrorM += errorRot;
        
        fileString = (boost::format(gicpFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformGicp);
        transform_get_rotation_xyz_from_affine(transformGicp, &roll, &pitch, &yaw);
        transform_get_translation_from_affine(transformGicp, &tx, &ty, &tz);
        
        errorTrans = square(txGt-tx) + square(tyGt-ty) + square(tzGt-tz);
        errorTrans = sqrt(errorTrans);
        transGicpErrorM += errorTrans;
        
        errorRot = square(rollGt-roll) + square(pitchGt-pitch) + square(yawGt-yaw);
        errorRot = sqrt(errorRot);
        rotGicpErrorM += errorRot;
        
        // step = 15
        scanNo2 = scanNo1 + scanStepL;
        fileString = (boost::format(gtFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformGT);
        
        transform_get_rotation_xyz_from_affine(transformGT, &rollGt, &pitchGt, &yawGt);
        transform_get_translation_from_affine(transformGT, &txGt, &tyGt, &tzGt);
        
        errorTrans = square(txGt) + square(tyGt) + square(tzGt);
        errorTrans = sqrt(errorTrans);
        transBaseErrorL += errorTrans;
        
        errorRot = square(rollGt) + square(pitchGt) + square(yawGt);
        errorRot = sqrt(errorRot);
        rotBaseErrorL += errorRot;
        
        fileString = (boost::format(resultFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformResult);
        transform_get_rotation_xyz_from_affine(transformResult, &roll, &pitch, &yaw);
        transform_get_translation_from_affine(transformResult, &tx, &ty, &tz);
        
        errorTrans = square(txGt-tx) + square(tyGt-ty) + square(tzGt-tz);
        errorTrans = sqrt(errorTrans);
        transResultErrorL += errorTrans;
        
        errorRot = square(rollGt-roll) + square(pitchGt-pitch) + square(yawGt-yaw);
        errorRot = sqrt(errorRot);
        rotResultErrorL += errorRot;
        
        fileString = (boost::format(gicpFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformGicp);
        transform_get_rotation_xyz_from_affine(transformGicp, &roll, &pitch, &yaw);
        transform_get_translation_from_affine(transformGicp, &tx, &ty, &tz);
        
        errorTrans = square(txGt-tx) + square(tyGt-ty) + square(tzGt-tz);
        errorTrans = sqrt(errorTrans);
        transGicpErrorL += errorTrans;
        
        errorRot = square(rollGt-roll) + square(pitchGt-pitch) + square(yawGt-yaw);
        errorRot = sqrt(errorRot);
        rotGicpErrorL += errorRot;
        
        // step = 20
        scanNo2 = scanNo1 + scanStepXL;
        fileString = (boost::format(gtFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformGT);
        
        transform_get_rotation_xyz_from_affine(transformGT, &rollGt, &pitchGt, &yawGt);
        transform_get_translation_from_affine(transformGT, &txGt, &tyGt, &tzGt);
        
        errorTrans = square(txGt) + square(tyGt) + square(tzGt);
        errorTrans = sqrt(errorTrans);
        transBaseErrorXL += errorTrans;
        
        errorRot = square(rollGt) + square(pitchGt) + square(yawGt);
        errorRot = sqrt(errorRot);
        rotBaseErrorXL += errorRot;
        
        fileString = (boost::format(resultFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformResult);
        transform_get_rotation_xyz_from_affine(transformResult, &roll, &pitch, &yaw);
        transform_get_translation_from_affine(transformResult, &tx, &ty, &tz);
        
        errorTrans = square(txGt-tx) + square(tyGt-ty) + square(tzGt-tz);
        errorTrans = sqrt(errorTrans);
        transResultErrorXL += errorTrans;
        
        errorRot = square(rollGt-roll) + square(pitchGt-pitch) + square(yawGt-yaw);
        errorRot = sqrt(errorRot);
        rotResultErrorXL += errorRot;
        
        fileString = (boost::format(gicpFileName)%dataDir%scanNo1%scanNo2).str();
        loadTransform(fileString, transformGicp);
        transform_get_rotation_xyz_from_affine(transformGicp, &roll, &pitch, &yaw);
        transform_get_translation_from_affine(transformGicp, &tx, &ty, &tz);
        
        errorTrans = square(txGt-tx) + square(tyGt-ty) + square(tzGt-tz);
        errorTrans = sqrt(errorTrans);
        transGicpErrorXL += errorTrans;
        
        errorRot = square(rollGt-roll) + square(pitchGt-pitch) + square(yawGt-yaw);
        errorRot = sqrt(errorRot);
        rotGicpErrorXL += errorRot;
        
        //
        scanNo1 += 10;
        scanCounter++;
    }
    
    std::cout << "Total Scans: " << scanCounter << std::endl;
    
    transBaseErrorS /= scanCounter;
    rotBaseErrorS /= scanCounter;
    transBaseErrorM /= scanCounter;
    rotBaseErrorM /= scanCounter;
    transBaseErrorL /= scanCounter;
    rotBaseErrorL /= scanCounter;
    transBaseErrorXL /= scanCounter;
    rotBaseErrorXL /= scanCounter;
    
    transGicpErrorS /= scanCounter;
    rotGicpErrorS /= scanCounter;
    transGicpErrorM /= scanCounter;
    rotGicpErrorM /= scanCounter;
    transGicpErrorL /= scanCounter;
    rotGicpErrorL /= scanCounter;
    transGicpErrorXL /= scanCounter;
    rotGicpErrorXL /= scanCounter;
    
    transResultErrorS /= scanCounter;
    rotResultErrorS /= scanCounter;
    transResultErrorM /= scanCounter;
    rotResultErrorM /= scanCounter;
    transResultErrorL /= scanCounter;
    rotResultErrorL /= scanCounter;
    transResultErrorXL /= scanCounter;
    rotResultErrorXL /= scanCounter;
    
    errorTransBase.push_back(std::pair<double, double>(5, transBaseErrorS));
    errorRotBase.push_back(std::pair<double, double>(5, rotBaseErrorS));
    errorTransBase.push_back(std::pair<double, double>(10, transBaseErrorM));
    errorRotBase.push_back(std::pair<double, double>(10, rotBaseErrorM));
    errorTransBase.push_back(std::pair<double, double>(15, transBaseErrorL));
    errorRotBase.push_back(std::pair<double, double>(15, rotBaseErrorL));
    errorTransBase.push_back(std::pair<double, double>(20, transBaseErrorXL));
    errorRotBase.push_back(std::pair<double, double>(20, rotBaseErrorXL));
    
    errorTransResult.push_back(std::pair<double, double>(5, transResultErrorS));
    errorRotResult.push_back(std::pair<double, double>(5, rotResultErrorS));
    errorTransResult.push_back(std::pair<double, double>(10, transResultErrorM));
    errorRotResult.push_back(std::pair<double, double>(10, rotResultErrorM));
    errorTransResult.push_back(std::pair<double, double>(15, transResultErrorL));
    errorRotResult.push_back(std::pair<double, double>(15, rotResultErrorL));
    errorTransResult.push_back(std::pair<double, double>(20, transResultErrorXL));
    errorRotResult.push_back(std::pair<double, double>(20, rotResultErrorXL));
    
    errorTransGicp.push_back(std::pair<double, double>(5, transGicpErrorS));
    errorRotGicp.push_back(std::pair<double, double>(5, rotGicpErrorS));
    errorTransGicp.push_back(std::pair<double, double>(10, transGicpErrorM));
    errorRotGicp.push_back(std::pair<double, double>(10, rotGicpErrorM));
    errorTransGicp.push_back(std::pair<double, double>(15, transGicpErrorL));
    errorRotGicp.push_back(std::pair<double, double>(15, rotGicpErrorL));
    errorTransGicp.push_back(std::pair<double, double>(20, transGicpErrorXL));
    errorRotGicp.push_back(std::pair<double, double>(20, rotGicpErrorXL));
    
    
    plotterTranslation->addPlotData(errorTransBase, "Base Error", vtkChart::POINTS);
    plotterTranslation->addPlotData(errorTransResult, "Result Error", vtkChart::POINTS);
    plotterTranslation->addPlotData(errorTransGicp, "GICP Error", vtkChart::POINTS);
    
    plotterRotation->addPlotData(errorRotBase, "Base Error", vtkChart::POINTS);
    plotterRotation->addPlotData(errorRotResult, "Result Error", vtkChart::POINTS);
    plotterRotation->addPlotData(errorRotGicp, "GICP Error", vtkChart::POINTS);
    
    plotterRotation->plot();
    plotterTranslation->plot();
    
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
        
        double mi = calculateMI(A, transformedScanB, 1);
        
        std::cout << "X: " << xVar << " MI: " << mi << std::endl;
        
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
        
        double mi = calculateMI(A, transformedScanB, 1);
        std::cout << "Y: " << yVar << " MI: " << mi << std::endl;
        
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
        
        double mi = calculateMI(A, transformedScanB, 1);
        std::cout << "Z: " << zVar << " MI: " << mi << std::endl;
        
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
        
        double mi = calculateMI(A, transformedScanB, 1);
        std::cout << "Roll: " << rollVar << " MI: " << mi << std::endl;
        
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
        
        double mi = calculateMI(A, transformedScanB, 1);
        std::cout << "Pitch: " << pitchVar << " MI: " << mi << std::endl;
        
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
        
        double mi = calculateMI(A, transformedScanB, 1);
        std::cout << "Yaw: " << yawVar << " MI: " << mi << std::endl;
        
        data.push_back(std::pair<double, double> (yawDegree, mi));
        yawDegree += plotStepAngle;
        
    }
    
    std::string labelString = (boost::format("Yaw: %f")%actualYawDeg).str();
    plotter->addPlotData(data, labelString.c_str());
    
}
