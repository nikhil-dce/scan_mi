#ifndef MI_PLOTTER_H
#define MI_PLOTTER_H

#include "mi_util.h"
#include "mi_cost.h"

static float plotRangeTrans = 1.;
static float plotStepTrans = 0.1;

static float plotRangeAngle = 30;
static float plotStepAngle = 1;

void
plotErrorGraph();

void
plotForX (Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
plotForY (Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
plotForZ (Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
plotForRoll (Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
plotForPitch (Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
plotForYaw (Eigen::Affine3d& transform, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B, pcl::visualization::PCLPlotter *plotter);

void
statisticalPlot ();

#endif
