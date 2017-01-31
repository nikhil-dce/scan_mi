#ifndef MI_COST_H
#define MI_COST_H

#include "mi_util.h"

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/octree/octree_pointcloud.h>
#include <pcl/visualization/pcl_plotter.h>

#define CUDA_VERSION 0;

double
calculateMI (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanA, pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanB, float resolution);

void
getCombinedScan (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A,
                      pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B,
                      pcl::PointCloud<pcl::PointXYZRGBL>::Ptr combinedScan,
                      pcl::octree::OctreePointCloud<pcl::PointXYZRGBL>::Ptr octree,
                      int& minXo, int& minYo,
                      int& minZo, int& maxXo,
                      int& maxYo, int& maxZo,
                      float resolution);

void
createMap (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A,
           pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B,
           pcl::PointCloud<pcl::PointXYZRGBL>::Ptr combinedScan,
           float resolution);

double
calculateMIFromMap (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanA,
                    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanB,
                    float resolution);

static long P1 = 5923;
static long P2 = 9109;

long
getKeyForPoint (float x, float y, float z, float resolution);

#endif
