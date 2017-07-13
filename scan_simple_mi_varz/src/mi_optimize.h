#ifndef MI_OPTIMIZE_H
#define MI_OPTIMIZE_H

#include "mi_cost.h"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

struct numDerivativeData {
    
    int derivWrt;
    void* data;
    const gsl_vector* pose;
    
};

struct optiData {
    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanA;
    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanB;
};

double
numericalGradient (double p, void* params);

double
f (const gsl_vector*, void*);

void
df (const gsl_vector*, void*, gsl_vector*);

void
fdf (const gsl_vector *pose, void *params,
     double *fValue, gsl_vector *g);

#endif
