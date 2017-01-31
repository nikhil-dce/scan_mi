#ifndef MI_UTIL_H
#define MI_UTIL_H

#include <string>
#include <Eigen/Geometry>

struct MIPoint_Set {
    
    double x, y, z;
    unsigned char r, g, b;
    unsigned char reflectivity;
    
};

void
loadTransform (std::string filename, Eigen::Affine3d& transform);

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

template <typename Type>
inline Type square(Type x)
{
    return x * x; // will only work on types for which there is the * operator.
}


#endif
