#ifndef MI_UTIL_H
#define MI_UTIL_H

#include <string>
#include <math.h>
#include <vector>

struct MI_ScanPoint {

float x, y, z;
//unsigned char r, g, b;
unsigned char refc;

};

typedef double* mi_transform_t; 
// Array of 12 

void
mi_transform_create (mi_transform_t result, double x, double y, double z, double roll, double pitch, double yaw);

void
mi_transform_left_multiply (mi_transform_t left, mi_transform_t right);

void
mi_transform_identity (mi_transform_t t);

void 
mi_transform_translate (mi_transform_t t, double &x, double &y, double &z);

void
mi_transform_rotate_x (mi_transform_t t, double &theta);

void
mi_transform_rotate_y (mi_transform_t t, double &theta);

void
mi_transform_rotate_z (mi_transform_t t, double &theta);

void
printTransform (mi_transform_t t);

void
loadTransform (std::string filename, mi_transform_t t);

void
saveTransform (std::string filename, mi_transform_t t, float time);

void
loadScan (std::string filename, std::vector<MI_ScanPoint>& scan);

void inline
transform_copy (mi_transform_t dest, mi_transform_t src) {

  for(int r = 0; r < 12; r++)
      dest[r] = src[r];  
}

void inline
transform_inverse (mi_transform_t in, mi_transform_t out) {

    double temp, t1, t2, t3;

    transform_copy(out, in);

    temp = out[1];
    out[1] = out[4];
    out[4] = temp;

    temp = out[2];
    out[2] = out[8];
    out[8] = temp;

    temp = out[6];
    out[6] = out[9];
    out[9] = temp;

    t1 = 
    -out[0] * out[3]
    -out[1] * out[7]
    -out[2] * out[11];
    t2 = 
    -out[4] * out[3]
    -out[5] * out[7]
    -out[6] * out[11];
    t3 = 
    -out[8] * out[3]
    -out[9] * out[7]
    -out[10] * out[11];

    out[3] = t1;
    out[7] = t2;
    out[11] = t3;

}

void
saveHistogram (std::string filename, int *hist, int size);

void
savePredicate (std::string filename, bool *hist, int size);

void
saveFloatArray (std::string filename, float *data, int size);

void inline
transform_get_translation_from_affine(mi_transform_t t, double *x, double *y, double *z) {
    *x = t[3];
    *y = t[7];
    *z = t[11];
}

void inline
transform_get_rotation_xyz_from_affine(mi_transform_t t, double *x, double *y, double *z) {
    double a = t[0]; // cycz
    double b = t[1]; // -cysz
    double c = t[2]; // sy
    double d = t[6]; // -sxcy
    double e = t[10]; // cxcy
    
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
