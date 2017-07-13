#include "mi_util.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

void
loadTransform(std::string filename, mi_transform_t t) {
    
    //std::cout << "Loading Transformation: " << filename << std::endl;
    
    std::ifstream in(filename.c_str());
    if (!in) {
        std::stringstream err;
        err << "Error loading transformation " << filename.c_str() << std::endl;
        std::cerr << err.str();
    }
    
    std::string line;
    for (int i = 0; i < 3; i++) {
        std::getline(in,line);
        std::istringstream sin(line);
        for (int j = 0; j < 4; j++) {
        	double a(0);
            sin >> a;
            t[i*4+j] = a;
        }
    }

    in.close();
}

void
saveTransform (std::string filename, mi_transform_t t, float time) {

    std::ofstream fout(filename.c_str());
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            fout << t[i*4+j] << ' ';
        }
        fout << '\n';
    }
    
    fout << "0 0 0 1 \n" << time << '\n';
    
    fout.close();

}

void
loadScan (std::string filename, std::vector<MI_ScanPoint>& scan) {

    std::ifstream in(filename.c_str());
    if (!in) {
        std::stringstream err;
        err << "Error loading scan " << filename.c_str() << std::endl;
        std::cerr << err.str();
    }

    std::string line;    
    unsigned char r, g, b, ref;
    while (std::getline(in, line)) {
        std::istringstream sin(line);
        MI_ScanPoint p;


        sin >> p.x >> p.y >> p.z >> r >> g >> b >> p.refc;
        scan.push_back(p);
    }
    
    in.close();

}

void
printTransform (mi_transform_t t) {

	std::cout << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3] << std::endl;
	std::cout << t[4] << ' ' << t[5] << ' ' << t[6] << ' ' << t[7] << std::endl;
	std::cout << t[8] << ' ' << t[9] << ' ' << t[10] << ' ' << t[11] << std::endl;
	std::cout << 0 << ' ' << 0 << ' ' << 0 << ' ' << 1 << std::endl;

}

void
mi_transform_create (mi_transform_t t, double x, double y, double z, double roll, double pitch, double yaw) {
    
    // translation
    t[3] = x;
    t[7] = y; 
    t[11]= z;

    // rotation
    double cx = cos(roll), sx = sin(roll), cy = cos(pitch), sy = sin(pitch), cz = cos(yaw), sz = sin(yaw);
    t[0] = cy*cz;
    t[1] = -cy * sz;
    t[2] = sy;
    t[4] = cx*sz + sx*sy*cz;
    t[5] = cx*cz - sx*sy*sz;
    t[6] = -sx*cy;
    t[8] = sx*sz - cx*sy*cz;
    t[9] = cx*sy*sz + sx*cz;
    t[10] = cx*cy;
    
}

void
mi_transform_left_multiply (mi_transform_t left, mi_transform_t right) {
    
    mi_transform_t result = (mi_transform_t) calloc (12, sizeof(double));

    result[0] = right[0] * left[0] + right[1] * left[4] + right[2] * left[8];
    result[1] = right[0] * left[1] + right[1] * left[5] + right[2] * left[9];
    result[2] = right[0] * left[2] + right[1] * left[6] + right[2] * left[10];
    result[3] = right[0] * left[3] + right[1] * left[7] + right[2] * left[11] + right[3];

    result[4] = right[4] * left[0] + right[5] * left[4] + right[6] * left[8];
    result[5] = right[4] * left[1] + right[5] * left[5] + right[6] * left[9];
    result[6] = right[4] * left[2] + right[5] * left[6] + right[6] * left[10];
    result[7] = right[4] * left[3] + right[5] * left[7] + right[6] * left[11] + right[7];

    result[8] = right[8] * left[0] + right[9] * left[4] + right[10] * left[8];  
    result[9] = right[8] * left[1] + right[9] * left[5] + right[10] * left[9];
    result[10] = right[8] * left[2] + right[9] * left[6] + right[10] * left[10];
    result[11] = right[8] * left[3] + right[9] * left[7] + right[10] * left[11] + right[11];
    
    transform_copy (left, result);    
}

void
mi_transform_identity (mi_transform_t t) {
    t[0] = 1;
    t[1] = 0;
    t[2] = 0;
    t[3] = 0;
    t[4] = 0;
    t[5] = 1;
    t[6] = 0;
    t[7] = 0;
    t[8] = 0;
    t[9] = 0;
    t[10]= 1;
    t[11]= 0;
}

void 
mi_transform_translate (mi_transform_t t, double &x, double &y, double &z) {

    t[3] += x;
    t[7] += y;
    t[11]+= z;

}

void
mi_transform_rotate_x (mi_transform_t t, double &theta) {
    
    double ctheta = cos(theta), stheta = sin(theta);
    
    mi_transform_t temp = (mi_transform_t) calloc(12, sizeof(double));
    mi_transform_identity(temp);    
    temp[5] = ctheta;
    temp[6] = -stheta;
    temp[9] = stheta;
    temp[10] = ctheta;

    mi_transform_left_multiply(t,temp);
}

void
mi_transform_rotate_y (mi_transform_t t, double &theta) {
    
    double ctheta = cos(theta), stheta = sin(theta);
    
    mi_transform_t temp = (mi_transform_t) calloc(12, sizeof(double));
    mi_transform_identity(temp);    
    temp[0] = ctheta;
    temp[2] = stheta;
    temp[8] = -stheta;
    temp[10] = ctheta;

    mi_transform_left_multiply(t,temp);
}

void
mi_transform_rotate_z (mi_transform_t t, double &theta) {
    
    double ctheta = cos(theta), stheta = sin(theta);
    
    mi_transform_t temp = (mi_transform_t) calloc(12, sizeof(double));
    mi_transform_identity(temp);    
    temp[0] = ctheta;
    temp[1] = -stheta;
    temp[4] = stheta;
    temp[5] = ctheta;

    mi_transform_left_multiply(t,temp);
}

void
saveHistogram (std::string filename, int *hist, int size) {

    std::ofstream fout(filename.c_str());

    for (int i=0; i<size; i++) {
        fout << hist[i] << std::endl;
    }
    
    fout.close();
}

void
savePredicate (std::string filename, bool *hist, int size) {

    std::ofstream fout(filename.c_str());

    for (int i=0; i<size; i++) {
        fout << hist[i] << std::endl;
    }
    
    fout.close();
}

void
saveFloatArray (std::string filename, float *data, int size) {

    std::ofstream fout(filename.c_str());

    for (int i=0; i<size; i++) {
        fout << data[i] << std::endl;
    }
    
    fout.close();
       
}

