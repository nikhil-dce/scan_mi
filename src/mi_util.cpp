#include "mi_util.h"

#include <iostream>
#include <fstream>
#include <string>

void
loadTransform(std::string filename, Eigen::Affine3d& transform) {
    
//    std::cout << "Loading Transformation: " << filename << std::endl;
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
