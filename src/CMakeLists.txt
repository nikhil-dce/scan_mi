cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(scan_mi)

find_package(PCL 1.7 REQUIRED COMPONENTS common io visualization octree features)
find_package(Boost COMPONENTS program_options timer REQUIRED)
find_package(GSL REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (mi mi_optimize.cpp mi_util.cpp mi_cost.cpp mi_plotter.cpp main.cpp)
 
target_link_libraries (mi ${PCL_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES})
