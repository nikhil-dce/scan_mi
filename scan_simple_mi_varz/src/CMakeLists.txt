cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(scan_mi_varz_cuda)

include_directories(/media/nikhil/hdv/gsl-build/build/include)
include_directories(/usr/local/cuda-8.0/samples/common/inc/)
include_directories(/media/nikhil/hdv/scan_mi_util)

link_directories(/usr/local/cuda-8.0/samples/lib /usr/local/cuda-8.0/lib64 /media/nikhil/hdv/gsl-build/build/lib ../build)

find_package(Boost COMPONENTS system program_options timer REQUIRED)
set (CUDA_LIBS cudart cuda)
set (GSL_LIBS gsl gslcblas)


set (LIBS ${LIBS} ${Boost_LIBRARIES})
set (LIBS ${LIBS} ${CUDA_LIBS})
set (LIBS ${LIBS} ${GSL_LIBS})

#file(GLOB util_sources RELATIVE "/home/nikhil/scan_mi_util/" "*.cpp")
add_executable (mi_cuda /media/nikhil/hdv/scan_mi_util/mi_util.cpp main.cpp mi_cost.o)

target_link_libraries (mi_cuda ${LIBS})
