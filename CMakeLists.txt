project(RHGPE)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
SET(CMAKE_BUILD_TYPE Debug)
 
FIND_PACKAGE (Threads)
FIND_PACKAGE(OpenCV REQUIRED)
find_package(PCL 1.7 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(RHGPE Geomentry.cpp RHGPE.cpp main.cpp)
target_link_libraries (RHGPE ${PCL_LIBRARIES} ${OpenCV_LIBS} ANN pthread)

