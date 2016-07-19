#define PCL_NO_PRECOMPILE
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/visualization/histogram_visualizer.h>
#include <pcl/registration/correspondence_estimation.h>
#include <pcl/kdtree/kdtree_flann.h>
#include<pcl/visualization/pcl_plotter.h>
#include <iostream>
using namespace std;

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <time.h>
#include <math.h>
#include <ANN/ANN.h>
#include <fstream>

#include <cv.h>
#include <highgui.h>

class RHGPE
{
private:
  
public:
  RHGPE();
  ~RHGPE();
};