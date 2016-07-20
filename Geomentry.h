#ifndef GEOMENTRY_H
#define GEOMENTRY_H
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

typedef pcl::PointXYZRGBNormal PointType;
class Geomentry
{
public:
Geomentry();
~Geomentry();

  pcl::PointCloud<PointType> readPointCloud(string pointcloudFileName);
  
  float VectorAngle(Eigen::Vector3f veca, 
		   Eigen::Vector3f vecb);
  
  Eigen::Vector3f CrossProduct(Eigen::Vector3f veca, 
		   Eigen::Vector3f vecb);
  
  float DotProduct(Eigen::Vector3f veca, 
		   Eigen::Vector3f vecb);
  
  double computeCloudResolution(const pcl::PointCloud<PointType>::ConstPtr& cloud);
  
  Eigen::Matrix<float, 3, 3> RotationAboutVector(Eigen::Vector3f rotationAxis, float theta);
  
  pcl::PointCloud< PointType > calculateNormal(pcl::PointCloud< PointType > pointcloud, 
						      float radious);

  cv::Mat findRTfromS(cv::Mat Xcorre, cv::Mat Ycorre, float sinit);
  
  cv::Mat computeS(const cv::Mat &pointsX, const cv::Mat &pointsY, float sinit);
  
  float computescale_s(const cv::Mat &pointsX, const cv::Mat &pointsY, const cv::Mat &R_h);
  
  cv::Mat findCenter(cv::Mat points);
  
  pcl::PointCloud<PointType> computeBarycenterC(const pcl::PointCloud<PointType>::ConstPtr& cloud);
  
  float computeRadius(pcl::PointCloud<PointType> BarycenterC, pcl::PointCloud<PointType>::ConstPtr cloud);
  
};

#endif // GEOMENTRY_H
