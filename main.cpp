#include "RHGPE.h"
#include "Geomentry.h"
#include <boost/graph/graph_concepts.hpp>

int main(int argc, char ** argv)
{
  if(argc < 2)
  {
    cerr << "RHGPE cup.ply" << endl;
  }
  string pointcloudstr = argv[1];
  cout << "Loading point cloud File name: " << pointcloudstr << endl;
  Geomentry geomentry;
  pcl::PointCloud<PointType> pointcloud = geomentry.readPointCloud(pointcloudstr);
  pcl::PointCloud<pcl::Normal> pointcloudnormal;
  pointcloudnormal.resize(pointcloud.size()); 
  //pointcloud = geomentry.calculateNormal(pointcloud, 30);
  for(int i = 0; i < pointcloud.size(); ++i)
  {
    pointcloudnormal.at(i).normal_x = pointcloud.at(i).normal_x;
    pointcloudnormal.at(i).normal_y = pointcloud.at(i).normal_y;
    pointcloudnormal.at(i).normal_z = pointcloud.at(i).normal_z;
    //cerr << pointcloudnormal.at(i).normal_x << "\t" <<pointcloudnormal.at(i).normal_y << "\t" <<pointcloudnormal.at(i).normal_z << endl;
  }
  //pcl::io::savePLYFileASCII<PointType>("point1.ply", pointcloud);
  pcl::visualization::PCLVisualizer mainview("pointcloud");  
  mainview.setPosition(0,0);	
  mainview.setBackgroundColor(0.9, 0.9, 0.9);
  pcl::visualization::PointCloudColorHandlerCustom<PointType> single_color(pointcloud.makeShared(), 108, 166, 205);
  mainview.addPointCloud (pointcloud.makeShared(), single_color, "newpointcloud1");
  mainview.addPointCloudNormals<PointType, pcl::Normal>(pointcloud.makeShared(), pointcloudnormal.makeShared(), 10, 30, "normal1");

  while (!mainview.wasStopped ())
  {
    mainview.spinOnce ();
  }
  
  
  
  return 1;
}