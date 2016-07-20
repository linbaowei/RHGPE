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
  pcl::PointCloud<PointType> barycenterC = geomentry.computeBarycenterC(pointcloud.makeShared());
  float pointcloudradius = geomentry.computeRadius(barycenterC, pointcloud.makeShared());
  pcl::PointCloud<PointType> sphere1 = geomentry.generateSphere(barycenterC, pointcloudradius*2);
  pcl::PointCloud<PointType> sphere2 = geomentry.generateSphere(barycenterC, pointcloudradius*1.2);
  
    
  pcl::visualization::PCLVisualizer mainview("pointcloud");  
  mainview.setPosition(0,0);	
  mainview.setBackgroundColor(0.9, 0.9, 0.9);
  pcl::visualization::PointCloudColorHandlerCustom<PointType> single_color(pointcloud.makeShared(), 108, 166, 205);
  mainview.addPointCloud (pointcloud.makeShared(), single_color, "newpointcloud1");
  mainview.addPointCloudNormals<PointType, pcl::Normal>(pointcloud.makeShared(), pointcloudnormal.makeShared(), 10, 30, "normal1");
  mainview.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "newpointcloud1");
  pcl::visualization::PointCloudColorHandlerCustom<PointType> single_colorBarycenter(barycenterC.makeShared(), 255, 0, 0);
  mainview.addPointCloud (barycenterC.makeShared(), single_colorBarycenter, "single_colorBarycenter");
  mainview.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "single_colorBarycenter");
  pcl::visualization::PointCloudColorHandlerCustom<PointType> single_colorsphere1(sphere1.makeShared(), 0, 255, 0);
  mainview.addPointCloud (sphere1.makeShared(), single_colorsphere1, "single_colorsphere1");
  mainview.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "single_colorsphere1");
  pcl::visualization::PointCloudColorHandlerCustom<PointType> single_colorsphere2(sphere2.makeShared(), 0, 255, 0);
  mainview.addPointCloud (sphere2.makeShared(), single_colorsphere2, "single_colorsphere2");
  mainview.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "single_colorsphere2");
  
  
  
  srand (time (NULL));
  for(int i = 0; i < 10000; ++ i)
  {
    int j1 = sphere1.size() * rand () / (RAND_MAX + 1.0f);
    int j2 = sphere2.size() * rand () / (RAND_MAX + 1.0f);
    int j3 = sphere2.size() * rand () / (RAND_MAX + 1.0f);
    Eigen::Vector3f p1 = sphere1.at(j1).getArray3fMap();
    Eigen::Vector3f p2 = sphere2.at(j2).getArray3fMap();
    Eigen::Vector3f p3 = sphere2.at(j3).getArray3fMap();
    Eigen::Vector3f p1p2 = (p2-p1);
    Eigen::Vector3f p2p3 = (p3-p2);
    Eigen::Vector3f p1p3 = (p3-p1);
    float anglep1p2p3 = geomentry.VectorAngle(-p1p2, p2p3);
    float anglep2p3p1 = geomentry.VectorAngle(-p2p3,-p1p3);
    float anglep3p1p2 = geomentry.VectorAngle(p1p2, p1p3);
    float disp1p2 = p1p2.norm();
    float disp2p3 = p2p3.norm();
    float disp1p3 = p1p3.norm();
    if(disp1p2 > 2*pointcloudradius && disp2p3 > 2*pointcloudradius && disp1p3 > 2*pointcloudradius && anglep1p2p3 < 3.14/2 && anglep2p3p1 < 3.14/2 && anglep3p1p2 < 3.14/2)
    {
//       bool intersaction = false;
//       for(int j = 0; j < pointcloud.size(); ++ j)
//       {
// 	Eigen::Vector3f tmppoint;
// 	tmppoint = pointcloud.at(j).getArray3fMap();
// 	Eigen::Vector3f pointp1 = p1 - tmppoint;
// 	Eigen::Vector3f pointp2 = p2 - tmppoint;
// 	Eigen::Vector3f pointp3 = p3 - tmppoint;
// 	float anglep1pointp2 = geomentry.VectorAngle(pointp1, pointp2);
// 	float anglep1pointp3 = geomentry.VectorAngle(pointp1, pointp3);
// 	float anglep2pointp3 = geomentry.VectorAngle(pointp2, pointp3);
// 	if(abs(anglep1pointp2+anglep1pointp3+anglep2pointp3 - 3.14*2) < 0.001)
// 	{
// 	  intersaction = true;
// 	  continue;
// 	}
//       }
//       if(intersaction)
//       {
	mainview.addLine(sphere1.at(j1), sphere2.at(j2), "line1"+i);
	mainview.addLine(sphere1.at(j1), sphere2.at(j3), "line2"+i);
	mainview.addLine(sphere2.at(j2), sphere2.at(j3), "line3"+i);
	mainview.spinOnce(2000);
	mainview.removeShape("line1"+i);
	mainview.removeShape("line2"+i);
	mainview.removeShape("line3"+i);
	
//       }
    }
    cout << i << endl;
  }
  
  while (!mainview.wasStopped ())
  {
    mainview.spinOnce ();
  }
  
  
  
  return 1;
}