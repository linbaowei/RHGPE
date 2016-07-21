#include "RHGPE.h"
#include "Geomentry.h"
#include <boost/graph/graph_concepts.hpp>

pcl::PointCloud<PointType> holdingPoints(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3, pcl::PointCloud<PointType> inner)
{
  Eigen::Vector3f p1p2 = (p2-p1);
  Eigen::Vector3f p1p3 = (p3-p1);
  float minang1 = 10e10;
  float minang2 = 10e10;
  Geomentry geo;
  pcl::PointCloud<PointType> holdpoints;
  holdpoints.resize(2);
  for(int i = 0; i < inner.size(); ++ i)
  {
    Eigen::Vector3f tmppoint = inner.at(i).getArray3fMap();
    Eigen::Vector3f p1tmppoint = tmppoint - p1;
    float angle1 = geo.VectorAngle(p1p2, p1tmppoint);
    if(angle1 < minang1)
    {
      minang1 = angle1;
      holdpoints.at(0) = inner.at(i);
    }
    float angle2 = geo.VectorAngle(p1p3, p1tmppoint);
    if(angle2 < minang2)
    {
      minang2 = angle2;
      holdpoints.at(1) = inner.at(i);
    }
  }
  return holdpoints;
}

pcl::PointCloud<PointType> ajustHoldingPoints(pcl::PointCloud<PointType> pointcloud, pcl::PointCloud<PointType> holdingPoints)
{
  pcl::PointCloud<PointType> ajustHoldingPoints;
  ajustHoldingPoints.resize(2);
  Geomentry geo;
  Eigen::Vector3f holdpointleft = holdingPoints.at(0).getArray3fMap();  
  Eigen::Vector3f holdpointright = holdingPoints.at(1).getArray3fMap();      
  int stepNUM = 10;
  
  pcl::KdTreeFLANN<PointType> kdtree;
  kdtree.setInputCloud (pointcloud.makeShared());
  PointType searchPoint1, searchPoint2;
  int K = 1;
  std::vector<int> pointIdxNKNSearch(K);
  std::vector<float> pointNKNSquaredDistance(K);
  
  for(int i = 0; i< stepNUM; ++ i)
  {
    Eigen::Vector3f holdpointleftnormal = holdingPoints.at(0).getNormalVector3fMap();
    Eigen::Vector3f holdpointrightnormal = holdingPoints.at(1).getNormalVector3fMap();
    Eigen::Vector3f handline = holdpointright - holdpointleft;
    float stepangle1 = geo.VectorAngle(holdpointleftnormal, -handline);
    float stepangle2 = geo.VectorAngle(holdpointrightnormal, handline);
    Eigen::Vector3f rotationVec1 = geo.CrossProduct(holdpointleftnormal, -handline);
    Eigen::Vector3f rotationVec2 = geo.CrossProduct(holdpointrightnormal, handline);
    Eigen::Matrix<float, 3, 3> rotationMat1 = geo.RotationAboutVector(rotationVec1, stepangle1/stepNUM);
    Eigen::Matrix<float, 3, 3> rotationMat2 = geo.RotationAboutVector(rotationVec2, stepangle2/stepNUM);
    holdpointleft = rotationMat1 * holdpointleft;
    holdpointright = rotationMat2 * holdpointright;
    searchPoint1.x = holdpointleft[0];
    searchPoint1.y = holdpointleft[1];
    searchPoint1.z = holdpointleft[2];
    if ( kdtree.nearestKSearch (searchPoint1, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0 )
    {
      for (size_t i = 0; i < pointIdxNKNSearch.size (); ++i)
      {
	holdpointleft[0] = pointcloud.points[ pointIdxNKNSearch[i] ].x;
	holdpointleft[1] = pointcloud.points[ pointIdxNKNSearch[i] ].y;
	holdpointleft[2] = pointcloud.points[ pointIdxNKNSearch[i] ].z; 
	ajustHoldingPoints.at(0) = pointcloud.points[ pointIdxNKNSearch[i] ];
      }
    }
    
    searchPoint2.x = holdpointright[0];
    searchPoint2.y = holdpointright[1];
    searchPoint2.z = holdpointright[2];
    if ( kdtree.nearestKSearch (searchPoint2, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0 )
    {
      for (size_t i = 0; i < pointIdxNKNSearch.size (); ++i)
      {
	holdpointright[0] = pointcloud.points[ pointIdxNKNSearch[i] ].x;
	holdpointright[1] = pointcloud.points[ pointIdxNKNSearch[i] ].y;
	holdpointright[2] = pointcloud.points[ pointIdxNKNSearch[i] ].z; 
	ajustHoldingPoints.at(1) = pointcloud.points[ pointIdxNKNSearch[i] ];
      }
    }    
  }
  
  return ajustHoldingPoints;
}

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
    pcl::PointCloud<PointType> cord, inner;
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
      int count1 = 0;
      for(int j = 0; j < pointcloud.size(); ++ j)
      {
	Eigen::Vector3f tmppoint;
	tmppoint = pointcloud.at(j).getArray3fMap();
	Eigen::Vector3f pointp1 = p1 - tmppoint;
	Eigen::Vector3f pointp2 = p2 - tmppoint;
	Eigen::Vector3f pointp3 = p3 - tmppoint;
	Eigen::Vector3f pointp1pointp2normal = geomentry.CrossProduct(pointp1, pointp2);
	Eigen::Vector3f p1p2p2p3normal = geomentry.CrossProduct(p1p2,p2p3);
	float normalangle = geomentry.VectorAngle(pointp1pointp2normal, p1p2p2p3normal);
	
	if(abs(normalangle - 3.14) < 0.01 || normalangle < 0.01)
	{
	  count1 ++;
	  float anglep1pointp2 = geomentry.VectorAngle(pointp1, pointp2);
	  float anglep1pointp3 = geomentry.VectorAngle(pointp1, pointp3);
	  float anglep2pointp3 = geomentry.VectorAngle(pointp2, pointp3);
	  if(abs(3.14*2 - (anglep1pointp2+anglep1pointp3+anglep2pointp3)) < 0.01)
	  {
	    inner.push_back(pointcloud.at(j));
	  }
	}
      }

      if(inner.size() == count1)
      {
	pcl::PointCloud<PointType> holdingpoints = holdingPoints(p1, p2, p3, inner);
	pcl::PointCloud<PointType> ajustedpoints = ajustHoldingPoints(pointcloud, holdingpoints);
	
	
	
	pcl::visualization::PointCloudColorHandlerCustom<PointType> single_colorinner(inner.makeShared(), 255, 0, 0);
	mainview.addPointCloud (inner.makeShared(), single_colorinner, "single_colorinner");
	mainview.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "single_colorinner");
	pcl::visualization::PointCloudColorHandlerCustom<PointType> single_colorholdingpoints(holdingpoints.makeShared(), 0, 0, 0);
	mainview.addPointCloud (holdingpoints.makeShared(), single_colorholdingpoints, "single_colorholdingpoints");
	mainview.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "single_colorholdingpoints");
	pcl::visualization::PointCloudColorHandlerCustom<PointType> single_colorajustedpoints(ajustedpoints.makeShared(), 0, 0, 0);
	mainview.addPointCloud (ajustedpoints.makeShared(), single_colorajustedpoints, "single_colorajustedpoints");
	mainview.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "single_colorajustedpoints");
	
	mainview.addLine(sphere1.at(j1), sphere2.at(j2), "line1"+i);
	mainview.addLine(sphere1.at(j1), sphere2.at(j3), "line2"+i);
	mainview.addLine(sphere2.at(j2), sphere2.at(j3), "line3"+i);
	mainview.spinOnce(10000);
	mainview.removeShape("line1"+i);
	mainview.removeShape("line2"+i);
	mainview.removeShape("line3"+i);
	mainview.removePointCloud("single_colorinner");
	mainview.removePointCloud("single_colorholdingpoints");
	mainview.removePointCloud("single_colorajustedpoints");
      }
    }
  }
  
  while (!mainview.wasStopped ())
  {
    mainview.spinOnce ();
  }
  
  
  
  return 1;
}