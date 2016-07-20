#include "Geomentry.h"

Geomentry::Geomentry()
{

}

Geomentry::~Geomentry()
{

}


cv::Mat Geomentry::findCenter(cv::Mat points)
{
    cv::Mat center = cv::Mat::zeros(3, 1, CV_32F);
    for(int i = 0; i < points.cols; i ++)
    {
	center.at<float>(0,0) += points.at<float>(0,i);
	center.at<float>(1,0) += points.at<float>(1,i);
	center.at<float>(2,0) += points.at<float>(2,i);
    }
    
    center /= points.cols;
    return center;
}

float Geomentry::computescale_s(const cv::Mat &pointsX, const cv::Mat &pointsY, const cv::Mat &R_h)
{
    if(pointsX.cols != pointsY.cols)
    {
	cout << "pointsX.size is not equal to pointsY.size, system exit!!" << endl;
	exit(1);
    }
    cv::Mat centerY = findCenter(pointsY);
    cv::Mat centerX = findCenter(pointsX);
    cv::Mat pointsXtmp = cv::Mat(pointsX.rows, pointsX.cols, CV_32F);
    cv::Mat pointsYtmp = cv::Mat(pointsY.rows, pointsY.cols, CV_32F);
    for(int i = 0; i < pointsX.cols; i ++)
    {
	pointsXtmp.at<float>(0,i) = pointsX.at<float>(0,i);
	pointsXtmp.at<float>(1,i) = pointsX.at<float>(1,i);
	pointsXtmp.at<float>(2,i) = pointsX.at<float>(2,i);
	
	pointsYtmp.at<float>(0,i) = pointsY.at<float>(0,i);
	pointsYtmp.at<float>(1,i) = pointsY.at<float>(1,i);
	pointsYtmp.at<float>(2,i) = pointsY.at<float>(2,i);
    }
    for(int i = 0; i < pointsY.cols; i ++)
    {
	pointsYtmp.col(i) -= centerY;
	pointsXtmp.col(i) -= centerX;
    }
    
    cv::Mat R_htmp = R_h.t();
    
    float a = cv::trace(pointsYtmp * pointsXtmp.t() * R_htmp )[0] ;
    
    float c = cv::trace(pointsXtmp * pointsXtmp.t())[0];

    float s = a / c;    
    return s;
}

cv::Mat Geomentry::computeS(const cv::Mat &pointsX, const cv::Mat &pointsY, float sinit)
{
    if(pointsX.cols != pointsY.cols)
    {
	cout << "pointsX.size is not equal to pointsY.size, system exit!!" << endl;
	exit(1);
    }
    cv::Mat centerY = findCenter(pointsY);
    cv::Mat centerX = findCenter(pointsX);
    cv::Mat pointsXtmp = cv::Mat(pointsX.rows, pointsX.cols, CV_32F);
    cv::Mat pointsYtmp = cv::Mat(pointsY.rows, pointsY.cols, CV_32F);
    for(int i = 0; i < pointsX.cols; i ++)
    {
	pointsXtmp.at<float>(0,i) = pointsX.at<float>(0,i);
	pointsXtmp.at<float>(1,i) = pointsX.at<float>(1,i);
	pointsXtmp.at<float>(2,i) = pointsX.at<float>(2,i);
	
	pointsYtmp.at<float>(0,i) = pointsY.at<float>(0,i);
	pointsYtmp.at<float>(1,i) = pointsY.at<float>(1,i);
	pointsYtmp.at<float>(2,i) = pointsY.at<float>(2,i);
    }
    for(int i = 0; i < pointsY.cols; i ++)
    {
	pointsYtmp.col(i) -= centerY;
	pointsXtmp.col(i) -= centerX;
    }
    pointsXtmp *=sinit;
    cv::Mat S = pointsXtmp*pointsYtmp.t();
    return S;
}

cv::Mat Geomentry::findRTfromS(cv::Mat Xcorre, cv::Mat Ycorre, float sinit)
{
    cv::Mat S = computeS(Xcorre, Ycorre, sinit);
    
    float Sxx = S.at<float>(0,0);
    float Sxy = S.at<float>(0,1);
    float Sxz = S.at<float>(0,2);
    float Syx = S.at<float>(1,0);
    float Syy = S.at<float>(1,1);
    float Syz = S.at<float>(1,2);
    float Szx = S.at<float>(2,0);
    float Szy = S.at<float>(2,1);
    float Szz = S.at<float>(2,2);
    
    cv::Mat K = cv::Mat(4,4,CV_32F);
    K.at<float>(0,0) = Sxx+Syy+Szz;
    K.at<float>(0,1) = Syz-Szy;
    K.at<float>(0,2) = Szx-Sxz;
    K.at<float>(0,3) = Sxy-Syx;
    K.at<float>(1,0) = Syz-Szy;
    K.at<float>(1,1) = Sxx-Syy-Szz;
    K.at<float>(1,2) = Sxy+Syx;
    K.at<float>(1,3) = Szx+Sxz;
    K.at<float>(2,0) = Szx-Sxz;
    K.at<float>(2,1) = Sxy+Syx;
    K.at<float>(2,2) = Syy-Sxx-Szz;
    K.at<float>(2,3) = Syz+Szy;
    K.at<float>(3,0) = Sxy-Syx;
    K.at<float>(3,1) = Szx+Sxz;
    K.at<float>(3,2) = Syz+Szy;
    K.at<float>(3,3) = Szz-Sxx-Syy;
    
    cv::Mat eigenvalues, eigenvectors;
    cv::eigen(K, eigenvalues, eigenvectors);
    //     cout << eigenvectors << endl;
    //     exit(1);
    //cout << "size of eigenvector of K:" << eigenvectors.rows << " " << eigenvectors.cols << endl;
    float q0 = eigenvectors.at<float>(0, 0);
    float qx = eigenvectors.at<float>(0, 1);
    float qy = eigenvectors.at<float>(0, 2);
    float qz = eigenvectors.at<float>(0, 3);
    
    cv::Mat Rtmp = cv::Mat(3,3,CV_32F);
    // quaternion to rotation matrix
    Rtmp.at<float>(0,0) = q0*q0 + qx*qx - qy*qy - qz*qz;
    Rtmp.at<float>(0,1) = 2 * (qx*qy - q0*qz);
    Rtmp.at<float>(0,2) = 2 * (qx*qz + q0*qy);
    Rtmp.at<float>(1,0) = 2 * (qy*qx + q0*qz);
    Rtmp.at<float>(1,1) = q0*q0 - qx*qx + qy*qy - qz*qz;
    Rtmp.at<float>(1,2) = 2 * (qy*qz - q0*qx);
    Rtmp.at<float>(2,0) = 2 * (qz*qx - q0*qy);
    Rtmp.at<float>(2,1) = 2 * (qz*qy + q0*qx);
    Rtmp.at<float>(2,2) = q0*q0 - qx*qx - qy*qy + qz*qz;
    
    //     cout << centerY << " " << centerXcorres << endl;
    cv::Mat centerXcorres = findCenter(Xcorre);
    cv::Mat centerY = findCenter(Ycorre);
    cv::Mat Ttmp =  centerY/sinit - Rtmp*centerXcorres;
    //computescale_s(Xcorre, Ycorre, Rtmp);
    cv::Mat RT = cv::Mat::zeros(3,4,CV_32F);  
    for(int i = 0; i < 3; i ++)
    {
	for(int j = 0; j < 3; j++)
	{
	    RT.at<float>(i,j) = Rtmp.at<float>(i,j);
	}
	RT.at<float>(i,3) = Ttmp.at<float>(i, 0);
    }        
    return RT;
}

Eigen::Matrix<float, 3, 3> Geomentry::RotationAboutVector(Eigen::Vector3f rotationAxis, float theta)
{
  Eigen::Matrix<float, 3, 3> rotationMatrix;
  float norm = rotationAxis.norm();
  Eigen::Vector3f normlizedRotationAxis = rotationAxis / norm;
  rotationMatrix(0, 0) = cos(theta)+normlizedRotationAxis(0)*normlizedRotationAxis(0)*(1-cos(theta));
  rotationMatrix(0, 1) = normlizedRotationAxis(0)*normlizedRotationAxis(1)*(1-cos(theta))-normlizedRotationAxis(2)*sin(theta);
  rotationMatrix(0, 2) = normlizedRotationAxis(1)*sin(theta)+normlizedRotationAxis(0)*normlizedRotationAxis(2)*(1-cos(theta));
  
  rotationMatrix(1, 0) = normlizedRotationAxis(2)*sin(theta)+normlizedRotationAxis(0)*normlizedRotationAxis(1)*(1-cos(theta));
  rotationMatrix(1, 1) = cos(theta)+normlizedRotationAxis(1)*normlizedRotationAxis(1)*(1-cos(theta));
  rotationMatrix(1, 2) = -normlizedRotationAxis(0)*sin(theta)+normlizedRotationAxis(1)*normlizedRotationAxis(2)*(1-cos(theta));
  
  rotationMatrix(2, 0) = -normlizedRotationAxis(1)*sin(theta)+normlizedRotationAxis(0)*normlizedRotationAxis(2)*(1-cos(theta));
  rotationMatrix(2, 1) = normlizedRotationAxis(0)*sin(theta)+normlizedRotationAxis(1)*normlizedRotationAxis(2)*(1-cos(theta));
  rotationMatrix(2, 2) = cos(theta)+normlizedRotationAxis(2)*normlizedRotationAxis(2)*(1-cos(theta));
  
  if(abs(rotationMatrix.determinant() - 1) < 0.0001)
    return rotationMatrix;
  else 
  {
    cerr << "the determinat of rotation matirx is not equal to 1, system exit(1)!!" << endl;
    exit(1);
  }
}


float Geomentry::VectorAngle(Eigen::Vector3f veca, Eigen::Vector3f vecb)
{
  float noma = sqrt(veca[0]*veca[0]+veca[1]*veca[1]+veca[2]*veca[2]);
  float nomb = sqrt(vecb[0]*vecb[0]+vecb[1]*vecb[1]+vecb[2]*vecb[2]);
  float costheta = (veca[0]*vecb[0] + veca[1]*vecb[1] + veca[2]*vecb[2])/(noma*nomb);
  if(costheta > 1.0000000)
  costheta = 0.999999999;
  if(costheta < -1.0000000)
  costheta = -0.99999999;	
  return acos(costheta);
}

float Geomentry::DotProduct(Eigen::Vector3f veca, Eigen::Vector3f vecb)
{
  float noma = sqrt(veca[0]*veca[0]+veca[1]*veca[1]+veca[2]*veca[2]);
  float nomb = sqrt(vecb[0]*vecb[0]+vecb[1]*vecb[1]+vecb[2]*vecb[2]);
  float dotproduct = (veca[0]*vecb[0] + veca[1]*vecb[1] + veca[2]*vecb[2]);	
  return dotproduct;
}


Eigen::Vector3f Geomentry::CrossProduct(Eigen::Vector3f veca, Eigen::Vector3f vecb)
{
  Eigen::Vector3f crossVec;
  crossVec[0] = veca[1]*vecb[2] - vecb[1]*veca[2];
  crossVec[1] = veca[2]*vecb[0] - vecb[2]*veca[0];
  crossVec[2] = veca[0]*vecb[1] - vecb[0]*veca[1];
  return crossVec;
}

pcl::PointCloud< PointType > Geomentry::readPointCloud(string pointcloudFileName)
{
  pcl::PointCloud<PointType> pointcloud;
  if(pointcloudFileName.find(".ply") == pointcloudFileName.length() - 4)
  {
    if (pcl::io::loadPLYFile<PointType> (pointcloudFileName, pointcloud) == -1) //* load the file
    {
      cerr << "can not read file " << pointcloudFileName << endl;
      exit(1);
    }
  }
  if(pointcloudFileName.find(".pcd") == pointcloudFileName.length() - 4)
  {
    if (pcl::io::loadPCDFile<PointType> (pointcloudFileName, pointcloud) == -1) //* load the file
    {
      cerr << "can not read file " << pointcloudFileName << endl;
      exit(1);
    }
  }
  Eigen::Vector4f sensor_origin;  
  Eigen::Quaternion<float> sensor_orientation;  
  sensor_origin = Eigen::Vector4f::Zero();  
  sensor_orientation = Eigen::Quaternionf::Identity ();
  pointcloud.sensor_origin_ = sensor_origin;
  pointcloud.sensor_orientation_ = sensor_orientation;
  return pointcloud;
}

pcl::PointCloud< PointType > Geomentry::calculateNormal(pcl::PointCloud< PointType > pointcloud, 
						      float radious)
{
  
  pcl::PointXYZ centerpoint1;
  centerpoint1.x = 0;
  centerpoint1.y = 0;
  centerpoint1.z = 0;
  #pragma omp parallel for
  for(int i = 0; i < pointcloud.size(); ++i)
  {
    centerpoint1.x += pointcloud.at(i).x;
    centerpoint1.y += pointcloud.at(i).y;
    centerpoint1.z += pointcloud.at(i).z;
  }
  centerpoint1.x /= pointcloud.size();
  centerpoint1.y /= pointcloud.size();
  centerpoint1.z /= pointcloud.size();
		
  pcl::PointCloud<pcl::Normal> normal;
  cout << "calculating normal" << endl;
  pcl::NormalEstimationOMP<PointType, pcl::Normal> ne1;
  ne1.setInputCloud (pointcloud.makeShared());
  ne1.setNumberOfThreads(8);
  // Create an empty kdtree representation, and pass it to the normal estimation object.
  // Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
  pcl::search::KdTree<PointType>::Ptr tree (new pcl::search::KdTree<PointType> ());
  ne1.setSearchMethod (tree);
  // Use all neighbors in a sphere of radius 3cm
  ne1.setRadiusSearch (radious);
  ne1.setViewPoint (centerpoint1.x, centerpoint1.y, centerpoint1.z);
  // Compute the features
  ne1.compute (normal);
  for(int i = 0; i < pointcloud.size(); ++i)
  {
    pointcloud.at(i).normal_x = normal.at(i).normal_x;
    pointcloud.at(i).normal_y = normal.at(i).normal_y;
    pointcloud.at(i).normal_z = normal.at(i).normal_z;
  }
  return pointcloud;
}


double Geomentry::computeCloudResolution(const pcl::PointCloud<PointType>::ConstPtr& cloud)
{
	double resolution = 0.0;
	int numberOfPoints = 0;
	int nres;
	std::vector<int> indices(2);
	std::vector<float> squaredDistances(2);
	pcl::search::KdTree<PointType> tree;
	tree.setInputCloud(cloud);

	for (size_t i = 0; i < cloud->size(); ++i)
	{
		if (! pcl_isfinite((*cloud)[i].x))
			continue;

		// Considering the second neighbor since the first is the point itself.
		nres = tree.nearestKSearch(i, 2, indices, squaredDistances);
		if (nres == 2)
		{
			resolution += sqrt(squaredDistances[1]);
			++numberOfPoints;
		}
	}
	if (numberOfPoints != 0)
		resolution /= numberOfPoints;

	return resolution;
}

pcl::PointCloud< PointType > Geomentry::computeBarycenterC(const pcl::PointCloud< PointType >::ConstPtr& cloud)
{
  pcl::PointCloud< PointType > barycenterC;
  PointType tmppoint;
  tmppoint.x = 0;
  tmppoint.y = 0;
  tmppoint.z = 0;
  for(int i = 0; i < cloud->size(); ++i)
  {    
    tmppoint.x += cloud->at(i).x;
    tmppoint.y += cloud->at(i).y;
    tmppoint.z += cloud->at(i).z;
  }
  tmppoint.x /= cloud->size();
  tmppoint.y /= cloud->size();
  tmppoint.z /= cloud->size();
  barycenterC.push_back(tmppoint);
  return barycenterC;
}

float Geomentry::computeRadius(pcl::PointCloud< PointType > BarycenterC, pcl::PointCloud< PointType >::ConstPtr cloud)
{
  float PointcloudRadius = 0;
  float minidis = 0;
  for(int i = 0; i < cloud->size(); ++ i)
  {
    Eigen::Vector3f veca = cloud->at(i).getArray3fMap()-BarycenterC.at(0).getArray3fMap();
    float tmpdis = veca.norm();
    if(tmpdis > minidis)
    {
      minidis = tmpdis;
    }
  }
  PointcloudRadius = minidis;
  return PointcloudRadius;
}

pcl::PointCloud< PointType > Geomentry::generateSphere(pcl::PointCloud< PointType > BarycenterC, float radius)
{
  pcl::PointCloud< PointType > spherepoints;
  //rotation about z axis
  Eigen::Vector3f rotationaxis;
  rotationaxis[0] = 0;
  rotationaxis[1] = 0;
  rotationaxis[2] = 1;
  int step = 30;
  float theta = 3.14 / step;
  Eigen::Vector3f BarycenterEigen = BarycenterC.at(0).getArray3fMap();
  Eigen::Vector3f startpoint;
  startpoint[0] = -radius;
  startpoint[1] = 0;
  startpoint[2] = 0;
  for(int i = 0; i <= step; ++ i)
  {
    Eigen::Matrix<float, 3, 3> rotationMatrix = this->RotationAboutVector(rotationaxis, theta*i);
    Eigen::Vector3f newpoint = rotationMatrix * startpoint;
    PointType point;
    point.x = newpoint[0];
    point.y = newpoint[1];
    point.z = newpoint[2];
    spherepoints.push_back(point);
  }
  pcl::PointCloud< PointType > sphere;
  //rotation about x axis
  rotationaxis[0] = 1;
  rotationaxis[1] = 0;
  rotationaxis[2] = 0;
  for(int i = 0; i <= 2*step; ++ i)
  {
    Eigen::Matrix<float, 3, 3> rotationMatrix = this->RotationAboutVector(rotationaxis, theta*i);
    for(int j = 0; j < spherepoints.size(); ++ j)
    {
      Eigen::Vector3f onepoint = spherepoints.at(j).getArray3fMap();
      Eigen::Vector3f newpoint = rotationMatrix * onepoint;
      PointType point;
      point.x = newpoint[0];
      point.y = newpoint[1];
      point.z = newpoint[2];
      sphere.push_back(point);
    }
  }
  
  for(int i = 0; i < sphere.size(); ++ i)
  {
    sphere.at(i).x += BarycenterC.at(0).x;
    sphere.at(i).y += BarycenterC.at(0).y;
    sphere.at(i).z += BarycenterC.at(0).z;
  }
  return sphere;
}
