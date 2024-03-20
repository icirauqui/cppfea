#include "vis_fem.hpp"


std::pair<Eigen::Vector3d, Eigen::Vector3d> QuaternionLine(
  Eigen::Vector4d qvec, Eigen::Vector3d point, double radius) {

  Eigen::Quaterniond q(qvec(0), qvec(1), qvec(2), qvec(3));
  Eigen::Matrix3d rotation_matrix = q.normalized().toRotationMatrix();

  Eigen::Vector3d reference_vector(1.0, 0.0, 0.0);

  Eigen::Vector3d direction = rotation_matrix * reference_vector;
  direction.normalize();

  Eigen::Vector3d endpoint = point + direction;

  return std::make_pair(point, endpoint);
}


void ViewMesh(bool extrusion,
                   int wait,
                   FEM* fems1,
                   FEM* fems2) {

  std::cout << " View mesh in " << std::endl;

  pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));

  //pcl::visualization::PCLVisualizer viewer;
  std::cout << " Viewer created" << std::endl;
  viewer->setBackgroundColor(1, 1, 1);
  std::cout << " Background color set" << std::endl;

  if (fems1 == nullptr) {
    std::cout << " No FEMs to visualize" << std::endl;
    return;
  }

  AddMesh(extrusion, {0.0, 0.7, 0.0}, viewer, fems1, "_0_");

  if (fems2 != nullptr) {
    AddMesh(extrusion, {0.7, 0.0, 0.0}, viewer, fems2, "_1_");
  }


  std::cout << " All meshes added" << std::endl;
  
  //viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2);
  //viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 0.0);

  viewer->addCoordinateSystem(0.5);

  std::cout << " before spin" << std::endl;
  if (wait > 0) {
    viewer->spinOnce(wait*1000);
    return;
  } else {
    viewer->spin();
  }

  std::cout << " after spin" << std::endl;

  if (viewer->wasStopped()) {
    viewer->close();
  }
  //viewer->close();

  std::cout << " Viewer closed" << std::endl;
  // Remove viewer shared pointer
  viewer.reset();
}



void AddMesh(bool extrusion,
                  std::vector<double> color,
                  pcl::visualization::PCLVisualizer::Ptr viewer,
                  FEM* fem,
                  std::string fem_id) {

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ> (fem->GetCloud()));
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZ> (fem->GetCloud2()));
  
  //viewer->addPointCloud(cloud);
  //viewer->addPointCloud(cloud2);

  std::vector<std::vector<unsigned int>> triangles = fem->GetTriangles();
  std::vector<bool> points_alive = fem->GetPointsAlive();

  for (unsigned int i=0; i<triangles.size(); i++) {
    pcl::PointXYZ p0 = cloud->points[triangles[i][0]];
    pcl::PointXYZ p1 = cloud->points[triangles[i][1]];
    pcl::PointXYZ p2 = cloud->points[triangles[i][2]];
    std::string name = "triangle_1" + fem_id + std::to_string(i);

    // Add line
    viewer->addLine<pcl::PointXYZ>(p0, p1, color[0], color[1], color[2], name + "a");
    viewer->addLine<pcl::PointXYZ>(p1, p2, color[0], color[1], color[2], name + "b");
    viewer->addLine<pcl::PointXYZ>(p2, p0, color[0], color[1], color[2], name + "c");
  }

  for (unsigned int i=0; i<points_alive.size(); i++) {
    if (points_alive[i]) {
      continue;
    }
    viewer->addSphere(cloud->points[i], 0.01, 0.0, 0.7, 0.0, "point_not_alive" + fem_id + std::to_string(i));
  }

  if (extrusion) {
    std::vector<double> color2;
    for (unsigned int j=0; j<color.size(); j++) {
      if (color[j] < 0.5) {
        color2.push_back(color[j] + 0.4);
      } else {
        color2.push_back(1.0);
      }
    }

    for (unsigned int i=0; i<triangles.size(); i++) {
      pcl::PointXYZ p0 = cloud2->points[triangles[i][0]];
      pcl::PointXYZ p1 = cloud2->points[triangles[i][1]];
      pcl::PointXYZ p2 = cloud2->points[triangles[i][2]];
      std::string name = "triangle_1" + fem_id + std::to_string(i) + "_extrusion";

      // Add line
      viewer->addLine<pcl::PointXYZ>(p0, p1, color2[0], color2[1], color2[2], name + "_a");
      viewer->addLine<pcl::PointXYZ>(p1, p2, color2[0], color2[1], color2[2], name + "_b");
      viewer->addLine<pcl::PointXYZ>(p2, p0, color2[0], color2[1], color2[2], name + "_c");
    }

    for (unsigned int i=0; i<cloud->points.size(); i++) {
      if (!points_alive[i]) {
        continue;
      }

      pcl::PointXYZ p0 = cloud->points[i];
      pcl::PointXYZ p1 = cloud2->points[i];
      std::string name = "line_1" + fem_id + std::to_string(i);
      viewer->addLine<pcl::PointXYZ>(p0, p1, color2[0], color2[1], color2[2], name);
    }
  }

  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose = fem->GetPose();

  // If pose not zero, then add as a thick point
  if (pose.second.norm() > 0) {
    pcl::PointXYZ p1(pose.second(0), pose.second(1), pose.second(2));
    viewer->addSphere(p1, 0.1, 0.0, 0.7, 0.0, "pose" + fem_id);

    // Add a line in the direction of the quaternion in pose2.first
    std::pair<Eigen::Vector3d, Eigen::Vector3d> line_pts = QuaternionLine(pose.first, pose.second);
    pcl::PointXYZ p1_a(line_pts.first(0), line_pts.first(1), line_pts.first(2));
    pcl::PointXYZ p1_b(line_pts.second(0), line_pts.second(1), line_pts.second(2));
    viewer->addLine<pcl::PointXYZ>(p1_a, p1_b, 0.0, 0.7, 0.0, "pose_line" + fem_id);
  }
}