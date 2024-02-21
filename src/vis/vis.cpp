#include "vis.hpp"

PCLViewer::PCLViewer() {
    initializeViewer();
}

void PCLViewer::initializeViewer() {
  viewer_.reset(new pcl::visualization::PCLVisualizer("3D Viewer"));
  viewer_->setBackgroundColor(1.0, 1.0, 1.0);
  viewer_->addCoordinateSystem(1.0);
  viewer_->initCameraParameters();
  viewer_->setSize(1070, 820);
}

void PCLViewer::AddNodes(std::vector<Eigen::Vector3d> &pts, std::string name, Eigen::Vector3d color) {
  cloud_names_.push_back(name);
  cloud_colors_.push_back(color);

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<pcl::PointXYZ>);
  cloud->width = pts.size();
  cloud->height = 1;
  cloud->is_dense = false;
  cloud->points.resize(cloud->width * cloud->height);

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_contour = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<pcl::PointXYZ>);
  cloud_contour->width = pts.size();
  cloud_contour->height = 1;
  cloud_contour->is_dense = false;
  cloud_contour->points.resize(cloud_contour->width * cloud_contour->height);

  contour_bound_min_ = Eigen::Vector3d(1e10, 1e10, 1e10);
  contour_bound_max_ = Eigen::Vector3d(-1e10, -1e10, -1e10);
  contour_nodes_.resize(pts.size());

  for (size_t i = 0; i < cloud->points.size(); ++i) {
    cloud->points[i].x = pts[i](0);
    cloud->points[i].y = pts[i](1);
    cloud->points[i].z = pts[i](2);

    for (unsigned int j=0; j<3; j++) {
      if (pts[i](j) < contour_bound_min_(j)) contour_bound_min_(j) = pts[i](j);
      if (pts[i](j) > contour_bound_max_(j)) contour_bound_max_(j) = pts[i](j);
    }
  }

  for (size_t i = 0; i < cloud->points.size(); ++i) {
    if (cloud->points[i].x == contour_bound_min_(0) || cloud->points[i].x == contour_bound_max_(0) ||
        cloud->points[i].y == contour_bound_min_(1) || cloud->points[i].y == contour_bound_max_(1) ||
        cloud->points[i].z == contour_bound_min_(2) || cloud->points[i].z == contour_bound_max_(2)) {
      cloud_contour->points[i].x = cloud->points[i].x;
      cloud_contour->points[i].y = cloud->points[i].y;
      cloud_contour->points[i].z = cloud->points[i].z;
      contour_nodes_[i] = true;
    } else {
      cloud_contour->points[i].x = NAN;
      cloud_contour->points[i].y = NAN;
      cloud_contour->points[i].z = NAN;
      contour_nodes_[i] = false;
    }
  }

  clouds_.push_back(cloud);
  clouds_countour_.push_back(cloud_contour);

  //if (contours_only_)
  //  viewer_->addPointCloud<pcl::PointXYZ>(cloud_contour, name);
  //else
  //  viewer_->addPointCloud<pcl::PointXYZ>(cloud, name);
  //viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, name);
  //viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, color(0), color(1), color(2), name);
}

void PCLViewer::AddEdges(std::vector<std::vector<unsigned int>> &elts, std::string elt) {

  // Create a vector of tuples with pairs of point indices
  std::vector<std::tuple<unsigned int, unsigned int>> edges;
  
  if (elt == "C2D4") {
    edges = {
      std::make_tuple(0, 1),
      std::make_tuple(1, 2),
      std::make_tuple(2, 3),
      //std::make_tuple(3, 0)
    };
  } else if (elt == "C3D6") {
    edges = {
      std::make_tuple(0, 1),
      std::make_tuple(1, 2),
      std::make_tuple(2, 0),
      std::make_tuple(3, 4),
      std::make_tuple(4, 5),
      std::make_tuple(5, 3),
      std::make_tuple(0, 3),
      std::make_tuple(1, 4),
      std::make_tuple(2, 5)
    };
  } else if (elt == "C3D8") {
    edges = {
      std::make_tuple(0, 1),
      std::make_tuple(1, 2),
      std::make_tuple(2, 3),
      std::make_tuple(3, 0),
      std::make_tuple(4, 5),
      std::make_tuple(5, 6),
      std::make_tuple(6, 7),
      std::make_tuple(7, 4),
      std::make_tuple(0, 4),
      std::make_tuple(1, 5),
      std::make_tuple(2, 6),
      std::make_tuple(3, 7)
    };
  }

  for (unsigned int e=0; e<elts.size(); e++) {
    //std::string name = "edge_" + name + "_" + std::to_string(e) + "_";
    for (auto edge : edges) {
      unsigned int i0 = std::get<0>(edge);
      unsigned int i1 = std::get<1>(edge);
      unsigned int node0 = elts[e][i0];
      unsigned int node1 = elts[e][i1];

      std::tuple<unsigned int, unsigned int> edge_pair = std::make_tuple(node0, node1);
      if (std::find(edges_.begin(), edges_.end(), edge_pair) == edges_.end()) {
        edges_.push_back(edge_pair);
      }

      if (!contour_nodes_[node0] || !contour_nodes_[node1]) {
        if (std::find(edges_contour_.begin(), edges_contour_.end(), edge_pair) == edges_contour_.end()) {
          edges_contour_.push_back(edge_pair);
        } 
      }
    }
  }
}

void PCLViewer::AddBCs(std::vector<bool> &nodes, std::vector<std::vector<double>> &mag, bool isLoad) {
  for (unsigned int i=0; i<nodes.size(); i++) {
    if (!nodes[i]) continue;

    load_nodes_.push_back(i);
    load_mags_.push_back(mag[i]);
    is_load_.push_back(isLoad);
  }
}

void PCLViewer::Render(double scale, bool contours_only) {
  if (contours_only) {
    for (unsigned int i=0; i<clouds_countour_.size(); i++) {
      viewer_->addPointCloud<pcl::PointXYZ>(clouds_countour_[i], cloud_names_[i]);
      
      for (auto edge : edges_contour_) {
        unsigned int node0 = std::get<0>(edge);
        unsigned int node1 = std::get<1>(edge);
        viewer_->addLine(clouds_countour_[i]->points[node0], 
                         clouds_countour_[i]->points[node1], 
                         cloud_colors_[i](0), cloud_colors_[i](1), cloud_colors_[i](2),
                         cloud_names_[i] + std::to_string(node0)+"_"+std::to_string(node1));
      }
    }
  } else {
    for (unsigned int i=0; i<clouds_.size(); i++) {
      viewer_->addPointCloud<pcl::PointXYZ>(clouds_[i], cloud_names_[i]);
     
      for (auto edge : edges_) {
        unsigned int node0 = std::get<0>(edge);
        unsigned int node1 = std::get<1>(edge);
        viewer_->addLine(clouds_[i]->points[node0], 
                         clouds_[i]->points[node1], 
                         cloud_colors_[i](0), cloud_colors_[i](1), cloud_colors_[i](2),
                         cloud_names_[i] + std::to_string(node0)+"_"+std::to_string(node1));
      }
    }
  }

  for (unsigned int i=0; i<clouds_.size(); i++) {
    viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, cloud_names_[i]);
    viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, cloud_colors_[i](0), cloud_colors_[i](1), cloud_colors_[i](2), cloud_names_[i]);
  }

  for (unsigned int i=0; i<load_nodes_.size(); i++) {
    unsigned int node = load_nodes_[i];
    std::vector<double> mag = load_mags_[i];
    if (mag[0] == 0.0 && mag[1] == 0.0 && mag[2] == 0.0) {
      viewer_->addSphere<pcl::PointXYZ>(clouds_[0]->points[node], 0.1, 0.4, 0.4, 0.4, "load_" + std::to_string(node));
    } else {
      if (is_load_[i]) {
        viewer_->addArrow<pcl::PointXYZ, pcl::PointXYZ>(pcl::PointXYZ(clouds_[0]->points[node].x + scale * mag[0],
                                                                      clouds_[0]->points[node].y + scale * mag[1],
                                                                      clouds_[0]->points[node].z + scale * mag[2]),
                                                        clouds_[0]->points[node],
                                                        1.0, 0.0, 0.0, false, "load_" + std::to_string(node));
      }
      else {
        viewer_->addSphere<pcl::PointXYZ>(clouds_[0]->points[node], 0.1, 0.8, 0.8, 0.0, "displ_node_" + std::to_string(node));
        viewer_->addArrow<pcl::PointXYZ, pcl::PointXYZ>(pcl::PointXYZ(clouds_[0]->points[node].x + scale * mag[0],
                                                                      clouds_[0]->points[node].y + scale * mag[1],
                                                                      clouds_[0]->points[node].z + scale * mag[2]),
                                                        clouds_[0]->points[node],
                                                        0.8, 0.8, 0.0, false, "displ_mag_" + std::to_string(node));
      }
    }
  }
  
  viewer_->spin();
}