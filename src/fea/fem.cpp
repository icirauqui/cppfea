#include "fem.hpp"



std::pair<Eigen::Vector4d, Eigen::Vector3d> ApproximatePose(std::vector<Eigen::Vector3d> pts) {
  Eigen::Vector3d centroid(0.0, 0.0, 0.0);
  for (int i = 0; i < pts.size(); i++) {
    centroid += pts[i];
  }
  centroid /= pts.size();

  double dist = sqrt(pow(centroid(0), 2) + pow(centroid(1), 2) + pow(centroid(2), 2));
  Eigen::Vector3d tvec = centroid + Eigen::Vector3d(0, 0, dist/2);

  // Compute rotation, from t to centroid, as a quaternion
  Eigen::Vector3d direction = centroid - tvec;
  direction.normalize();
  Eigen::Quaterniond quaternion;
  quaternion.setFromTwoVectors(-Eigen::Vector3d::UnitX(), direction);
  Eigen::Vector4d qvec = quaternion.coeffs();

  return std::make_pair(qvec, tvec);
}


FEM::FEM(std::string element): element_(element) {
  colorsf_.push_back({0.0, 0.7, 0.0});
  colorsf_.push_back({0.7, 0.0, 0.0});
  colorsf_.push_back({0.0, 0.0, 0.7});
  colorsb_.push_back({0.7, 0.9, 0.7});
  colorsb_.push_back({0.9, 0.7, 0.7});
  colorsb_.push_back({0.7, 0.7, 0.9});

}



FEM::FEM(FEM& fem) {
  element_ = fem.element_;
  points_ = fem.points_;
  points_alive_ = fem.points_alive_;
  pc_ = fem.pc_;
  points_indices_ = fem.points_indices_;
  mls_indices_ = fem.mls_indices_;
  normals_ = fem.normals_;
  mesh_ = fem.mesh_;
  triangles_ = fem.triangles_;
  quadrilaterals_ = fem.quadrilaterals_;
  indices_not_triangulated_ = fem.indices_not_triangulated_;
  points2_ = fem.points2_;
  pc2_ = fem.pc2_;
  elements_ = fem.elements_;
  transform_into_quads_ = fem.transform_into_quads_;
  element_height_ = fem.element_height_;
  pose_ = fem.pose_;
  colorsf_ = fem.colorsf_;
  colorsb_ = fem.colorsb_;
}


void FEM::Replicate(FEM* fem) {
  
  
  InitCloud();

  // copy mls_indices_ and mls_points_
  mls_indices_ = fem->mls_indices_;
  points_indices_ = fem->points_indices_;
  //pc_.width = fem->pc_.width;
  //pc_.height = fem->pc_.height;
  //pc_.points.resize(pc_.width * pc_.height);
  //for (size_t i = 0; i < pc_.points.size(); ++i) {
  //  pc_.points[i].x = fem->pc_.points[i].x;
  //  pc_.points[i].y = fem->pc_.points[i].y;
  //  pc_.points[i].z = fem->pc_.points[i].z;
  //}

  // copy normals_, triangles_, points_indices_, quadrilaterals_
  normals_.clear();
  for (unsigned int i=0; i<fem->normals_.size(); i++) {
    normals_.push_back(Eigen::Vector3d(fem->normals_[i](0), fem->normals_[i](1), fem->normals_[i](2)));
  }
  triangles_.clear();
  for (unsigned int i=0; i<fem->triangles_.size(); i++) {
    std::vector<unsigned int> triangle;
    for (unsigned int j=0; j<fem->triangles_[i].size(); j++) {
      triangle.push_back(fem->triangles_[i][j]);
    }
    triangles_.push_back(triangle);
  }
  points_indices_ = fem->points_indices_;
  pc_.width = fem->pc_.width;
  pc_.height = fem->pc_.height;
  pc_.points.resize(pc_.width * pc_.height);
  for (size_t i = 0; i < pc_.points.size(); ++i) {
    pc_.points[i].x = fem->pc_.points[i].x;
    pc_.points[i].y = fem->pc_.points[i].y;
    pc_.points[i].z = fem->pc_.points[i].z;
  }
  quadrilaterals_.clear();
  for (unsigned int i=0; i<fem->quadrilaterals_.size(); i++) {
    std::vector<unsigned int> quadrilateral;
    for (unsigned int j=0; j<fem->quadrilaterals_[i].size(); j++) {
      quadrilateral.push_back(fem->quadrilaterals_[i][j]);
    }
    quadrilaterals_.push_back(quadrilateral);
  }



  // copy element_height_, pose_
  points2_.clear();
  pc2_.width = fem->pc2_.width;
  pc2_.height = fem->pc2_.height;
  pc2_.points.resize(pc2_.width * pc2_.height);


  // Average normal
  Eigen::Vector3d normal = Eigen::Vector3d::Zero();
  for (Eigen::Vector3d n : normals_) {
    normal -= n;
  }
  normal /= normals_.size();    

  // truncate normals to 3 decimals
  double scale = std::pow(10.0, 3);
  for (unsigned int i=0; i<3; i++) {
    normal(i) = std::round(normal(i) * scale) / scale;
  }

  element_height_ = fem->element_height_;

  //std::cout << "Element height: " << element_height_ << std::endl;
  //int cloud_idx = 0;
  //for (unsigned int i=0; i<points_.size(); i++) {
  //  Eigen::Vector3d point2 = points_[i] - element_height_/2 * normal;
  //  points2_.push_back(point2);
  //  if (points_alive_[i]) {
  //    pc2_.points[cloud_idx].x = point2(0);
  //    pc2_.points[cloud_idx].y = point2(1);
  //    pc2_.points[cloud_idx].z = point2(2);
  //    cloud_idx++;
  //  }
  //}
  //std::cout << "points2_ built" << std::endl;

  //std::cout << "points size: " << points_.size() << std::endl;
  //std::cout << pc2_.width << " " << pc2_.height << std::endl;

  int cloud_idx = 0;
  for (size_t i = 0; i < pc_.points.size(); ++i) {
    if (!points_alive_[i]) {
      continue;
    }
    Eigen::Vector3d point2 = points_[points_indices_[i]] - element_height_/2 * normal;
    points2_.push_back(point2);
    pc2_.points[cloud_idx].x = point2(0);
    pc2_.points[cloud_idx].y = point2(1);
    pc2_.points[cloud_idx].z = point2(2);
    cloud_idx++;

    //pc2_.points[i].x = fem->pc2_.points[i].x;
    //pc2_.points[i].y = fem->pc2_.points[i].y;
    //pc2_.points[i].z = fem->pc2_.points[i].z;
  }


  elements_.clear();
  for (unsigned int i=0; i<fem->elements_.size(); i++) {
    std::vector<unsigned int> element;
    for (unsigned int j=0; j<fem->elements_[i].size(); j++) {
      element.push_back(fem->elements_[i][j]);
    }
    elements_.push_back(element);
  }

  points_alive_ = fem->points_alive_;

  std::cout << " - Model replicated" << std::endl;


}



void FEM::AddPoint(Eigen::Vector3d point) {
  points_.push_back(point);
  points_alive_.push_back(true);
}


bool FEM::InitCloud() {
  pc_.width = points_.size();
  pc_.height = 1;
  pc_.points.resize(pc_.width * pc_.height);
  points_indices_.clear();

  // simulation
  //points_alive_[4] = false;

  for (size_t i = 0; i < pc_.points.size(); ++i) {
    if (!points_alive_[i]) {
      continue;
    }

    points_indices_.push_back(i);
    pc_.points[i].x = points_[i](0);
    pc_.points[i].y = points_[i](1);
    pc_.points[i].z = points_[i](2);
  }

  // simulation
  //for (int i=0; i<points_indices_.size(); i++) {
  //  std::cout << " " << points_indices_[i];
  //}
  //std::cout << std::endl;

  if (pc_.width == 0) {} else {
    pose_ = ApproximatePose(points_);
    std::cout << "   Approximated pose: " << pose_.first.transpose() << " | " << pose_.second.transpose() << std::endl;
  }

  return pc_.width > 0;
}

bool FEM::MovingLeastSquares(bool simulation) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
    new pcl::PointCloud<pcl::PointXYZ> (pc_));

  // Create a KD-Tree
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);

  // Output has the PointNormal type in order to store the normals calculated by MLS
  pcl::PointCloud<pcl::PointNormal> mls_points;

  // Init object (second point type is for the normals, even if unused)
  pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
 
  mls.setComputeNormals(true);
  mls.setInputCloud (cloud);


  int mls_attempts = 0;
  double search_radious = 0.1;
  
  while (mls_attempts < 10) {
    // Set parameters
    mls.setPolynomialOrder (5);
    mls.setSearchMethod (tree);
    mls.setSearchRadius (search_radious);
    //mls.setUpsamplingMethod (pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal>::NONE);
    //mls.setUpsamplingRadius (search_radious);
    //mls.setUpsamplingStepSize (10);
    // Reconstruct
    mls.process (mls_points);

    if (mls_points.size() >= 0.5 * pc_.points.size()) {
      break;
    } else {
      mls_attempts++;
      search_radious *= 2.0;
      mls_points.clear();
      std::cout << "   MLS failed, trying again with search radius: " << search_radious << std::endl;
    }
  }



  if (mls_points.size() == 0) {
    std::cout << "   MLS failed" << std::endl;
    return false;
  }

  // Get corresponding indexes: for each output point, returns the index of the
  // input one.
  mls_indices_.clear();

  pcl::PointIndicesPtr pIdx1 = mls.getCorrespondingIndices();
  for (unsigned int i = 0; i < pIdx1->indices.size(); i++)
    mls_indices_.push_back(pIdx1->indices[i]);

  //std::cout << std::endl;
  //for (int i = 0; i<points_alive_.size(); i++) {
  //  std::cout << " " << points_alive_[i];
  //}
  //std::cout << std::endl;

  // simulation: erase position 0 in mls_indices_
  if (simulation) {
    std::cout << "   mls_points.size(): " << mls_points.size() << std::endl;
    mls_indices_.erase(mls_indices_.begin() + 0);
    mls_points.erase(mls_points.begin() + 0);
  }
    
  std::vector<int> points_indices;
  int idxit = 0;
  for (unsigned int i = 0; i < points_indices_.size(); i++) {
    if (i == mls_indices_[idxit]) {
      idxit++;
      points_indices.push_back(points_indices_[i]);
    } else {
      points_alive_[points_indices_[i]] = false;
    }
  }
    points_indices_ = points_indices;

  if (simulation) {
    // simulation
    std::cout << "   Simulation MLS:";
    for (int i = 0; i<points_indices_.size(); i++) {
      std::cout << " " << points_indices_[i];
    }
    std::cout << std::endl;
    std::cout << "   Simulation MLS:";
    for (int i = 0; i<points_alive_.size(); i++) {
      std::cout << " " << points_alive_[i];
    }
    std::cout << std::endl;
  }

  std::cout << "   Point Cloud size: " << pc_.points.size() << " -> " << mls_points.size() << std::endl;


  // Replace pc_ with mls_points
  pc_.width = mls_points.size();
  pc_.height = 1;
  pc_.points.resize(pc_.width * pc_.height);

  for (size_t i = 0; i < pc_.points.size(); ++i) {
    pc_.points[i].x = mls_points[i].x;
    pc_.points[i].y = mls_points[i].y;
    pc_.points[i].z = mls_points[i].z;
  }

  return true;
}


bool FEM::Triangulate(bool simulation) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ> (pc_));

  //std::cout << "TRIANGULATE" << std::endl;
  //std::cout << "pc_.points.size() = " << pc_.points.size() << std::endl;

  // Normal estimation
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
  pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
  tree->setInputCloud (cloud);
  n.setInputCloud(cloud);
  n.setSearchMethod(tree);
  n.setKSearch(20);
  n.compute(*normals);

  // Store normals in normals_ as Eigen::Vector3d
  normals_.clear();
  for (unsigned int i = 0; i < normals->size(); i++) {
    Eigen::Vector3d normal;
    normal << normals->points[i].normal_x,
              normals->points[i].normal_y,
              normals->points[i].normal_z;
    normals_.push_back(normal);
  }

  // Concatenate the XYZ and normal fields*
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::PointNormal>);
  pcl::concatenateFields (*cloud, *normals, *cloud_with_normals);

    // Create search tree*
  pcl::search::KdTree<pcl::PointNormal>::Ptr tree2 (new pcl::search::KdTree<pcl::PointNormal>);
  tree2->setInputCloud (cloud_with_normals);

  // Initialize objects
  pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;

  // Set the maximum distance between connected points (maximum edge length)
  gp3.setSearchRadius (1.5);

  // Set typical values for the parameters
  gp3.setMu (2.5);
  gp3.setMaximumNearestNeighbors (100);
  gp3.setMaximumSurfaceAngle(M_PI/4); // 45 degrees
  gp3.setMinimumAngle(M_PI/18); // 10 degrees
  gp3.setMaximumAngle(2*M_PI/3); // 120 degrees
  gp3.setNormalConsistency(false);

  // Get result
  gp3.setInputCloud (cloud_with_normals);
  gp3.setSearchMethod (tree2);
  gp3.reconstruct (mesh_);

  // Additional vertex information
  // To which part (cloud) does each vertex belong
  std::vector<int> parts = gp3.getPartIDs();
  // Whether the vertex status is [-1,0,1,2,3] = [NONE,FREE,FRINGE,BOUNDARY,COMPLETED]
  std::vector<int> states = gp3.getPointStates();

  // Get largest part
  std::unordered_map<int, int> map;
  int max_val = -1;
  int max_key = -1;
  for (unsigned int i=0; i<parts.size(); i++) {
    if (map.find(parts[i]) == map.end()) {
      map[parts[i]] = 1;
      if (max_val < 1) {
        max_val = 1;
        max_key = parts[i];
      }
    }
    else {
      map[parts[i]]++;
      if (map[parts[i]] > max_val) {
        max_val = map[parts[i]];
        max_key = parts[i];
      }
    }
  }

  // Get triangles of largest part only
  triangles_.clear();
  for (unsigned int i=0; i<mesh_.polygons.size(); i++) {
    unsigned int nver0 = mesh_.polygons[i].vertices[0];
    unsigned int nver1 = mesh_.polygons[i].vertices[1];
    unsigned int nver2 = mesh_.polygons[i].vertices[2];

    if (parts[nver0] == max_key && parts[nver1] == max_key && parts[nver2] == max_key) {
      std::vector<unsigned int> triangle;
      triangle.push_back(nver0);
      triangle.push_back(nver1);
      triangle.push_back(nver2);
      triangles_.push_back(triangle);
    }
  }

  // simulation: remove the last triangle
  if (simulation) {
    std::cout << "   Number of triangles: " << triangles_.size() << std::endl;
    for (unsigned int i=0; i<3; i++) {
      std::cout << " " << triangles_[triangles_.size() - 3][i];
    }
    std::cout << std::endl;
    for (unsigned int i=0; i<3; i++) {
      std::cout << " " << triangles_[triangles_.size() - 2][i];
    }
    std::cout << std::endl;
    for (unsigned int i=0; i<3; i++) {
      std::cout << " " << triangles_[triangles_.size() - 1][i];
    }
    std::cout << std::endl;
    triangles_.erase(triangles_.end() - 1);
    triangles_.erase(triangles_.end() - 1);
    triangles_.erase(triangles_.end() - 1);
    std::cout << "   Number of triangles: " << triangles_.size() << std::endl;
  }

  // Update points_alive_ and points_indices_
  std::vector<unsigned int> triangle_indices;
  for (auto t: triangles_) {
    for (auto i: t) {
      auto it = std::find(triangle_indices.begin(), triangle_indices.end(), i);
      if (it == triangle_indices.end()) {
        triangle_indices.push_back(i);
      }
    }
  }

  // sort ascending triangle_indices
  std::sort(triangle_indices.begin(), triangle_indices.end());

  //for (unsigned int i=0; i<triangle_indices.size(); i++) {
  //  std::cout << " " << triangle_indices[i];
  //}
  //std::cout << std::endl;



  std::vector<int> missing_indices, points_indices, points_indices_internal;
  int internal_idx = 0;
  for (unsigned int i=0; i<pc_.points.size(); i++) {
    auto it = std::find(triangle_indices.begin(), triangle_indices.end(), i);
    if (it == triangle_indices.end()) {
      missing_indices.push_back(i);
      points_alive_[points_indices_[i]] = false;
      points_indices_internal.push_back(-1);
    } else {
      points_indices.push_back(points_indices_[i]);
      points_indices_internal.push_back(internal_idx);
      internal_idx++;
    }
  }

  //std::cout << "Missing indices: ";
  //for (unsigned int i=0; i<missing_indices.size(); i++) {
  //  std::cout << " " << missing_indices[i];
  //}
  //std::cout << std::endl;

  //std::cout << "Points indices: ";
  //for (unsigned int i=0; i<points_indices.size(); i++) {
  //  std::cout << " " << points_indices[i];
  //}
  //std::cout << std::endl;

  //std::cout << "Triangulation: ";
  //for (unsigned int i=0; i<points_alive_.size(); i++) {
  //  std::cout << " " << points_alive_[i];
  //}
  //std::cout << std::endl;
  std::cout << "   Number of indices / points: " << points_indices.size() << " / " << points_.size() << std::endl;
  std::cout << "   Number of triangles: " << triangles_.size() << std::endl;

  points_indices_ = points_indices;

  // Replace pc_ with triangulated points
  pc_.width = points_indices_.size();
  pc_.height = 1;
  pc_.points.resize(pc_.width * pc_.height);

  for (size_t i = 0; i < pc_.points.size(); ++i) {
    pc_.points[i].x = points_[points_indices_[i]](0);
    pc_.points[i].y = points_[points_indices_[i]](1);
    pc_.points[i].z = points_[points_indices_[i]](2);
  }

  // Finally update triangles_ with internal indices
  for (unsigned int i=0; i<triangles_.size(); i++) {
    for (unsigned int j=0; j<triangles_[i].size(); j++) {
      triangles_[i][j] = points_indices_internal[triangles_[i][j]];
    }
  }

  //std::cout << "pc_.points.size() = " << pc_.points.size() << std::endl;


  
  return triangles_.size() > 0 ? true : false;
}



bool FEM::TriangulateByProjection(bool simulation) {

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ> (pc_));


  // Use pose_ to project points_ onto a plane
  // pose is composed of (Eigen::Vector4d, Eigen::Vector3d) -> rotation quaternion and translation vector
  pcl::PointCloud<pcl::PointXYZ>::Ptr projected_cloud(new pcl::PointCloud<pcl::PointXYZ>);

  Eigen::Vector4d pose_q = pose_.first;
  Eigen::Quaterniond q(pose_q(0), pose_q(1), pose_q(2), pose_q(3));
  Eigen::Matrix3d rotation_matrix = q.toRotationMatrix();
  Eigen::Vector3d translation = pose_.second;

  for (unsigned int i=0; i<pc_.points.size(); i++) {
    Eigen::Vector3d point(pc_.points[i].x, pc_.points[i].y, pc_.points[i].z);
    Eigen::Vector3d projected_point = rotation_matrix * point + translation;
    double scale = 1.0 / projected_point(2);
    projected_point *= scale;

    projected_cloud->points.push_back(pcl::PointXYZ(projected_point(0), projected_point(1), projected_point(2)));
  }


  // Normal estimation
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
  pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
  tree->setInputCloud (cloud);
  n.setInputCloud(cloud);
  n.setSearchMethod(tree);
  n.setKSearch(20);
  n.compute(*normals);

  // Store normals in normals_ as Eigen::Vector3d
  normals_.clear();
  for (unsigned int i = 0; i < normals->size(); i++) {
    Eigen::Vector3d normal;
    normal << normals->points[i].normal_x,
              normals->points[i].normal_y,
              normals->points[i].normal_z;
    normals_.push_back(normal);
  }


  // Recalculate with projected points, but don't save normals
  tree->setInputCloud (projected_cloud);
  n.setInputCloud(projected_cloud);
  n.setSearchMethod(tree);
  n.setKSearch(20);
  n.compute(*normals);





  // Concatenate the XYZ and normal fields*
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::PointNormal>);
  pcl::concatenateFields (*projected_cloud, *normals, *cloud_with_normals);

    // Create search tree*
  pcl::search::KdTree<pcl::PointNormal>::Ptr tree2 (new pcl::search::KdTree<pcl::PointNormal>);
  tree2->setInputCloud (cloud_with_normals);

  // Initialize objects
  pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;

  // Set the maximum distance between connected points (maximum edge length)
  gp3.setSearchRadius (1.5);

  // Set typical values for the parameters
  gp3.setMu (2.5);
  gp3.setMaximumNearestNeighbors (100);
  gp3.setMaximumSurfaceAngle(M_PI/4); // 45 degrees
  gp3.setMinimumAngle(M_PI/18); // 10 degrees
  gp3.setMaximumAngle(2*M_PI/3); // 120 degrees
  gp3.setNormalConsistency(false);

  // Get result
  gp3.setInputCloud (cloud_with_normals);
  gp3.setSearchMethod (tree2);
  gp3.reconstruct (mesh_);





  // Additional vertex information
  // To which part (cloud) does each vertex belong
  std::vector<int> parts = gp3.getPartIDs();
  // Whether the vertex status is [-1,0,1,2,3] = [NONE,FREE,FRINGE,BOUNDARY,COMPLETED]
  std::vector<int> states = gp3.getPointStates();

  // Get largest part
  std::unordered_map<int, int> map;
  int max_val = -1;
  int max_key = -1;
  for (unsigned int i=0; i<parts.size(); i++) {
    if (map.find(parts[i]) == map.end()) {
      map[parts[i]] = 1;
      if (max_val < 1) {
        max_val = 1;
        max_key = parts[i];
      }
    }
    else {
      map[parts[i]]++;
      if (map[parts[i]] > max_val) {
        max_val = map[parts[i]];
        max_key = parts[i];
      }
    }
  }

  // Get triangles of largest part only
  triangles_.clear();
  for (unsigned int i=0; i<mesh_.polygons.size(); i++) {
    unsigned int nver0 = mesh_.polygons[i].vertices[0];
    unsigned int nver1 = mesh_.polygons[i].vertices[1];
    unsigned int nver2 = mesh_.polygons[i].vertices[2];

    if (parts[nver0] == max_key && parts[nver1] == max_key && parts[nver2] == max_key) {
      std::vector<unsigned int> triangle;
      triangle.push_back(nver0);
      triangle.push_back(nver1);
      triangle.push_back(nver2);
      triangles_.push_back(triangle);
    }
  }

  // simulation: remove the last triangle
  if (simulation) {
    std::cout << "   Number of triangles: " << triangles_.size() << std::endl;
    for (unsigned int i=0; i<3; i++) {
      std::cout << " " << triangles_[triangles_.size() - 3][i];
    }
    std::cout << std::endl;
    for (unsigned int i=0; i<3; i++) {
      std::cout << " " << triangles_[triangles_.size() - 2][i];
    }
    std::cout << std::endl;
    for (unsigned int i=0; i<3; i++) {
      std::cout << " " << triangles_[triangles_.size() - 1][i];
    }
    std::cout << std::endl;
    triangles_.erase(triangles_.end() - 1);
    triangles_.erase(triangles_.end() - 1);
    triangles_.erase(triangles_.end() - 1);
    std::cout << "   Number of triangles: " << triangles_.size() << std::endl;
  }

  // Update points_alive_ and points_indices_
  std::vector<unsigned int> triangle_indices;
  for (auto t: triangles_) {
    for (auto i: t) {
      auto it = std::find(triangle_indices.begin(), triangle_indices.end(), i);
      if (it == triangle_indices.end()) {
        triangle_indices.push_back(i);
      }
    }
  }

  // sort ascending triangle_indices
  std::sort(triangle_indices.begin(), triangle_indices.end());

  //for (unsigned int i=0; i<triangle_indices.size(); i++) {
  //  std::cout << " " << triangle_indices[i];
  //}
  //std::cout << std::endl;



  std::vector<int> missing_indices, points_indices, points_indices_internal;
  int internal_idx = 0;
  for (unsigned int i=0; i<pc_.points.size(); i++) {
    auto it = std::find(triangle_indices.begin(), triangle_indices.end(), i);
    if (it == triangle_indices.end()) {
      missing_indices.push_back(i);
      points_alive_[points_indices_[i]] = false;
      points_indices_internal.push_back(-1);
    } else {
      points_indices.push_back(points_indices_[i]);
      points_indices_internal.push_back(internal_idx);
      internal_idx++;
    }
  }

  //std::cout << "Missing indices: ";
  //for (unsigned int i=0; i<missing_indices.size(); i++) {
  //  std::cout << " " << missing_indices[i];
  //}
  //std::cout << std::endl;

  //std::cout << "Points indices: ";
  //for (unsigned int i=0; i<points_indices.size(); i++) {
  //  std::cout << " " << points_indices[i];
  //}
  //std::cout << std::endl;

  //std::cout << "Triangulation: ";
  //for (unsigned int i=0; i<points_alive_.size(); i++) {
  //  std::cout << " " << points_alive_[i];
  //}
  //std::cout << std::endl;
  std::cout << "   Number of indices / points: " << points_indices.size() << " / " << points_.size() << std::endl;
  std::cout << "   Number of triangles: " << triangles_.size() << std::endl;

  points_indices_ = points_indices;

  // Replace pc_ with triangulated points
  pc_.width = points_indices_.size();
  pc_.height = 1;
  pc_.points.resize(pc_.width * pc_.height);

  for (size_t i = 0; i < pc_.points.size(); ++i) {
    pc_.points[i].x = points_[points_indices_[i]](0);
    pc_.points[i].y = points_[points_indices_[i]](1);
    pc_.points[i].z = points_[points_indices_[i]](2);
  }

  // Finally update triangles_ with internal indices
  for (unsigned int i=0; i<triangles_.size(); i++) {
    for (unsigned int j=0; j<triangles_[i].size(); j++) {
      triangles_[i][j] = points_indices_internal[triangles_[i][j]];
    }
  }

  //std::cout << "pc_.points.size() = " << pc_.points.size() << std::endl;


  
  return triangles_.size() > 0 ? true : false;
}



void FEM::SimulateFailedTriangulation() {
  std::cout << " - SimulateFailedTriangulation" << std::endl;
  std::cout << "   Original number of triangles: " << triangles_.size() << std::endl;
  // Remove triangle 0
  if (!triangles_.empty()) {
    triangles_.erase(triangles_.begin());
  }
  std::cout << "   Revised number of triangles: " << triangles_.size() << std::endl;

}


bool FEM::ClearNotTriangulated() {
  std::cout << " - ClearNotTriangulated" << std::endl;

  // Transform triangles accounting for points_active_
  //for (unsigned int i=0; i<triangles_.size(); i++) {
  //  for (unsigned int j=0; j<triangles_[i].size(); j++) {
  //    triangles_[i][j] = points_indices_[triangles_[i][j]];
  //  }
  //}

  std::vector<unsigned int> triangle_indices;
  for (auto t: triangles_) {
    for (auto i: t) {
      auto it = std::find(triangle_indices.begin(), triangle_indices.end(), i);
      if (it == triangle_indices.end()) {
        triangle_indices.push_back(i);
      }
    }
  }


  // Sum of vector of bools
  int points_alive_0 = std::accumulate(points_alive_.begin(), points_alive_.end(), 0);

  indices_not_triangulated_.clear();

  std::vector<int> points_indices;
  for (unsigned int i=0; i<points_indices_.size(); i++) {
    auto it = std::find(triangle_indices.begin(), triangle_indices.end(), i);
    if (it == triangle_indices.end()) {
      indices_not_triangulated_.push_back(i);
      points_alive_[points_indices_[i]] = false;
    } else {
      points_indices.push_back(points_indices_[i]);
    }
  }

  points_indices_ = points_indices;

  // Sum of vector of bools
  int points_alive_1 = std::accumulate(points_alive_.begin(), points_alive_.end(), 0);
  std::cout << "   Original point cloud size: " << points_alive_0 << std::endl;
  std::cout << "   Revised point cloud size: " << points_alive_1 << std::endl;

  if (points_alive_0 == points_alive_1) {
    std::cout << "   No points removed" << std::endl;
    return true;
  } 

  return InitCloud();
}

Eigen::Vector3d FEM::calculateMidpoint(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) {
    return (p1 + p2) / 2.0;
}

Eigen::Vector3d FEM::calculateOrthocenter(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) {
    return (p1 + p2 + p3) / 3.0;
}

bool FEM::TransformIntoQuads() {

    std::unordered_map<std::string, unsigned int> edgeMidpointIndexMap;
    
    for (const auto& triangle : triangles_) {
        std::vector<unsigned int> quadPoints;
        std::vector<unsigned int> midPointsIndices(3);

        // Calculate midpoints for each side and the orthocenter.
        for (int i = 0; i < 3; ++i) {
            unsigned int startIdx = triangle[i];
            unsigned int endIdx = triangle[(i + 1) % 3];
            std::string edgeKey = std::to_string(std::min(startIdx, endIdx)) + "-" + std::to_string(std::max(startIdx, endIdx));
            
            if (edgeMidpointIndexMap.find(edgeKey) == edgeMidpointIndexMap.end()) {
                Eigen::Vector3d midpoint = calculateMidpoint(points_[startIdx], points_[endIdx]);
                points_.push_back(midpoint);
                unsigned int newPointIndex = points_.size() - 1;
                edgeMidpointIndexMap[edgeKey] = newPointIndex;
                midPointsIndices[i] = newPointIndex;
            } else {
                midPointsIndices[i] = edgeMidpointIndexMap[edgeKey];
            }
        }

        Eigen::Vector3d orthocenter = calculateOrthocenter(points_[triangle[0]], points_[triangle[1]], points_[triangle[2]]);
        points_.push_back(orthocenter);
        unsigned int orthocenterIndex = points_.size() - 1;

        // For each side of the triangle, form a quadrilateral using one original vertex, two midpoints, and the orthocenter.
        for (int i = 0; i < 3; ++i) {
            std::vector<unsigned int> quadrilateral = {
                triangle[i],
                midPointsIndices[i],
                orthocenterIndex,
                midPointsIndices[(i + 2) % 3] // Use the previous side's midpoint by cycling backwards.
            };
            quadrilaterals_.push_back(quadrilateral);
        }
    }

  std::cout << " Number of Triangles: " << triangles_.size() << std::endl;
  std::cout << " Number of Quadrilaterals: " << quadrilaterals_.size() << std::endl;

  return true;
}






int FEM::CheckNodeOrderConsistency() {
  // Check if the triangles are ordered in a consistent way
  // If not, reorder them
  // Correct order is such that the normal of the triangle points outwards
  // and that the vertices are ordered in a counter-clockwise direction
  int nodes_reversed = 0;


  if (transform_into_quads_) {
    for (unsigned int i=0; i<quadrilaterals_.size(); i++) {
      pcl::PointXYZ p0 = pc_.points[quadrilaterals_[i][0]];
      pcl::PointXYZ p1 = pc_.points[quadrilaterals_[i][1]];
      pcl::PointXYZ p2 = pc_.points[quadrilaterals_[i][2]];
      pcl::PointXYZ p3 = pc_.points[quadrilaterals_[i][3]];

      Eigen::Vector3d v0(p0.x, p0.y, p0.z);
      Eigen::Vector3d v1(p1.x, p1.y, p1.z);
      Eigen::Vector3d v2(p2.x, p2.y, p2.z);
      Eigen::Vector3d v3(p3.x, p3.y, p3.z);

      Eigen::Vector3d normal1 = (v1 - v0).cross(v2 - v0);
      normal1.normalize();
      Eigen::Vector3d normal2 = (v2 - v0).cross(v3 - v0);
      normal2.normalize();

      if (normal1(2) < 0 || normal2(2) < 0) {
        unsigned int temp = quadrilaterals_[i][0];
        quadrilaterals_[i][0] = quadrilaterals_[i][1];
        quadrilaterals_[i][1] = quadrilaterals_[i][2];
        quadrilaterals_[i][2] = quadrilaterals_[i][3];
        quadrilaterals_[i][3] = temp;
        nodes_reversed++;
      }
    }

  } else {
    for (unsigned int i=0; i<triangles_.size(); i++) {
      pcl::PointXYZ p0 = pc_.points[triangles_[i][0]];
      pcl::PointXYZ p1 = pc_.points[triangles_[i][1]];
      pcl::PointXYZ p2 = pc_.points[triangles_[i][2]];

      Eigen::Vector3d v0(p0.x, p0.y, p0.z);
      Eigen::Vector3d v1(p1.x, p1.y, p1.z);
      Eigen::Vector3d v2(p2.x, p2.y, p2.z);

      Eigen::Vector3d normal = (v1 - v0).cross(v2 - v0);
      normal.normalize();

      if (normal(2) < 0) {
        unsigned int temp = triangles_[i][0];
        triangles_[i][0] = triangles_[i][1];
        triangles_[i][1] = temp;
        nodes_reversed++;
      }
    }

  }





  return nodes_reversed;
}



bool FEM::Compute(bool moving_least_squares, bool triangulate_planar, bool simulation) {

  std::cout << " - Compute " << element_ << " FEM model" << std::endl;

  bool ok = InitCloud();
  std::cout << "   Result: " << ok << std::endl;

  if (element_ == "C3D8") {
    transform_into_quads_ = true;
  } else {
    transform_into_quads_ = false;
  }

  if (ok && moving_least_squares) {
    std::cout << " - Moving Least Squares:" << std::endl;
    ok = MovingLeastSquares(simulation);
    std::cout << "   Result: " << ok << std::endl;
  }

  if (ok) {
    std::cout << " - Triangulation: " << ok << std::endl;
    if (triangulate_planar) {
      ok = TriangulateByProjection(simulation);
    } else {
      ok = Triangulate(simulation);
    }
    std::cout << "   Result: " << ok << std::endl;
  }

  if (ok && transform_into_quads_) {
    std::cout << " - TransformIntoQuads:" << std::endl;
    ok = TransformIntoQuads();
    std::cout << "   Result: " << ok << std::endl;
    ok = InitCloud();
    std::cout << "   Cloud Re-Initialization: " << ok << std::endl;
  }

  if (ok) {
    std::cout << " - CheckNodeOrderConsistency:" << std::endl;
    int nodes_reversed = CheckNodeOrderConsistency();
    std::cout << "   " << nodes_reversed << " nodes reversed" << std::endl;
    ok = nodes_reversed != -1;
  }
  

  if (ok and element_ != "C2D4") {
    std::cout << " - ComputeExtrusion:" << std::endl;
    ok = ComputeExtrusion();
  }



  return ok;
}

bool FEM::ComputeExtrusion() {

  //std::cout << "EXTRUSION" << std::endl;
  //std::cout << " - Number of triangles: " << triangles_.size() << std::endl;

  std::vector<double> distances;
  for (std::vector<unsigned int> triangle : triangles_) {
    pcl::PointXYZ p0 = pc_.points[triangle[0]];
    pcl::PointXYZ p1 = pc_.points[triangle[1]];
    pcl::PointXYZ p2 = pc_.points[triangle[2]];

    double d01 = pcl::euclideanDistance(p0, p1);
    double d12 = pcl::euclideanDistance(p1, p2);
    double d20 = pcl::euclideanDistance(p2, p0);

    distances.push_back(d01);
    distances.push_back(d12);
    distances.push_back(d20);
  }

  // Get median
  std::sort(distances.begin(), distances.end());
  element_height_ = distances[distances.size()/2];
  std::cout << "   Computed element height: mag = " << element_height_ << std::endl;
  
  // Compute second layer at a distance of 1/2 element height and with the direction of the normal vector
  points2_.clear();
  pc2_.width = pc_.width;
  pc2_.height = pc_.height;
  pc2_.points.resize(pc2_.width * pc2_.height);

  // Average normal
  Eigen::Vector3d normal = Eigen::Vector3d::Zero();
  for (Eigen::Vector3d n : normals_) {
    normal -= n;
  }
  normal /= normals_.size();    

  // truncate normals to 3 decimals
  double scale = std::pow(10.0, 3);
  for (unsigned int i=0; i<3; i++) {
    normal(i) = std::round(normal(i) * scale) / scale;
  }
  

  int cloud_idx = 0;
  for (unsigned int i=0; i<points_.size(); i++) {
    Eigen::Vector3d point2 = points_[i] - element_height_/2 * normal;
    points2_.push_back(point2);
    if (points_alive_[i]) {
      pc2_.points[cloud_idx].x = point2(0);
      pc2_.points[cloud_idx].y = point2(1);
      pc2_.points[cloud_idx].z = point2(2);
      cloud_idx++;
    }
  }

  int points_alive = std::accumulate(points_alive_.begin(), points_alive_.end(), 0);
  //std::cout << " - points_alive: " << points_alive << std::endl;
  //std::cout << "pc_size: " << pc_.points.size() << std::endl;
  //std::cout << "pc2_size: " << pc2_.points.size() << std::endl;
  //for (unsigned int i=0; i<pc2_.points.size(); i++) {
  //  std::cout << " " << i << " " << pc_.points[i].x << " " << pc_.points[i].y << " " << pc_.points[i].z << " -> "
  //            << " " << pc2_.points[i].x << " " << pc2_.points[i].y << " " << pc2_.points[i].z << std::endl;
  //}

  for (std::vector<unsigned int> triangle : triangles_) {
    std::vector<unsigned int> element;
    element.push_back(triangle[0] + points_alive);
    element.push_back(triangle[1] + points_alive);
    element.push_back(triangle[2] + points_alive);
    element.push_back(triangle[0]);
    element.push_back(triangle[1]);
    element.push_back(triangle[2]);
    elements_.push_back(element);
  }

  return true;
}

std::vector<Eigen::Vector3d> FEM::GetExtrusionDelta() {
  std::vector<Eigen::Vector3d> delta;
  for (unsigned int i=0; i<points_.size(); i++) {
    delta.push_back(points2_[i] - points_[i]);
  }
  return delta;
}

std::vector<Eigen::Vector3d> FEM::GetExtrusion() {
  return points2_;
}

std::vector<unsigned int> FEM::GetExtrusionIndices() {
  int points_alive = std::accumulate(points_alive_.begin(), points_alive_.end(), 0);
  std::vector<unsigned int> indices;
  for (unsigned int i=0; i<points_alive; i++) {
    indices.push_back(i + points_alive);
  }
  //for (unsigned int i=0; i<points_.size(); i++) {
  //  if (points_alive_[i]) {
  //    int idx = indices.size() + points_alive;
  //    indices.push_back(idx);
  //  }
  //}
  return indices;
}

void FEM::SetExtrusion(std::vector<Eigen::Vector3d> extrusion_delta, double element_height) {
  pc2_.width = pc_.width;
  pc2_.height = pc_.height;
  pc2_.points.resize(pc2_.width * pc2_.height);
  
  int cloud_idx = 0;
  for (unsigned int i=0; i<points_.size(); i++) {
    Eigen::Vector3d pt2 = points_[i] + extrusion_delta[i];
    points2_.push_back(pt2);
    if (points_alive_[i]) {
      pc2_.points[cloud_idx].x = pt2(0);
      pc2_.points[cloud_idx].y = pt2(1);
      pc2_.points[cloud_idx].z = pt2(2);
      cloud_idx++;
    }
  }
  element_height_ = element_height;



}

void FEM::SetAliveNodes(std::vector<bool> alive_nodes) {
  points_alive_ = alive_nodes;
}

double FEM::GetElementHeight() {
  return element_height_;
}

std::vector<unsigned int> FEM::GetIndicesNotTriangulated() {
  return indices_not_triangulated_;
}

std::vector<Eigen::Vector3d> FEM::GetPoints(bool alive_only) {
  if (!alive_only) {
    return points_;
  } else {
    std::vector<Eigen::Vector3d> points_alive;
    for (unsigned int i=0; i<points_.size(); i++) {
      if (points_alive_[i]) {
        points_alive.push_back(points_[i]);
      }
    }
    return points_alive;
  }
}













void FEM::ViewMesh(bool extrusion,
                   int wait,
                   pcl::visualization::PCLVisualizer viewer) {

  std::cout << "VIEW MESH" << std::endl;

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudf(new pcl::PointCloud<pcl::PointXYZ> (pc_));
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudb(new pcl::PointCloud<pcl::PointXYZ> (pc2_));

  //pcl::visualization::PCLVisualizer viewer;
  viewer.setBackgroundColor(1, 1, 1);
  viewer.addPointCloud(cloudf);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 0.0);
  
  int points_alive = std::accumulate(points_alive_.begin(), points_alive_.end(), 0);
  std::cout << "Number of triangles: " << triangles_.size() << std::endl;
  std::cout << "Number of points: " << cloudf->points.size() << std::endl;
  std::cout << "Number of points alive: " << points_alive << std::endl;

  for (unsigned int i=0; i<triangles_.size(); i++) {

    //std::cout << "Triangle: " << i << ": ";
    //for (unsigned int j=0; j<triangles_[i].size(); j++) {
    //  std::cout << " " << triangles_[i][j];
    //}
    //std::cout << std::endl;

    pcl::PointXYZ p0 = cloudf->points[triangles_[i][0]];
    pcl::PointXYZ p1 = cloudf->points[triangles_[i][1]];
    pcl::PointXYZ p2 = cloudf->points[triangles_[i][2]];
    std::string name = "triangle_1_" + std::to_string(i);

    // Add line
    viewer.addLine<pcl::PointXYZ>(p0, p1, 0.0, 0.7, 0.0, name + "a");
    viewer.addLine<pcl::PointXYZ>(p1, p2, 0.0, 0.7, 0.0, name + "b");
    viewer.addLine<pcl::PointXYZ>(p2, p0, 0.0, 0.7, 0.0, name + "c");
  }


  if (extrusion) {
    for (unsigned int i=0; i<triangles_.size(); i++) {
      pcl::PointXYZ p0 = cloudb->points[triangles_[i][0]];
      pcl::PointXYZ p1 = cloudb->points[triangles_[i][1]];
      pcl::PointXYZ p2 = cloudb->points[triangles_[i][2]];
      std::string name = "triangle_1_" + std::to_string(i) + "_extrusion";

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.9, 0.7, name + "_a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, 0.7, 0.9, 0.7, name + "_b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, 0.7, 0.9, 0.7, name + "_c");
    }

    for (unsigned int i=0; i<cloudf->points.size(); i++) {
      pcl::PointXYZ p0 = cloudf->points[i];
      pcl::PointXYZ p1 = cloudb->points[i];
      std::string name = "line_1_" + std::to_string(i);
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.9, 0.7, name);
    }
  }

  // If pose_ not zero, then add as a thick point
  if (pose_.second.norm() > 0) {
    pcl::PointXYZ p1(pose_.second(0), pose_.second(1), pose_.second(2));
    viewer.addSphere(p1, 0.1, 0.0, 0.7, 0.0, "pose_");

    // Add a line in the direction of the quaternion in pose2.first
    std::pair<Eigen::Vector3d, Eigen::Vector3d> line_pts = QuaternionLine(pose_.first, pose_.second);
    pcl::PointXYZ p1_a(line_pts.first(0), line_pts.first(1), line_pts.first(2));
    pcl::PointXYZ p1_b(line_pts.second(0), line_pts.second(1), line_pts.second(2));
    viewer.addLine<pcl::PointXYZ>(p1_a, p1_b, 0.0, 0.7, 0.0, "pose_line");
  }

  viewer.addCoordinateSystem(0.5);

  if (wait > 0) {
    viewer.spinOnce(wait*1000);
    return;
  } else {
    viewer.spin();
  }
  viewer.close();
}





void FEM::ViewMesh(bool extrusion, 
                   std::vector<Eigen::Vector3d> points2,
                   std::vector<Eigen::Vector3d> points2extrusion,
                   std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2,
                   int wait,
                   pcl::visualization::PCLVisualizer viewer) {

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudf(new pcl::PointCloud<pcl::PointXYZ> (pc_));
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudb(new pcl::PointCloud<pcl::PointXYZ> (pc2_));

  //pcl::visualization::PCLVisualizer viewer;
  viewer.setBackgroundColor(1, 1, 1);
  viewer.addPointCloud(cloudf);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 0.0);
  


  // Draw first model
  
  for (unsigned int i=0; i<triangles_.size(); i++) {
    pcl::PointXYZ p0 = cloudf->points[triangles_[i][0]];
    pcl::PointXYZ p1 = cloudf->points[triangles_[i][1]];
    pcl::PointXYZ p2 = cloudf->points[triangles_[i][2]];
    std::string name = "triangle_1_" + std::to_string(i);

    // Add line
    viewer.addLine<pcl::PointXYZ>(p0, p1, 0.0, 0.7, 0.0, name + "a");
    viewer.addLine<pcl::PointXYZ>(p1, p2, 0.0, 0.7, 0.0, name + "b");
    viewer.addLine<pcl::PointXYZ>(p2, p0, 0.0, 0.7, 0.0, name + "c");
  }

  for (unsigned int i=0; i<points_alive_.size(); i++) {
    if (points_alive_[i]) {
      continue;
    }
    viewer.addSphere(cloudf->points[i], 0.01, 0.0, 0.7, 0.0, "point_not_alive_" + std::to_string(i));
  }

  if (extrusion) {
    for (unsigned int i=0; i<triangles_.size(); i++) {
      pcl::PointXYZ p0 = cloudb->points[triangles_[i][0]];
      pcl::PointXYZ p1 = cloudb->points[triangles_[i][1]];
      pcl::PointXYZ p2 = cloudb->points[triangles_[i][2]];
      std::string name = "triangle_1_" + std::to_string(i) + "_extrusion";

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.9, 0.7, name + "_a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, 0.7, 0.9, 0.7, name + "_b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, 0.7, 0.9, 0.7, name + "_c");
    }

    for (unsigned int i=0; i<cloudf->points.size(); i++) {
      if (!points_alive_[i]) {
        continue;
      }

      pcl::PointXYZ p0 = cloudf->points[i];
      pcl::PointXYZ p1 = cloudb->points[i];
      std::string name = "line_1_" + std::to_string(i);
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.9, 0.7, name);
    }
  }

  // If pose_ not zero, then add as a thick point
  if (pose_.second.norm() > 0) {
    pcl::PointXYZ p1(pose_.second(0), pose_.second(1), pose_.second(2));
    viewer.addSphere(p1, 0.1, 0.0, 0.7, 0.0, "pose_");

    // Add a line in the direction of the quaternion in pose2.first
    std::pair<Eigen::Vector3d, Eigen::Vector3d> line_pts = QuaternionLine(pose_.first, pose_.second);
    pcl::PointXYZ p1_a(line_pts.first(0), line_pts.first(1), line_pts.first(2));
    pcl::PointXYZ p1_b(line_pts.second(0), line_pts.second(1), line_pts.second(2));
    viewer.addLine<pcl::PointXYZ>(p1_a, p1_b, 0.0, 0.7, 0.0, "pose_line");
  }




  // Draw second model


  pcl::PointCloud<pcl::PointXYZ> cloud2f, cloud2b;
  cloud2f.width = cloud2b.width = pc2_.width;
  cloud2f.height = cloud2b.height = pc2_.height;
  cloud2f.points.resize(cloud2f.width * cloud2f.height);
  cloud2b.points.resize(cloud2b.width * cloud2b.height);

  int idx = 0;
  for (unsigned int i=0; i<points_.size(); i++) {
    if (points_alive_[i]) {
      pcl::PointXYZ p(points2[idx](0), points2[idx](1), points2[idx](2));
      cloud2f.points[i] = p;

      pcl::PointXYZ p2(points2extrusion[idx](0), points2extrusion[idx](1), points2extrusion[idx](2));
      cloud2b.points[i] = p2;

      idx++;
    }
  }



  for (unsigned int i=0; i<triangles_.size(); i++) {
    pcl::PointXYZ p0 = cloud2f.points[triangles_[i][0]];
    pcl::PointXYZ p1 = cloud2f.points[triangles_[i][1]];
    pcl::PointXYZ p2 = cloud2f.points[triangles_[i][2]];
    std::string name = "triangle_2_" + std::to_string(i);

    // Add line
    viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.0, 0.0, name + "a");
    viewer.addLine<pcl::PointXYZ>(p1, p2, 0.7, 0.0, 0.0, name + "b");
    viewer.addLine<pcl::PointXYZ>(p2, p0, 0.7, 0.0, 0.0, name + "c");
  }


  if (extrusion) {
    for (unsigned int i=0; i<triangles_.size(); i++) {
      pcl::PointXYZ p0 = cloud2b.points[triangles_[i][0]];
      pcl::PointXYZ p1 = cloud2b.points[triangles_[i][1]];
      pcl::PointXYZ p2 = cloud2b.points[triangles_[i][2]];
      std::string name = "triangle_2_" + std::to_string(i) + "_extrusion";

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.9, 0.7, 0.7, name + "_a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, 0.9, 0.7, 0.7, name + "_b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, 0.9, 0.7, 0.7, name + "_c");
    }

    for (unsigned int i=0; i<cloud2f.points.size(); i++) {
      if (!points_alive_[i]) {
        continue;
      }

      pcl::PointXYZ p0 = cloud2f.points[i];
      pcl::PointXYZ p1 = cloud2b.points[i];
      
      std::string name = "line_2_" + std::to_string(i);

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.9, 0.7, 0.7, name);
    }

  }

  // If pose2 not zero, then add as a thick point
  if (pose2.second.norm() > 0) {
    pcl::PointXYZ p2(pose2.second(0), pose2.second(1), pose2.second(2));
    viewer.addSphere(p2, 0.1, 0.7, 0.0, 0.0, "pose2");
    // Add a line in the direction of the quaternion in pose2.first
    std::pair<Eigen::Vector3d, Eigen::Vector3d> line_pts = QuaternionLine(pose2.first, pose2.second);
    pcl::PointXYZ p2_a(line_pts.first(0), line_pts.first(1), line_pts.first(2));
    pcl::PointXYZ p2_b(line_pts.second(0), line_pts.second(1), line_pts.second(2));
    viewer.addLine<pcl::PointXYZ>(p2_a, p2_b, 0.7, 0.0, 0.0, "pose2_line");
  }


  viewer.addCoordinateSystem(0.5);
  

  if (wait > 0) {
    viewer.spinOnce(wait*1000);
    return;
  } else {
    viewer.spin();
  }
  viewer.close();
}




std::vector<std::vector<float>> FEM::GetNodes() {
  std::vector<std::vector<float>> points;
  for (Eigen::Vector3d pt: points_) {
    std::vector<float> point;
    point.push_back(pt[0]);
    point.push_back(pt[1]);
    point.push_back(pt[2]);
    points.push_back(point);
  }
  for (Eigen::Vector3d pt: points2_) {
    std::vector<float> point;
    point.push_back(pt[0]);
    point.push_back(pt[1]);
    point.push_back(pt[2]);
    points.push_back(point);
  }
  return points;
}


std::vector<Eigen::Vector3d> FEM::GetEigenNodesFront(bool active_only) {
  std::vector<Eigen::Vector3d> points;
  for (unsigned int i=0; i<points_.size(); i++) {
    if (active_only && !points_alive_[i]) {
      continue;
    }
    points.push_back(points_[i]);
  }
  return points;
}


std::vector<Eigen::Vector3d> FEM::GetEigenNodes(bool active_only) {
  std::vector<Eigen::Vector3d> points;
  for (unsigned int i=0; i<points_.size(); i++) {
    if (active_only && !points_alive_[i]) {
      continue;
    }
    points.push_back(points_[i]);
  }
  for (unsigned int i=0; i<points2_.size(); i++) {
    if (active_only && !points_alive_[i] && points_.size() == points2_.size()) {
      continue;
    }
    points.push_back(points2_[i]);
  }
  return points;
}


//std::vector<Eigen::Vector3d> FEM::GetEigenNodes(bool active_only, bool is_target) {
//  std::vector<Eigen::Vector3d> points;
//  for (unsigned int i=0; i<points_.size(); i++) {
//    if (active_only && !points_alive_[i]) {
//      continue;
//    }
//    points.push_back(points_[i]);
//  }
//
//  int points_alive = std::accumulate(points_alive_.begin(), points_alive_.end(), 0);
//
//  if (!is_target) {
//    for (unsigned int i=0; i<points2_.size(); i++) {
//      if (active_only && !points_alive_[i]) {
//        continue;
//      }
//      points.push_back(points2_[i]);
//    }
//  }
//  //for (unsigned int i=0; i<points2_.size(); i++) {
//  //  if (is_target && active_only && !points_alive_[i]) {
//  //    continue;
//  //  }
//  //  points.push_back(points2_[i]);
//  //}
//
//  std::cout << " - points alive / l1 / l2 / extracted / available: " << points_alive 
//            << " / " << points_.size() << " / " << points2_.size() 
//            << " / " << points.size() << " / " << points_.size() + points2_.size() 
//            << std::endl;
//
//  return points;
//}

std::vector<Eigen::Vector3d> FEM::GetEigenBaseNodes() {
  return points_;
}

std::vector<Eigen::Vector3d> FEM::GetEigenExtrudedNodes() {
  return points2_;
}

std::vector<bool> FEM::GetPointsAlive() {
  return points_alive_;
}

std::vector<std::vector<unsigned int>> FEM::GetTriangles() {
  return triangles_;
}

void FEM::SetTriangles(std::vector<std::vector<unsigned int>> triangles) {
  triangles_ = triangles;
}

std::vector<std::vector<unsigned int>> FEM::GetElements(bool alive_only) {
  if (!alive_only) {
    return elements_;
  }

  //std::cout << "points_indices (" << points_indices_.size() << "):";
  //for (unsigned int i=0; i<points_indices_.size(); i++) {
  //  std::cout << " " << points_indices_[i];
  //}
  //std::cout << std::endl;
  
  int points_alive = std::accumulate(points_alive_.begin(), points_alive_.end(), 0);

  std::vector<std::vector<unsigned int>> elements;
  for (std::vector<unsigned int> triangle : triangles_) {
    std::vector<unsigned int> element;
    unsigned int triangle_0 = points_indices_[triangle[0]];
    unsigned int triangle_1 = points_indices_[triangle[1]];
    unsigned int triangle_2 = points_indices_[triangle[2]];

    //std::cout << " " << triangle_0 << "-" << triangle[0] << "-" <<
    //                    triangle_1 << "-" << triangle[1] << "-" <<
    //                    triangle_2 << "-" << triangle[2] << std::endl;

    element.push_back(triangle[0] + points_alive);
    element.push_back(triangle[1] + points_alive);
    element.push_back(triangle[2] + points_alive);
    element.push_back(triangle[0]);
    element.push_back(triangle[1]);
    element.push_back(triangle[2]);
    
    elements.push_back(element);
  }

  return elements;
}


void FEM::SetElements(std::vector<std::vector<unsigned int>> elements) {
  elements_ = elements;
}

pcl::PointCloud<pcl::PointXYZ> FEM::GetCloud() {
  return pc_;
}

pcl::PointCloud<pcl::PointXYZ> FEM::GetCloud2() {
  return pc2_;
}

std::pair<Eigen::Vector4d, Eigen::Vector3d> FEM::GetPose() {
  return pose_;
}

void FEM::SetPose(std::pair<Eigen::Vector4d, Eigen::Vector3d> pose) {
  pose_ = pose;
}


void FEM::Transform(std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes,
                    std::pair<Eigen::Vector4d, Eigen::Vector3d> pose) {

  //std::cout << "Size pc_/points_/nodes.first :    "  << pc_.points.size() << " / " << points_.size() << " / " << nodes.first.size() << std::endl;
  //std::cout << "Size pc2_/points2_/nodes.second : "  << pc2_.points.size() << " / " << points2_.size() << " / " << nodes.second.size() << std::endl;

  if (nodes.first.size() == points_.size()) {
    for (unsigned int i=0; i<points_.size(); i++) {
      Eigen::Vector3d pt1 = nodes.first[i];
      points_[i] = pt1;
    }
  }
  if (nodes.second.size() == points2_.size()) {
    for (unsigned int i=0; i<points2_.size(); i++) {
      Eigen::Vector3d pt2 = nodes.second[i];
      points2_[i] = pt2;
    }
  }


//  int node_idx = 0;
//  for (unsigned int i=0; i<points_.size(); i++) {
//    std::cout << "i: " << i << " / " << points_.size() << std::endl;
//
//    if (!points_alive_[i]) {
//      continue;
//    }
//
//    Eigen::Vector3d pt1 = Eigen::Vector3d(nodes.first[node_idx][0], nodes.first[node_idx][1], nodes.first[node_idx][2]);
//    points_[i] = pt1;
//    pc_.points[i].x = points_[i](0);
//    pc_.points[i].y = points_[i](1);
//    pc_.points[i].z = points_[i](2);
//
//    //Eigen::Vector3d pt2 = Eigen::Vector3d(nodes.second[node_idx][0], nodes.second[node_idx][1], nodes.second[node_idx][2]);
//    //points2_[i] = pt2;
//    //pc2_.points[i].x = points2_[i](0);
//    //pc2_.points[i].y = points2_[i](1);
//    //pc2_.points[i].z = points2_[i](2);
//
//    node_idx++;
//  }

  

  for (unsigned int i=0; i<nodes.first.size(); i++) {
    Eigen::Vector3d pt2 = Eigen::Vector3d(nodes.first[i][0], nodes.first[i][1], nodes.first[i][2]);
    pc_.points[i].x = nodes.first[i][0];
    pc_.points[i].y = nodes.first[i][1];
    pc_.points[i].z = nodes.first[i][2];
  }

  for (unsigned int i=0; i<points2_.size(); i++) {
    Eigen::Vector3d pt2 = Eigen::Vector3d(nodes.second[i][0], nodes.second[i][1], nodes.second[i][2]);
    points2_[i] = pt2;
    pc2_.points[i].x = points2_[i](0);
    pc2_.points[i].y = points2_[i](1);
    pc2_.points[i].z = points2_[i](2);
  }



  //pose_ = pose;  
  pose_.first = Eigen::Vector4d(pose.first(0), pose.first(1), pose.first(2), pose.first(3));
  pose_.second = Eigen::Vector3d(pose.second(0), pose.second(1), pose.second(2));
}


std::pair<Eigen::Vector3d, Eigen::Vector3d> FEM::QuaternionLine2(
    Eigen::Vector4d qvec, Eigen::Vector3d point, double radius) {

  // Step 1: Normalize the quaternion
  Eigen::Quaterniond normalizedQuaternion(qvec.normalized());

  // Step 2: Convert the quaternion to a rotation matrix
  Eigen::Matrix3d rotationMatrix = normalizedQuaternion.toRotationMatrix();

  // Step 3: Define the line endpoint
  Eigen::Vector3d direction(radius, 0.0, 0.0); // Direction of the line
  Eigen::Vector3d endpoint = point + rotationMatrix * direction;

  return std::make_pair(point, endpoint);
}

std::pair<Eigen::Vector3d, Eigen::Vector3d> FEM::QuaternionLine(
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
                   std::vector<FEM*> fems) {

  pcl::visualization::PCLVisualizer viewer;
  viewer.setBackgroundColor(1, 1, 1);

  std::vector<std::vector<double>> colorsf, colorsb;
  colorsf.push_back({0.0, 0.7, 0.0});
  colorsf.push_back({0.7, 0.0, 0.0});
  colorsf.push_back({0.0, 0.0, 0.7});
  colorsb.push_back({0.7, 0.9, 0.7});
  colorsb.push_back({0.9, 0.7, 0.7});
  colorsb.push_back({0.7, 0.7, 0.9});

  for (unsigned int f=0; f<fems.size(); f++) {

    FEM* fem = fems[f];

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudf(new pcl::PointCloud<pcl::PointXYZ> (fem->GetCloud()));
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudb(new pcl::PointCloud<pcl::PointXYZ> (fem->GetCloud2()));
    std::vector<std::vector<unsigned int>> triangles = fem->GetTriangles();
    std::pair<Eigen::Vector4d, Eigen::Vector3d> pose = fem->GetPose();

    // add point cloud
    //viewer.addPointCloud(cloudf);
    //viewer.addPointCloud(cloudb);
    //viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2);
    //viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 0.0);

    
    for (unsigned int i=0; i<triangles.size(); i++) {
      pcl::PointXYZ p0 = cloudf->points[triangles[i][0]];
      pcl::PointXYZ p1 = cloudf->points[triangles[i][1]];
      pcl::PointXYZ p2 = cloudf->points[triangles[i][2]];

      std::string name = "triangle_" + std::to_string(f) + "_" + std::to_string(i);

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, colorsf[f][0], colorsf[f][1], colorsf[f][2], name + "_a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, colorsf[f][0], colorsf[f][1], colorsf[f][2], name + "_b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, colorsf[f][0], colorsf[f][1], colorsf[f][2], name + "_c");
    }


    if (extrusion) {
      for (unsigned int i=0; i<triangles.size(); i++) {
        pcl::PointXYZ p0 = cloudb->points[triangles[i][0]];
        pcl::PointXYZ p1 = cloudb->points[triangles[i][1]];
        pcl::PointXYZ p2 = cloudb->points[triangles[i][2]];
        std::string name = "triangle_" + std::to_string(f) + "_" + std::to_string(i) + "_extrusion";

        // Add line
        viewer.addLine<pcl::PointXYZ>(p0, p1, colorsb[f][0], colorsb[f][1], colorsb[f][2], name + "_a");
        viewer.addLine<pcl::PointXYZ>(p1, p2, colorsb[f][0], colorsb[f][1], colorsb[f][2], name + "_b");
        viewer.addLine<pcl::PointXYZ>(p2, p0, colorsb[f][0], colorsb[f][1], colorsb[f][2], name + "_c");
      }

      for (unsigned int i=0; i<cloudf->points.size(); i++) {
        pcl::PointXYZ p0 = cloudf->points[i];
        pcl::PointXYZ p1 = cloudb->points[i];
        std::string name = "line_" + std::to_string(f) + "_" + std::to_string(i);
        viewer.addLine<pcl::PointXYZ>(p0, p1, colorsb[f][0], colorsb[f][1], colorsb[f][2], name);
      }
    }

    // If pose_ not zero, then add as a thick point
    if (pose.second.norm() > 0) {
      pcl::PointXYZ p1(pose.second(0), pose.second(1), pose.second(2));
      viewer.addSphere(p1, 0.1, colorsf[f][0], colorsf[f][1], colorsf[f][2], "pose_" + std::to_string(f));

      // Add a line in the direction of the quaternion in pose2.first
      std::pair<Eigen::Vector3d, Eigen::Vector3d> line_pts = fem->QuaternionLine(pose.first, pose.second);
      pcl::PointXYZ p1_a(line_pts.first(0), line_pts.first(1), line_pts.first(2));
      pcl::PointXYZ p1_b(line_pts.second(0), line_pts.second(1), line_pts.second(2));
      viewer.addLine<pcl::PointXYZ>(p1_a, p1_b, colorsf[f][0], colorsf[f][1], colorsf[f][2], "pose_line_" + std::to_string(f));
    }
  }

  viewer.addCoordinateSystem(0.5);
  if (wait > 0) {
    viewer.spinOnce(wait*1000);
    return;
  } else {
    viewer.spin();
  }
  viewer.close();
}