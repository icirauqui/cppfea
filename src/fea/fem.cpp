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


FEM::FEM(std::string element): element_(element) {}


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

  if (pc_.width == 0) {
    std::cout << " - InitCloud failed" << std::endl;
  } else {
    std::cout << " - InitCloud success" << std::endl;
    pose_ = ApproximatePose(points_);
    std::cout << " - Pose: " << pose_.first.transpose() << " | " << pose_.second.transpose() << std::endl;
  }

  return pc_.width > 0;
}

bool FEM::MovingLeastSquares() {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
    new pcl::PointCloud<pcl::PointXYZ> (pc_));

  // Create a KD-Tree
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);

  // Output has the PointNormal type in order to store the normals calculated by MLS
  pcl::PointCloud<pcl::PointNormal> mls_points;

  // Init object (second point type is for the normals, even if unused)
  pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
 
  mls.setComputeNormals (true);

  // Set parameters
  mls.setInputCloud (cloud);
  mls.setPolynomialOrder (2);
  mls.setSearchMethod (tree);
  mls.setSearchRadius (2);

  // Reconstruct
  mls.process (mls_points);

  if (mls_points.size() == 0) {
    std::cout << "MLS failed" << std::endl;
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

  // simulation: erase position 6 in mls_indices_
  //mls_indices_.erase(mls_indices_.begin() + 6);
  
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

  // simulation
  //for (int i = 0; i<points_indices_.size(); i++) {
  //  std::cout << " " << points_indices_[i];
  //}
  //std::cout << std::endl;

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


bool FEM::Triangulate() {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
    new pcl::PointCloud<pcl::PointXYZ> (pc_));

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

    // Transform triangles accounting for points_active_
    //for (unsigned int i=0; i<triangles_.size(); i++) {
    //  for (unsigned int j=0; j<triangles_[i].size(); j++) {
    //    triangles_[i][j] = points_indices_[triangles_[i][j]];
    //  }
    //}

  }
  
  return triangles_.size() > 0 ? true : false;
}


/*
void FEM::SimulateFailedTriangulation() {
  // 5% of points are removed
  int num_points = points_.size();
  int num_removed = num_points * 0.05;
  if (num_removed < 1) {
    num_removed = 1;
  }

  std::vector<int> removed_indices;
  //for (int i = 0; i < num_removed; i++) {
  //  int idx = rand() % num_points;
  //  points_alive_[idx] = false;
  //  removed_indices.push_back(idx);
  //}
  removed_indices.push_back(0);
  points_alive_[0] = false;
  removed_indices.push_back(2);
  points_alive_[2] = false;

  std::cout << " - SimulateFailedTriangulation: " << num_removed << " points removed" << std::endl;
  std::cout << "   Indices:";
  for (auto i: removed_indices) {
    std::cout << " " << i;
  }
  std::cout << std::endl;

  // Copy triangles vector
  std::vector<std::vector<unsigned int>> old_triangles;
  for (auto t: triangles_) {
    std::vector<unsigned int> newt;
    for (auto i: t) {
      newt.push_back(i);
    }
    old_triangles.push_back(newt);
  }
  triangles_.clear();

  for (auto t: old_triangles) {
    std::vector<unsigned int> newt;
    for (auto i: t) {
      if (points_alive_[i]) {
        newt.push_back(i);
      }
    }
    if (newt.size() == 3) {
      triangles_.push_back(newt);
    }
  }

  std::cout << "   Number of original triangles: " << old_triangles.size() << std::endl;
  std::cout << "   Number of revised triangles:  " << triangles_.size() << std::endl;
}
*/

void FEM::SimulateFailedTriangulation() {
  std::cout << " - SimulateFailedTriangulation" << std::endl;
  std::cout << "   Original number of triangles: " << triangles_.size() << std::endl;
  // Remove triangle 0
  if (!triangles_.empty()) {
    triangles_.erase(triangles_.begin());
  }
  std::cout << "   Revised number of triangles: " << triangles_.size() << std::endl;

}

/*
bool FEM::ClearNotTriangulated() {
  std::cout << " - ClearNotTriangulated" << std::endl;

  std::vector<unsigned int> point_map = std::vector<unsigned int>(points_.size(), -1);
  std::vector<unsigned int> triangle_indices;
  for (auto t: triangles_) {
    for (auto i: t) {
      auto it = std::find(triangle_indices.begin(), triangle_indices.end(), i);
      if (it == triangle_indices.end()) {
        triangle_indices.push_back(i);
      }
    }
  }

  //std::cout << "   ";
  //for (unsigned int i=0; i<triangle_indices.size(); i++) {
  //  std::cout << " " << triangle_indices[i];
  //}
  //std::cout << std::endl;

  std::vector<Eigen::Vector3d> new_points;
  std::vector<bool> new_points_alive;
  indices_not_triangulated_.clear();

  for (unsigned int i=0; i<points_.size(); i++) {
    auto it = std::find(triangle_indices.begin(), triangle_indices.end(), i);
    if (it != triangle_indices.end()) {
      new_points.push_back(points_[i]);
      new_points_alive.push_back(true);
      point_map[i] = new_points.size() - 1;
    }
    else {
      indices_not_triangulated_.push_back(i);
    }
  }

  std::cout << "   Original point cloud size: " << points_.size() << std::endl;
  std::cout << "   Revised point cloud size: " << new_points.size() << std::endl;

  if (points_.size() == new_points.size()) {
    std::cout << "   No points removed" << std::endl;
    return true;
  } 

  points_ = new_points;
  points_alive_ = new_points_alive;

  // Modify triangle indices
  for (unsigned int i=0; i<triangles_.size(); i++) {
    for (unsigned int j=0; j<triangles_[i].size(); j++) {
      triangles_[i][j] = point_map[triangles_[i][j]];
    }
  }

  return InitCloud();
}
*/

bool FEM::ClearNotTriangulated() {
  std::cout << " - ClearNotTriangulated" << std::endl;

  //// Transform triangles accounting for points_active_
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



bool FEM::Compute(bool moving_least_squares, bool simulation) {
  bool ok = InitCloud();

  if (element_ == "C3D8") {
    transform_into_quads_ = true;
  } else {
    transform_into_quads_ = false;
  }

  if (ok && moving_least_squares) {
    ok = MovingLeastSquares();
    std::cout << " - MovingLeastSquares: " << ok << std::endl;
  }

  if (ok) {
    ok = Triangulate();
    std::cout << " - Triangulate: " << ok << std::endl;
  }

  if (ok && transform_into_quads_) {
    ok = TransformIntoQuads();
    std::cout << " - TransformIntoQuads: " << ok << std::endl;
    ok = InitCloud();
    std::cout << " - Cloud Re-Initialized: " << ok << std::endl;
  }

  if (ok) {
    int nodes_reversed = CheckNodeOrderConsistency();
    std::cout << " - CheckNodeOrderConsistency: " << nodes_reversed << " nodes reversed" << std::endl;
    ok = nodes_reversed != -1;
  }
  
  if (ok) {
    if (simulation) {
      SimulateFailedTriangulation();  
    }
    ok = ClearNotTriangulated();
    std::cout << " - ClearNotTriangulated: " << ok << std::endl;
  }

  return ok;
}

void FEM::ComputeExtrusion() {
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
  std::cout << " - Computed element height: mag = " << element_height_ << std::endl;
  
  // Compute second layer at a distance of 1/2 element height and with the direction of the normal vector
  points2_.clear();
  pc2_.width = pc_.width;
  pc2_.height = pc_.height;
  pc2_.points.resize(pc2_.width * pc2_.height);

  // Average normal
  Eigen::Vector3d normal;
  for (Eigen::Vector3d n : normals_) {
    normal -= n;
  }
  normal /= normals_.size();    

  // truncate normals to 3 decimals
  double scale = std::pow(10.0, 3);
  for (unsigned int i=0; i<3; i++) {
    normal(i) = std::round(normal(i) * scale) / scale;
  }
  
  //std::cout << "Element height: " << element_height_ << std::endl;
  
  for (unsigned int i=0; i<points_.size(); i++) {
    Eigen::Vector3d point2 = points_[i] - element_height_/2 * normal;
    points2_.push_back(point2);

    if (points_alive_[i]) {
      pc2_.points[i].x = point2(0);
      pc2_.points[i].y = point2(1);
      pc2_.points[i].z = point2(2);
    }

    //pc2_.points[i].x = point2(0);
    //pc2_.points[i].y = point2(1);
    //pc2_.points[i].z = point2(2);
  }

  int points_alive = std::accumulate(points_alive_.begin(), points_alive_.end(), 0);

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
  for (unsigned int i=0; i<points_.size(); i++) {
    if (points_alive_[i]) {
      int idx = indices.size() + points_alive;
      indices.push_back(idx);
    }
  }
  return indices;
}

void FEM::SetExtrusion(std::vector<Eigen::Vector3d> extrusion_delta, double element_height) {
  for (unsigned int i=0; i<points_.size(); i++) {
    Eigen::Vector3d pt2 = points_[i] + extrusion_delta[i];
    points2_.push_back(pt2);
  }
  element_height_ = element_height;
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
                   int wait) {

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ> (pc_));
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZ> (pc2_));

  pcl::visualization::PCLVisualizer viewer;
  viewer.setBackgroundColor(1, 1, 1);
  viewer.addPointCloud(cloud);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 0.0);
  
  for (unsigned int i=0; i<triangles_.size(); i++) {
    pcl::PointXYZ p0 = cloud->points[triangles_[i][0]];
    pcl::PointXYZ p1 = cloud->points[triangles_[i][1]];
    pcl::PointXYZ p2 = cloud->points[triangles_[i][2]];
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
    viewer.addSphere(cloud->points[i], 0.01, 0.0, 0.7, 0.0, "point_not_alive_" + std::to_string(i));
  }

  if (extrusion) {
    for (unsigned int i=0; i<triangles_.size(); i++) {
      pcl::PointXYZ p0 = cloud2->points[triangles_[i][0]];
      pcl::PointXYZ p1 = cloud2->points[triangles_[i][1]];
      pcl::PointXYZ p2 = cloud2->points[triangles_[i][2]];
      std::string name = "triangle_1_" + std::to_string(i) + "_extrusion";

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.9, 0.7, name + "_a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, 0.7, 0.9, 0.7, name + "_b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, 0.7, 0.9, 0.7, name + "_c");
    }

    std::cout << " Cloud points size = " << cloud->points.size() << std::endl;
    std::cout << " points alive size = " << points_alive_.size() << std::endl;
    for (unsigned int i=0; i<cloud->points.size(); i++) {
      if (!points_alive_[i]) {
        continue;
      }

      pcl::PointXYZ p0 = cloud->points[i];
      pcl::PointXYZ p1 = cloud2->points[i];
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
                   std::vector<Eigen::Vector3d> cloud2,
                   std::vector<Eigen::Vector3d> cloud2extrusion,
                   std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2,
                   int wait) {
                    
  pcl::PointCloud<pcl::PointXYZ> pc;
  pc.width = cloud2.size();
  pc.height = 1;
  pc.points.resize(pc_.width * pc_.height);
  for (unsigned int i=0; i<cloud2.size(); i++) {
    pc.points[i].x = cloud2[i](0);
    pc.points[i].y = cloud2[i](1);
    pc.points[i].z = cloud2[i](2);
  }

  ViewMesh(extrusion, pc, cloud2extrusion, pose2, wait);
}


void FEM::ViewMesh(bool extrusion, 
                   pcl::PointCloud<pcl::PointXYZ> cloud2,
                   std::vector<Eigen::Vector3d> cloud2extrusion,
                   std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2,
                   int wait) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ> (pc_));

  pcl::visualization::PCLVisualizer viewer;
  viewer.setBackgroundColor(1, 1, 1);
  viewer.addPointCloud(cloud);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 0.0);
  
  for (unsigned int i=0; i<triangles_.size(); i++) {
    pcl::PointXYZ p0 = cloud->points[triangles_[i][0]];
    pcl::PointXYZ p1 = cloud->points[triangles_[i][1]];
    pcl::PointXYZ p2 = cloud->points[triangles_[i][2]];
    std::string name = "triangle_1_" + std::to_string(i);

    // Add line
    viewer.addLine<pcl::PointXYZ>(p0, p1, 0.0, 0.7, 0.0, name + "a");
    viewer.addLine<pcl::PointXYZ>(p1, p2, 0.0, 0.7, 0.0, name + "b");
    viewer.addLine<pcl::PointXYZ>(p2, p0, 0.0, 0.7, 0.0, name + "c");
  }

  if (extrusion) {
    for (unsigned int i=0; i<triangles_.size(); i++) {
      Eigen::Vector3d p0e = points2_[triangles_[i][0]];
      Eigen::Vector3d p1e = points2_[triangles_[i][1]];
      Eigen::Vector3d p2e = points2_[triangles_[i][2]];

      pcl::PointXYZ p0(p0e(0), p0e(1), p0e(2));
      pcl::PointXYZ p1(p1e(0), p1e(1), p1e(2));
      pcl::PointXYZ p2(p2e(0), p2e(1), p2e(2));
      std::string name = "triangle_1_" + std::to_string(i) + "_extrusion";

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.9, 0.7, name + "_a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, 0.7, 0.9, 0.7, name + "_b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, 0.7, 0.9, 0.7, name + "_c");
    }

    for (unsigned int i=0; i<points2_.size(); i++) {
      pcl::PointXYZ p0 = cloud->points[i];
      pcl::PointXYZ p1(points2_[i](0), points2_[i](1), points2_[i](2));
      std::string name = "line_1_" + std::to_string(i);

      // Add line
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
      viewer.addLine<pcl::PointXYZ>(p1_a, p1_b, 0.0, 0.7, 0.0, "pose__line");
  }

  if (cloud2.points.size() > 0) {
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud2_color(cloud2.makeShared(), 0, 0, 255);
    viewer.addPointCloud(cloud2.makeShared(), cloud2_color, "cloud2");
    viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud2");

    for (unsigned int i=0; i<triangles_.size(); i++) {
      pcl::PointXYZ p0 = cloud2.points[triangles_[i][0]];
      pcl::PointXYZ p1 = cloud2.points[triangles_[i][1]];
      pcl::PointXYZ p2 = cloud2.points[triangles_[i][2]];
      std::string name = "triangle_2_" + std::to_string(i);

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.0, 0.0, name + "a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, 0.7, 0.0, 0.0, name + "b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, 0.7, 0.0, 0.0, name + "c");
    }

    if (extrusion) {
      for (unsigned int i=0; i<triangles_.size(); i++) {
        Eigen::Vector3d p0e = cloud2extrusion[triangles_[i][0]];
        Eigen::Vector3d p1e = cloud2extrusion[triangles_[i][1]];
        Eigen::Vector3d p2e = cloud2extrusion[triangles_[i][2]];

        pcl::PointXYZ p0(p0e(0), p0e(1), p0e(2));
        pcl::PointXYZ p1(p1e(0), p1e(1), p1e(2));
        pcl::PointXYZ p2(p2e(0), p2e(1), p2e(2));
        std::string name = "triangle_2_" + std::to_string(i) + "_extrusion";

        // Add line
        viewer.addLine<pcl::PointXYZ>(p0, p1, 0.9, 0.7, 0.7, name + "_a");
        viewer.addLine<pcl::PointXYZ>(p1, p2, 0.9, 0.7, 0.7, name + "_b");
        viewer.addLine<pcl::PointXYZ>(p2, p0, 0.9, 0.7, 0.7, name + "_c");
      }

      for (unsigned int i=0; i<points2_.size(); i++) {
        pcl::PointXYZ p0 = cloud2.points[i];
        pcl::PointXYZ p1(cloud2extrusion[i](0), cloud2extrusion[i](1), cloud2extrusion[i](2));
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
  }

  viewer.addCoordinateSystem(0.5);

  /*
  Eigen::Vector3d p(0.0,0.0,-1.0);
  double ang = 30*M_PI/180;
  Eigen::Vector3d axis(1,0,0);
  Eigen::Vector4d q = Eigen::Vector4d(cos(ang/2), axis(0)*sin(ang/2), axis(1)*sin(ang/2), axis(2)*sin(ang/2));
  
  double norm = q.norm();
  Eigen::Vector4d q_norm;
  if (norm == 0) {
    q_norm = Eigen::Vector4d(1.0, q(1), q(2), q(3));
  }
  else {
    q_norm = q / norm;
  }
  Eigen::Quaterniond q1(q_norm(0), q_norm(1), q_norm(2), q_norm(3));
  
  Eigen::Quaternion p_quat = Eigen::Quaternion(0.0, p(0), p(1), p(2));
  Eigen::Quaternion p1 = q1 * p_quat * q1.inverse();
  //std::cout << "Original " << p_quat << std::endl;
  //std::cout << "Rotated  " << p1 << " : " << p1.w() << " " << p1.x() << " " << p1.y() << " " << p1.z() << std::endl;
  Eigen::Vector3d p1_vec(p1.x(), p1.y(), p1.z());

  //std::cout << "Original " << p << std::endl;
  //std::cout << "Rotated  " << p1_vec << std::endl;
  pcl::PointXYZ p_o(0.0, 0.0, 0.0);
  pcl::PointXYZ p_a(p(0), p(1), p(2));
  pcl::PointXYZ p_b(p1_vec(0), p1_vec(1), p1_vec(2));
  // Add both points to viewer
  viewer.addSphere(p_a, 0.1, 0.0, 0.0, 0.0, "p_a");
  viewer.addSphere(p_b, 0.1, 0.0, 1.0, 0.0, "p_b");
  viewer.addLine<pcl::PointXYZ>(p_o, p_a, 0.0, 0.0, 0.0, "p_line_oa");
  viewer.addLine<pcl::PointXYZ>(p_o, p_b, 0.0, 0.0, 0.0, "p_line_ob");
  */

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

std::vector<Eigen::Vector3d> FEM::GetEigenNodes(bool active_only) {
  std::vector<Eigen::Vector3d> points;
  for (unsigned int i=0; i<points_.size(); i++) {
    if (!active_only || points_alive_[i]) {
      points.push_back(points_[i]);
    }
  }
  for (unsigned int i=0; i<points2_.size(); i++) {
    if (!active_only || points_alive_[i]) {
      points.push_back(points2_[i]);
    }
  }
  return points;
}

std::vector<Eigen::Vector3d> FEM::GetEigenBaseNodes() {
  return points_;
}

std::vector<Eigen::Vector3d> FEM::GetEigenExtrudedNodes() {
  return points2_;
}

std::vector<std::vector<unsigned int>> FEM::GetTriangles() {
  return triangles_;
}

void FEM::SetTriangles(std::vector<std::vector<unsigned int>> triangles) {
  triangles_ = triangles;
}

std::vector<std::vector<unsigned int>> FEM::GetElements() {
  return elements_;
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