#ifndef FEM_HPP
#define FEM_HPP


#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <unordered_map>
#include <algorithm>

#include <eigen3/Eigen/Dense>


// PCL Libraries
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_io.h>
#include <pcl/point_types.h>

#include <pcl/common/common.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/surface/mls.h>
#include <pcl/surface/impl/mls.hpp>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/visualization/cloud_viewer.h>

#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include <pcl/surface/gp3.h>

#include <pcl/visualization/vtk.h>

#include <pcl/console/parse.h>
#include <pcl/io/vtk_lib_io.h>


std::pair<Eigen::Vector4d, Eigen::Vector3d> ApproximatePose(std::vector<Eigen::Vector3d> pts);

class FEM {

public:

  FEM(std::string element);
  FEM(FEM& fem);

  void Replicate(FEM* fem);

  void AddPoint(Eigen::Vector3d point);

  bool InitCloud();

  bool MovingLeastSquares(bool simulation = false);
  
  bool Triangulate(bool simulation = false);
  bool TriangulateByProjection(bool simulation = false);

  void SimulateFailedTriangulation();

  bool ClearNotTriangulated();

  Eigen::Vector3d calculateMidpoint(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);

  Eigen::Vector3d calculateOrthocenter(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3);

  bool TransformIntoQuads();
  
  int CheckNodeOrderConsistency();

  bool Compute(bool moving_least_squares = true, bool triangulate_planar = false, bool simulation = false);

  bool ComputeExtrusion();
  std::vector<Eigen::Vector3d> GetExtrusionDelta();
  std::vector<Eigen::Vector3d> GetExtrusion();
  std::vector<unsigned int> GetExtrusionIndices();
  void SetExtrusion(std::vector<Eigen::Vector3d> extrusion_delta, double element_height);
  void SetAliveNodes(std::vector<bool> alive_nodes);
  double GetElementHeight();
  std::vector<unsigned int> GetIndicesNotTriangulated();
  std::vector<Eigen::Vector3d> GetPoints(bool alive_only = true);




  void ViewMesh(bool extrusion = false,
                int wait = 0,
                pcl::visualization::PCLVisualizer viewer = pcl::visualization::PCLVisualizer("3D Viewer"));
    
  void ViewMesh(bool extrusion = false,
                std::vector<Eigen::Vector3d> points2 = std::vector<Eigen::Vector3d>(),
                std::vector<Eigen::Vector3d> points2extrusion = std::vector<Eigen::Vector3d>(),
                std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2 = std::make_pair(Eigen::Vector4d(0.0, 0.0, 0.0, 0.0), Eigen::Vector3d(0.0, 0.0, 0.0)),
                int wait = 0,
                pcl::visualization::PCLVisualizer viewer = pcl::visualization::PCLVisualizer("3D Viewer"));


  std::vector<std::vector<float>> GetNodes();
  std::vector<Eigen::Vector3d> GetEigenNodesFront(bool active_only = true);
  std::vector<Eigen::Vector3d> GetEigenNodes(bool active_only = true);
  //std::vector<Eigen::Vector3d> GetEigenNodes(bool active_only = true, bool is_target = false);
  std::vector<Eigen::Vector3d> GetEigenBaseNodes();
  std::vector<Eigen::Vector3d> GetEigenExtrudedNodes();
  std::vector<bool> GetPointsAlive();

  std::vector<std::vector<unsigned int>> GetTriangles();
  void SetTriangles(std::vector<std::vector<unsigned int>> triangles);

  std::vector<std::vector<unsigned int>> GetElements(bool alive_only = true);
  void SetElements(std::vector<std::vector<unsigned int>> elements);

  pcl::PointCloud<pcl::PointXYZ> GetCloud();
  pcl::PointCloud<pcl::PointXYZ> GetCloud2();
  std::pair<Eigen::Vector4d, Eigen::Vector3d> GetPose();
  void SetPose(std::pair<Eigen::Vector4d, Eigen::Vector3d> pose);

  void Transform(std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes,
                 std::pair<Eigen::Vector4d, Eigen::Vector3d> pose);




  std::pair<Eigen::Vector3d, Eigen::Vector3d> QuaternionLine(
    Eigen::Vector4d qvec, Eigen::Vector3d point, double radius = 1.0);

private:

  std::pair<Eigen::Vector3d, Eigen::Vector3d> QuaternionLine2(
    Eigen::Vector4d qvec, Eigen::Vector3d point, double radius = 1.0);

  std::string element_;

  std::vector<Eigen::Vector3d> points_, points2_;
  std::vector<bool> points_alive_;
  std::vector<unsigned int> indices_not_triangulated_;

  pcl::PointCloud<pcl::PointXYZ> pc_, pc2_;

  std::vector<std::vector<unsigned int>> triangles_, elements_, quadrilaterals_;

  bool transform_into_quads_ = false;

  std::vector<int> mls_indices_, points_indices_;

  std::vector<Eigen::Vector3d> normals_;

  pcl::PolygonMesh mesh_;

  double element_height_ = 0.0;

  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_;

  std::vector<std::vector<double>> colorsf_, colorsb_;

  // Interface 
  float mls_search_radius_ = 1.0;
  int mls_polynomial_order_ = 3;
  float mesh_mu_ = 2.5;
  float mesh_search_radius_ = 10.0;
  int mesh_max_neighbours_ = 25;
  int mesh_surf_angle_ = 150;
  int mesh_min_angle_ = 5;
  int mesh_max_angle_ = 85;


};



void ViewMesh(bool extrusion,
              int wait,
              std::vector<FEM*> fems);












#endif