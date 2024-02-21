#ifndef VIS_HPP
#define VIS_HPP

#include <pcl/common/common_headers.h>
#include <pcl/console/parse.h>
#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/pcl_visualizer.h>

#if VTK_MAJOR_VERSION > 8
#include <vtk-9.2/QVTKOpenGLNativeWidget.h>
#include <vtk-9.2/vtkRenderWindow.h>
#else
#include <vtk-7.1/QVTKWidget.h>
#include <vtk-7.1/vtkRenderWindow.h>
#endif




class PCLViewer {

public:

  PCLViewer();

  void initializeViewer();

  void AddNodes(std::vector<Eigen::Vector3d> &pts, std::string name, Eigen::Vector3d color);
  void AddEdges(std::vector<std::vector<unsigned int>> &elts, std::string elt);
  void AddBCs(std::vector<bool> &nodes, std::vector<std::vector<double>> &mag, bool isLoad = false);

  void Render(double scale = 1.0, bool contours_only = false);

private:

  pcl::visualization::PCLVisualizer::Ptr viewer_;
  std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> clouds_;
  std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> clouds_countour_;
  std::vector<std::string> cloud_names_;
  std::vector<Eigen::Vector3d> cloud_colors_;

  std::vector<std::tuple<unsigned int, unsigned int>> edges_;
  std::vector<std::tuple<unsigned int, unsigned int>> edges_contour_;

  std::vector<int> load_nodes_;
  std::vector<std::vector<double>> load_mags_;
  std::vector<bool> is_load_;

  Eigen::Vector3d contour_bound_min_, contour_bound_max_;
  std::vector<bool> contour_nodes_;
};





#endif