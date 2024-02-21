#ifndef PCLVIEWER_H
#define PCLVIEWER_H

#include <pcl/common/common_headers.h>
#include <pcl/console/parse.h>
#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <eigen3/Eigen/Core>
#include <iostream>
#include <thread>
#include <vector>

#if VTK_MAJOR_VERSION > 8
#include <vtk-9.2/QVTKOpenGLNativeWidget.h>
#include <vtk-9.2/vtkRenderWindow.h>
#else
#include <vtk-7.1/QVTKWidget.h>
#include <vtk-7.1/vtkRenderWindow.h>
#endif

#include "../aux/pose.cpp"
#include "../core/image.h"

// #include <chrono>
using namespace std::chrono_literals;

class PCLViewer {
 public:
  struct PCLImage {
   public:
    Eigen::Vector4d q_;
    Eigen::Vector3d t_;
    int modelId_;

    pcl::PointCloud<pcl::PointXYZ>::Ptr poly_;

    PCLImage(Eigen::Vector4d q, Eigen::Vector3d t, int modelId)
        : q_(q),
          t_(t),
          modelId_(modelId),
          poly_(new pcl::PointCloud<pcl::PointXYZ>) {}

    Eigen::Vector4d get_q() {
      Eigen::Vector4d q = q_;
      return q;
    }

    Eigen::Vector3d get_t() {
      Eigen::Vector3d t = t_;
      return t;
    }

    double get_q(int idx) { return q_(idx); }

    double get_t(int idx) { return t_(idx); }
  };

  struct PCLTriangle {
   public:
    int modelId_;

    pcl::PointCloud<pcl::PointXYZ>::Ptr poly_;

    PCLTriangle(int modelId)
        : modelId_(modelId), poly_(new pcl::PointCloud<pcl::PointXYZ>) {}
  };

  inline PCLViewer(unsigned int nPointClouds = 3);

  inline ~PCLViewer();

  inline void initializeViewer();

  inline void registerWindowSize(int width, int height);

  // IMAGEs
  inline void addImages(std::vector<Image *> images, int modelId);
  inline void addImages2(std::vector<Image *> images, int modelId);
  inline void clearImages(int modelId);

  // POINTs3D
  inline void addPointCloud(std::vector<Eigen::Vector3d> points,
                            int cloudId = 0);
  inline void addFusedPointCloud(std::vector<Eigen::Vector3d> points,
                                 std::vector<Eigen::Vector3d> colors,
                                 std::vector<std::pair<double, double>> points_range,
                                 std::vector<std::vector<Eigen::Vector3d>> xyzs,
                                 int cloudId = 0);
  inline void delPointCloud(int cloudId);

  // POINTs3Dd
  inline void addPointCloudDef(std::vector<Eigen::Vector3d> points,
                               int cloudId = 0);
  inline void addFusedPointCloudDef(std::vector<Eigen::Vector3d> points,
                                    std::vector<std::pair<double, double>> points_range,
                                    std::vector<Eigen::Vector3d> xyzs,
                                    int cloudId = 0);
  inline void delPointCloudDef(int cloudId);

  // Mesh
  inline void AddMesh(std::vector<std::vector<int>> mesh, int cloudId = 1);
  inline void DelMesh(int cloudId);

  // Rendering
  inline void render();
  inline void render_images(pcl::visualization::PCLVisualizer::Ptr viewer);
  inline void render_clouds_p3D(pcl::visualization::PCLVisualizer::Ptr viewer);
  inline void render_clouds_p3Dd(pcl::visualization::PCLVisualizer::Ptr viewer);
  inline void render_meshes(pcl::visualization::PCLVisualizer::Ptr viewer);

  pcl::visualization::PCLVisualizer::Ptr viewer_;

  // Display mode
  inline void DisplayModeMesh();
  inline void DisplayModeAnalyze();
  inline void DisplayModeCompare();
  inline void DisplayModeAgg();

  bool display_images_1_ = true;
  bool display_images_2_ = true;
  bool display_images_agg_ = true;

  bool display_points3D_1_ = true;
  bool display_points3D_2_ = true;
  bool display_points3D_agg_ = true;

  bool display_points3Dd_1_ = false;
  bool display_points3Dd_2_ = false;
  bool display_points3Dd_agg_ = false;

  bool display_agg_spheres_ = false;
  bool display_agg_trajectories_ = false;

  bool display_meshes_ = false;
  bool display_analyze_ = true;
  bool display_compare_ = false;
  bool display_agg_ = false;

  bool display_mesh_tri_ = true;
  bool display_mesh_pts_ = true;

 private:
  pcl::visualization::PCLVisualizerInteractorStyle *style_;

  int win_width_1_ = 144;
  int win_height_1_ = 90;

  double sphere_transparency_avg_ = 0.10;
  double sphere_transparency_max_ = 0.05;

  std::vector<std::vector<PCLImage *>> images_;
  std::vector<std::vector<PCLTriangle *>> meshes_;

  std::vector<pcl::PointCloud<pcl::PointXYZ>::ConstPtr> clouds_images_;
  std::vector<pcl::PointCloud<pcl::PointXYZ>::ConstPtr> clouds_p3D_;
  std::vector<pcl::PointCloud<pcl::PointXYZ>::ConstPtr> clouds_p3Dd_;
  std::vector<pcl::PointCloud<pcl::PointXYZ>::ConstPtr> clouds_mesh_;

  std::vector<std::vector<Eigen::Vector3d>> clouds_p3D_colors_;
  std::vector<std::vector<Eigen::Vector3d>> clouds_p3Dd_colors_;

  std::vector<std::pair<double, double>> clouds_p3D_range_;
  std::vector<std::pair<double, double>> clouds_p3Dd_range_;

  std::vector<std::vector<Eigen::Vector3d>> clouds_p3D_xyzs_;
  std::vector<std::vector<Eigen::Vector3d>> clouds_p3Dd_xyzs_;
};

#endif
