#ifndef VIS_FEM_HPP 
#define VIS_FEM_HPP


#include <iostream>
#include <vector>

#include "fem.hpp"

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

std::pair<Eigen::Vector3d, Eigen::Vector3d> QuaternionLine(
  Eigen::Vector4d qvec, Eigen::Vector3d point, double radius = 1.0);

void ViewMesh(bool extrusion = false,
          int wait = 0,
          FEM* fems1 = nullptr,
          FEM* fems2 = nullptr);

void AddMesh(bool extrusion,
              std::vector<double> color,
              pcl::visualization::PCLVisualizer::Ptr viewer,
              FEM* fem,
              std::string fem_id);




#endif