#ifndef FEA_HPP
#define FEA_HPP

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <math.h>
#include <chrono>
#include <bits/stdc++.h>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/LU>
#include <limits>
#include <cmath>

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

#include <elts/element.hpp>
#include <elts/element2d.hpp>
#include <elts/c2d4.hpp>
#include <elts/element3d.hpp>
#include <elts/c3d6.hpp>
#include <elts/c3d8.hpp>

#include <bcs/boundary_conditions.hpp>
#include <bcs/boundary_conditions_2d.hpp>
#include <bcs/boundary_conditions_3d.hpp>

#include <loads/loads.hpp>
#include <loads/loads_2d.hpp>
#include <loads/loads_3d.hpp>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/LU>
#include <limits>
#include <cmath>

#include "solver.hpp"



class FEA {

public:


  FEA(std::string element_type,
      float young_modulus, float poisson_coefficient,
      bool debug_mode);
    
  void MatAssembly(std::vector<Eigen::Vector2d> &vpts, 
                   std::vector<std::vector<unsigned int>> &velts);
    
  void MatAssembly(std::vector<Eigen::Vector3d> &vpts, 
                   std::vector<std::vector<unsigned int>> &velts);

  void ApplyBoundaryConditions(BoundaryConditions &bc);
  void ApplyLoads(Loads &loads);

  /**
   * Solves the finite element system, using K and F to get U.
   *
   * @param method Solver to use: CG, BiCGSTAB, LU, LUFull.
   * @return <void>, U is stored in the object.
   */
  void Solve(std::string method = "LU");

  void PostProcess(std::vector<Eigen::Vector2d> &vpts, 
                   std::vector<std::vector<unsigned int>> &velts);

  void PostProcess(std::vector<Eigen::Vector3d> &vpts, 
                   std::vector<std::vector<unsigned int>> &velts);

  double ComputeStrainEnergy();
  double ComputeStrainEnergy(std::vector<Eigen::Vector3d> &u0,
                             std::vector<Eigen::Vector3d> &u1);

  
  // Legacy fea

  void EncastreBackLayer();
  void ImposeDirichletEncastre(std::vector<int> &dir);
  void ImposeDirichletEncastre(std::vector<std::vector<int>> &dir);

  void ComputeForces();
  void SetForces(std::vector<std::vector<float>> &vF);

  void ComputeDisplacements();



  // Reports

  void ReportNodes(std::string filename);
  void ExportAll(std::string filename);
  void ExportSystem(std::string filename);
  void ExportK(std::string filename);
  void ExportF(std::string filename);
  void ExportU(std::string filename);
  void PrintK();
  void PrintEigenvaluesK();
  void ReportFEAData(std::string filename);

  // Accessors

  Eigen::MatrixXd K() { return K_; }
  Eigen::MatrixXd F() { return F_; }
  Eigen::MatrixXd U() { return U_; }
  float StrainEnergy() { return sE_; }
  int NumDof() { return element_->getDofPerNode(); }





private:

  
  Element* element_;
  //BoundaryConditions* bcs_;
  //Loads* loads_;
  FEAData* fea_data_;

  Eigen::MatrixXd K_;
  Eigen::MatrixXd K1_;
  Eigen::MatrixXd Kei_;

  Eigen::MatrixXd F_, Fi_;
  Eigen::MatrixXd U_;

  std::vector<bool> bc_;
  std::vector<bool> loads_;

  float k_large_ = 1e8;

  float sE_ = 0.0;

  float E_ = 1.0;
  float nu_ = 0.499;

  bool debug_mode_ = false;
};

#endif