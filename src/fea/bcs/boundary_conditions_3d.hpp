#ifndef BOUNDARY_CONDITIONS_3D_HPP
#define BOUNDARY_CONDITIONS_3D_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

#include "boundary_conditions.hpp"

class BoundaryConditions3d : public BoundaryConditions {

public:

  BoundaryConditions3d(int num_dof, std::vector<Eigen::Vector3d>* nodes);

  // Setters
  void SetNodes(std::vector<Eigen::Vector3d>* nodes) { _nodes = nodes; }

private:

  std::vector<Eigen::Vector3d>* _nodes;

  void AddNodalByCoords(std::vector<double> coords, std::vector<bool> dof, std::vector<double> values);
  void EncastreInCoord(std::vector<double> coords, std::vector<bool> dof);

};



#endif