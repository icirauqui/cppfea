#ifndef BOUNDARY_CONDITIONS_2D_HPP
#define BOUNDARY_CONDITIONS_2D_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

#include "boundary_conditions.hpp"

class BoundaryConditions2d : public BoundaryConditions {

public:

  BoundaryConditions2d(int num_dof, std::vector<Eigen::Vector2d>* nodes);

  // Setters
  void SetNodes(std::vector<Eigen::Vector2d>* nodes) { _nodes = nodes; }

private:

  std::vector<Eigen::Vector2d>* _nodes;

  void AddNodalByCoords(std::vector<double> coords, std::vector<bool> dof, std::vector<double> values);
  void EncastreInCoord(std::vector<double> coords, std::vector<bool> dof);


  
};



#endif