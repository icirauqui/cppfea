#ifndef LOADS_2D_HPP
#define LOADS_2D_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

#include "loads.hpp"

class Loads2d : public Loads {

public:

  Loads2d(int num_dof, std::vector<Eigen::Vector2d>* nodes);


private:

  std::vector<Eigen::Vector2d>* _nodes;

  void AddNodalByCoords(std::vector<double> coords, std::vector<bool> dof, std::vector<double> values);
  
};



#endif