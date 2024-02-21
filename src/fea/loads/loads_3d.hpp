#ifndef LOADS_3D_HPP
#define LOADS_3D_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

#include "loads.hpp"

class Loads3d : public Loads {

public:

  Loads3d(int num_dof, std::vector<Eigen::Vector3d>* nodes);


private:

  std::vector<Eigen::Vector3d>* _nodes;
  
  void AddNodalByCoords(std::vector<double> coords, std::vector<bool> dof, std::vector<double> values);

};



#endif