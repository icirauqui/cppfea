#include "boundary_conditions_2d.hpp"


BoundaryConditions2d::BoundaryConditions2d(int num_dof, std::vector<Eigen::Vector2d>* nodes) 
  : BoundaryConditions(num_dof) {

  _num_dof_free = nodes->size();
  _nodes = nodes;
  _num_nodes = nodes->size();

  _node_ids = std::vector<bool>(_num_dof_free, false);
  _values = std::vector<std::vector<double>>(_num_dof_free, std::vector<double>(_num_dof, 0.0));
}


void BoundaryConditions2d::AddNodalByCoords(std::vector<double> coords, std::vector<bool> dof, std::vector<double> values) {
  std::vector<unsigned int> nodes;

  for (unsigned int node = 0; node < _nodes->size(); node++) {
    for (unsigned int d = 0; d < _num_dof; d++) {
      if (!dof[d])
        continue;

      if (abs((*_nodes)[node](d) - coords[d]) < _tolerance) {
        nodes.push_back(node);
        break;
      }
    }
  }

  AddNodalByNodeIds(nodes, values);
}


void BoundaryConditions2d::EncastreInCoord(std::vector<double> coords, std::vector<bool> dof) {
  std::vector<double> values = std::vector<double>(_num_dof, 0.0);
  std::vector<unsigned int> nodes;

  for (unsigned int node = 0; node < _num_nodes; node++) {
    for (unsigned int d = 0; d < _num_dof; d++) {
      if (!dof[d])
        continue;

      if (abs((*_nodes)[node](d) - coords[d]) < _tolerance) {
        nodes.push_back(node);
        break;
      }
    }
  }

  AddNodalByNodeIds(nodes, values);
}