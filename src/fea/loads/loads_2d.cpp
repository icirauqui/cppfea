#include "loads_2d.hpp"


Loads2d::Loads2d(int num_dof, std::vector<Eigen::Vector2d>* nodes) 
  : Loads(num_dof) {
  _nodes = nodes;
  _node_ids = std::vector<bool>(nodes->size(), false);
  _values = std::vector<std::vector<double>>(nodes->size(), std::vector<double>(_num_dof, 0.0));
}


void Loads2d::AddNodalByCoords(std::vector<double> coords, std::vector<bool> dof, std::vector<double> values) {
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


