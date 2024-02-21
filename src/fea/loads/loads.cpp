#include "loads.hpp"


Loads::Loads(int num_dof, bool bypass) {
  _num_dof = num_dof;
  _num_nodes_constrained = 0;
  _bypass = bypass;
}

void Loads::AddNodalByNodeIds(std::vector<unsigned int> &node_ids, std::vector<double> &values) {
  for (auto node : node_ids) {
    if (_node_ids[node]) {
      std::cout << "   Warning: Node " << node << " already has a load: [";
      for (unsigned int i = 0; i < _num_dof; i++) {
        std::cout << " " << _values[node][i];
      }
      std::cout << " ]" << std::endl;
      if (!_bypass) 
        continue;
    }
    _node_ids[node] = true;
    _values[node] = values;
  }
  std::cout << "   Added loads in " << node_ids.size() << " nodes" << std::endl;
}

void Loads::AddNodalX(float x, std::vector<double> values) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[0] = x;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[0] = true;
  AddNodalByCoords(coords, dof, values);
}

void Loads::AddNodalY(float y, std::vector<double> values) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[1] = y;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[1] = true;
  AddNodalByCoords(coords, dof, values);
}

void Loads::AddNodalZ(float z, std::vector<double> values) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[2] = z;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[2] = true;
  AddNodalByCoords(coords, dof, values);
}

void Loads::AddNodalXY(float x, float y, std::vector<double> values) {
  std::vector<double> coords = {x, y};
  std::vector<bool> dof = {true, true, false};
  AddNodalByCoords(coords, dof, values);
}

void Loads::AddNodalXYZ(float x, float y, float z, std::vector<double> values) {
  std::vector<double> coords = {x, y, z};
  std::vector<bool> dof = {true, true, true};
  AddNodalByCoords(coords, dof, values);
}

void Loads::Report() {
  std::cout << std::endl;
  std::cout << "Loads3d Report" << std::endl;
  std::cout << "  _num_dof: " << _num_dof << std::endl;
  std::cout << "  _node_ids.size(): " << _node_ids.size() << std::endl;
  std::cout << "  _num_nodes_constrained: " << _num_nodes_constrained << std::endl;
  std::cout << "  _num_nodes_free: " << _node_ids.size() - _num_nodes_constrained << std::endl;

  unsigned int num_bcs = 0;
  std::cout << "  Boundary Conditions:" << std::endl;
  for (unsigned int i = 0; i < _node_ids.size(); i++) {
    if (!_node_ids[i]) 
      continue;
    num_bcs++;
      
    std::cout << "    " << i << ":";
    std::cout << " values(";
    for (unsigned int j = 0; j < _num_dof; j++) {
      std::cout << " " << _values[i][j];
    }
    std::cout << " )" << std::endl;
  }
  std::cout << "  Number of nodes with BCs: " << num_bcs << std::endl;
}


