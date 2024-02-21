#include "boundary_conditions.hpp"



BoundaryConditions::BoundaryConditions(int num_dof, bool bypass) {
  _num_dof = num_dof;
  _num_nodes_constrained = 0;
  _bypass = bypass;
}

void BoundaryConditions::AddNodalByNodeIds(std::vector<unsigned int> &node_ids, std::vector<double> &values) {
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
    _num_nodes_constrained++;
  }
  std::cout << "   Added loads in " << node_ids.size() << " nodes" << std::endl;
}

void BoundaryConditions::AddNodalDisplacementByCoordX(float x, std::vector<double> values) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[0] = x;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[0] = true;
  AddNodalByCoords(coords, dof, values);
}

void BoundaryConditions::AddNodalDisplacementByCoordY(float y, std::vector<double> values) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[1] = y;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[1] = true;
  AddNodalByCoords(coords, dof, values);
}

void BoundaryConditions::AddNodalDisplacementByCoordZ(float z, std::vector<double> values) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[2] = z;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[2] = true;
  AddNodalByCoords(coords, dof, values);
}

void BoundaryConditions::AddNodalDisplacementByCoordsXY(float x, float y, std::vector<double> values) {
  std::vector<double> coords = {x, y};
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[0] = true;
  dof[1] = true;
  AddNodalByCoords(coords, dof, values);
}

void BoundaryConditions::AddNodalDisplacementByCoordsXYZ(float x, float y, float z, std::vector<double> values) {
  std::vector<double> coords = {x, y, z};
  std::vector<bool> dof = std::vector<bool>(_num_dof, true);
  AddNodalByCoords(coords, dof, values);
}


void BoundaryConditions::EncastreInCoordX(float x) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[0] = x;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[0] = true;
  EncastreInCoord(coords, dof);
}

void BoundaryConditions::EncastreInCoordY(float y) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[1] = y;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[1] = true;
  EncastreInCoord(coords, dof);
}

void BoundaryConditions::EncastreInCoordZ(float z) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[2] = z;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[2] = true;
  EncastreInCoord(coords, dof);
}

void BoundaryConditions::EncastreInCoordXY(float x, float y) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[0] = x;
  coords[1] = y;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[0] = true;
  dof[1] = true;
  EncastreInCoord(coords, dof);
}

void BoundaryConditions::EncastreInCoordXYZ(float x, float y, float z) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[0] = true;
  dof[1] = true;
  dof[2] = true;
  EncastreInCoord(coords, dof);
}



void BoundaryConditions::Report() {
  std::cout << std::endl;
  std::cout << "BoundaryConditions3d Report" << std::endl;
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
    
    std::cout << "    " << i << ":  values(";
    for (unsigned int j = 0; j < _num_dof; j++) {
      std::cout << " " << _values[i][j];
    }
    std::cout << " )" << std::endl;
  }
  std::cout << "  Number of nodes with BCs: " << num_bcs << std::endl;
}



