#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

class BoundaryConditions {

public:

  BoundaryConditions(int num_dof, bool bypass = false);


  void AddNodalByNodeIds(std::vector<unsigned int> &node_ids, std::vector<double> &values);

  virtual void AddNodalByCoords(std::vector<double> coords, std::vector<bool> dof, std::vector<double> values) = 0;

  void AddNodalDisplacementByCoordX(float x, std::vector<double> values);
  void AddNodalDisplacementByCoordY(float y, std::vector<double> values);
  void AddNodalDisplacementByCoordZ(float z, std::vector<double> values);
  void AddNodalDisplacementByCoordsXY(float x, float y, std::vector<double> values);
  void AddNodalDisplacementByCoordsXYZ(float x, float y, float z, std::vector<double> values);

  void AddEncastreByNodeIds(std::vector<unsigned int> &node_ids);

  virtual void EncastreInCoord(std::vector<double> coords, std::vector<bool> dof) = 0;

  void EncastreInCoordX(float x);
  void EncastreInCoordY(float y);
  void EncastreInCoordZ(float z);
  void EncastreInCoordXY(float x, float y);
  void EncastreInCoordXYZ(float x, float y, float z);


  void Report();



  // Accessors
  bool NodeIds(unsigned int idx) { return _node_ids[idx]; }
  std::vector<bool>& NodeIds() { return _node_ids; }
  int NumDof() { return _num_dof; }
  std::vector<double>& Values(unsigned int idx) { return _values[idx]; }
  std::vector<std::vector<double>>& Values() { return _values; }


protected:

  int _num_nodes;
  int _num_dof;
  int _num_nodes_constrained;
  int _num_dof_free;
  double _tolerance = 1e-6;
  double _bypass = false;

  std::vector<bool> _node_ids;
  std::vector<std::vector<double>> _values;
};



#endif