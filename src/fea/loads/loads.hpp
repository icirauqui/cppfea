#ifndef LOADS_HPP
#define LOADS_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

class Loads {

public:

  Loads(int num_dof, bool bypass = false);



  void AddNodalByNodeIds(std::vector<unsigned int> &node_ids, std::vector<double> &values);

  

  virtual void AddNodalByCoords(std::vector<double> coords, std::vector<bool> dof, std::vector<double> values) = 0;

  void AddNodalX(float x, std::vector<double> values);

  void AddNodalY(float y, std::vector<double> values);

  void AddNodalZ(float z, std::vector<double> values);

  void AddNodalXY(float x, float y, std::vector<double> values);

  void AddNodalXYZ(float x, float y, float z, std::vector<double> values);



  void Report();



  // Accessors
  bool NodeIds(unsigned int idx) { return _node_ids[idx]; }
  int NumDof() { return _num_dof; }
  std::vector<bool>& NodeIds() { return _node_ids; }
  std::vector<double>& Values(unsigned int idx) { return _values[idx]; }
  std::vector<std::vector<double>>& Values() { return _values; }


protected:
  int _num_dof;
  int _num_nodes_constrained = 0;
  double _tolerance = 1e-6;
  double _bypass = false;

  std::vector<bool> _node_ids;
  std::vector<std::vector<double>> _values;
};



#endif