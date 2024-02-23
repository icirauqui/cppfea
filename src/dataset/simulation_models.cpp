#include <iostream>
#include <vector>
#include "eigen3/Eigen/Dense"

class SimulationC3D8_1 {
public:
  SimulationC3D8_1(int z = 0): z0(z), z1(z) {
    std::cout << " - SimulationC3D8_1" << std::endl;

    for (int x=x0; x<=x1; x++) {
      for (int y=y0; y<=y1; y++) {
        _nodes.push_back(Eigen::Vector3d(x, y, z1));
      }
    }
  }

  void ReportNodes() {
    for (unsigned int i=0; i<_nodes.size(); i++) {
      std::cout << "   " << i << ": " << _nodes[i][0] << " " << _nodes[i][1] << " " << _nodes[i][2] << std::endl;
    }
  }

  std::vector<Eigen::Vector3d> GetNodes() {
    return _nodes;
  }


private:
  int x0 = -4;
  int x1 = 4;
  int y0 = -1;
  int y1 = 1;
  int z0 = 0;
  int z1 = 0;

  std::vector<Eigen::Vector3d> _nodes;
};