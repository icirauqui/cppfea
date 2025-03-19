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





class SimulationMeshmeld01 {
public:
  SimulationMeshmeld01(std::string path0, std::string path1) {
    std::cout << " - SimulationMeshmeld01" << std::endl;

    // Load nodes from path0
    std::ifstream file(path0);
    std::string line;
    double x, y, z;

    while (std::getline(file, line)) {
      std::istringstream iss(line);
      if (iss >> x >> y >> z) {
        _nodes0.push_back(Eigen::Vector3d(x, y, z));
      }
    }
    file.close();

    // Load nodes from path1
    file.open(path1);
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      if (iss >> x >> y >> z) {
        _nodes1.push_back(Eigen::Vector3d(x, y, z));
      }
    }
    file.close();
  }

  void ReportNodes0() {
    std::cout << std::endl << "Nodes0" << std::endl;
    for (unsigned int i=0; i<_nodes0.size(); i++) {
      std::cout << "   " << i << ": " << _nodes0[i][0] << " " << _nodes0[i][1] << " " << _nodes0[i][2] << std::endl;
    }
  }

  void ReportNodes1() {
    std::cout << std::endl << "Nodes1" << std::endl;
    for (unsigned int i=0; i<_nodes1.size(); i++) {
      std::cout << "   " << i << ": " << _nodes1[i][0] << " " << _nodes1[i][1] << " " << _nodes1[i][2] << std::endl;
    }
  }

  std::vector<Eigen::Vector3d> GetNodes0() { return _nodes0; }
  std::vector<Eigen::Vector3d> GetNodes1() { return _nodes1; }


private:

  std::vector<Eigen::Vector3d> _nodes0, _nodes1;
};