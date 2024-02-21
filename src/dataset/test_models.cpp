#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "eigen3/Eigen/Dense"


template <typename T>
void printVector(std::string title, std::vector<T> &v) {
  std::cout << title << ":";
  for (auto i : v) {
    std::cout << " " << i;
  }
  std::cout << std::endl;
}


class AbaqusC3D8_1 {
public:
  AbaqusC3D8_1() {
    unsigned int n = 0;
    
    for (int x=int(_x0); x>=0; x-=int(_w)) {
      for (int z=int(_z0); z>=0; z-=int(_w)) {
        for (int y=int(_y0); y>=0; y-=int(_w), n++) {
          //std::cout << "n = " << n << " : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
          _nodes.push_back(Eigen::Vector3d(double(x),double(y),double(z)));

          if (x==0 || y==0 || z==0) continue;

          unsigned int offset = (_y0+1)*(_z0+1);//*(_x0+1-x);
          //std::cout << "n (offset) = " << n << " (" << offset << ") : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
          
          std::vector<unsigned int> element = {
            n+1+(int(_y0)+1)+offset,  //11
            n+1+(int(_y0)+1),         //3
            n+(int(_y0)+1),           //2
            n+(int(_y0)+1)+offset,    //10

            n+1+offset,               //9
            n+1,                      //1
            n,                        //0
            n+offset,                 //8
          };

          //std::vector<unsigned int> element = {
          //  n,
          //  n+1,
          //  n+1+(int(_y0)+1),
          //  n+(int(_y0)+1),
          //  n+offset,
          //  n+1+offset,
          //  n+1+(int(_y0)+1)+offset,
          //  n+(int(_y0)+1)+offset
          //};

          _elements.push_back(element);
        }
      }
    } 
  }

  void ApplyDisplacements(std::vector<Eigen::Vector3d> &displacements, double scale = 1.0) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      _nodes_deformed.push_back(_nodes[n] + scale * displacements[n]);
    }
  }
    
  std::vector<Eigen::Vector3d> _nodes, _nodes_deformed;
  std::vector<std::vector<unsigned int>> _elements;

  double _x0 = 5.0;
  double _y0 = 5.0;
  double _z0 = 25.0;
  double _w = 1.0;

  std::string Name() { return "abaqus_c3d8_1"; }
  double LoadLocation() { return _z0; }
  std::string ElementType() { return "C3D8"; }
};



class AbaqusC3D8_2 {
public:
  AbaqusC3D8_2() {
    unsigned int n = 0;
    
    for (int x=int(_x0); x>=0; x-=int(_w)) {
      for (int z=int(_z0); z>=0; z-=int(_w)) {
        for (int y=int(_y0); y>=0; y-=int(_w), n++) {
          //std::cout << "n = " << n << " : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
          _nodes.push_back(Eigen::Vector3d(double(x),double(y),double(z)));

          if (x==0 || y==0 || z==0) continue;

          unsigned int offset = (_y0+1)*(_z0+1);//*(_x0+1-x);
          //std::cout << "n (offset) = " << n << " (" << offset << ") : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;

          std::vector<unsigned int> element = {
            n+1+(int(_y0)+1)+offset,  //11
            n+1+(int(_y0)+1),         //3
            n+(int(_y0)+1),           //2
            n+(int(_y0)+1)+offset,    //10

            n+1+offset,               //9
            n+1,                      //1
            n,                        //0
            n+offset,                 //8
          };

          //std::vector<unsigned int> element = {
          //  n,
          //  n+1,
          //  n+1+(int(_y0)+1),
          //  n+(int(_y0)+1),
          //  n+offset,
          //  n+1+offset,
          //  n+1+(int(_y0)+1)+offset,
          //  n+(int(_y0)+1)+offset
          //};

          _elements.push_back(element);
        }
      }
    } 
  }

  void ApplyDisplacements(std::vector<Eigen::Vector3d> &displacements, double scale = 1.0) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      _nodes_deformed.push_back(_nodes[n] + scale * displacements[n]);
    }
  }
    
  std::vector<Eigen::Vector3d> _nodes, _nodes_deformed;
  std::vector<std::vector<unsigned int>> _elements;

  double _x0 = 1.0;
  double _y0 = 1.0;
  double _z0 = 5.0;
  double _w = 1.0;

  std::string Name() { return "abaqus_c3d8_2"; }
  double LoadLocation() { return _z0; }
  std::string ElementType() { return "C3D8"; }
};



class AbaqusC3D8_3 {
public:
  AbaqusC3D8_3() {
    unsigned int n = 0;
    
    for (int x=int(_x0); x>=0; x-=int(_w)) {
      for (int z=int(_z0); z>=0; z-=int(_w)) {
        for (int y=int(_y0); y>=0; y-=int(_w), n++) {
          //std::cout << "n = " << n << " : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
          _nodes.push_back(Eigen::Vector3d(double(x),double(y),double(z)));

          if (x==0 || y==0 || z==0) continue;

          unsigned int offset = (_y0+1)*(_z0+1);//*(_x0+1-x);
          //std::cout << "n (offset) = " << n << " (" << offset << ") : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;

          std::vector<unsigned int> element = {
            n+1+(int(_y0)+1)+offset,  //11
            n+1+(int(_y0)+1),         //3
            n+(int(_y0)+1),           //2
            n+(int(_y0)+1)+offset,    //10

            n+1+offset,               //9
            n+1,                      //1
            n,                        //0
            n+offset,                 //8
          };

          /*
            n+1+offset,               //9
            n+1,                      //1
            n+1+(int(_y0)+1),         //3
            n+1+(int(_y0)+1)+offset,  //11
            
            n+offset,                 //8
            n,                        //0
            n+(int(_y0)+1),           //2
            n+(int(_y0)+1)+offset     //10
            */

          printVector("Element " + std::to_string(n), element);

          _elements.push_back(element);
        }
      }
    } 
  }

  void ApplyDisplacements(std::vector<Eigen::Vector3d> &displacements, double scale = 1.0) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      _nodes_deformed.push_back(_nodes[n] + scale * displacements[n]);
    }
  }
    
  std::vector<Eigen::Vector3d> _nodes, _nodes_deformed;
  std::vector<std::vector<unsigned int>> _elements;

  double _x0 = 1.0;
  double _y0 = 1.0;
  double _z0 = 3.0;
  double _w = 1.0;

  std::string Name() { return "abaqus_c3d8_3"; }
  double LoadLocation() { return _z0; }
  std::string ElementType() { return "C3D8"; }
};



class AbaqusC2D4_1 {
public:
  AbaqusC2D4_1(std::string filename) {
    read_inp_file(filename);
  }

  // Assuming `local_error` is a function that handles errors. You'll need to define it.
  void local_error(const std::string &message) {
      std::cerr << message << std::endl;
      // Possibly throw an exception or handle the error appropriately
  }

  void read_inp_file(const std::string& inpFileName) {
    std::cout << "\n** Read input file" << std::endl;

    std::ifstream inpFile(inpFileName);
    if (!inpFile.is_open()) {
        local_error("Unable to open file");
        return;
    }

    std::string line;
    int state = 0;
    while (getline(inpFile, line)) {
      // Trim leading and trailing whitespace
      line.erase(0, line.find_first_not_of(" \t"));
      line.erase(line.find_last_not_of(" \t") + 1);

      if (line.empty()) continue;
      if (line[0] == '*') {
        state = 0;
      }
      if (line.length() >= 5 && line.substr(0, 5).compare("*Node") == 0) {
          state = 1;
          continue;
      }
      if (line.length() >= 8 && line.substr(0, 8).compare("*Element") == 0) {
          state = 2;
          continue;
      }
      if (line.length() >= 9 && line.substr(0, 9).compare("*Boundary") == 0) {
          state = 3;
          continue;
      }

      //std::cout << "state: " << state << " : " << line << std::endl;

      if (state == 0) {
          continue;
      }

      std::istringstream iss(line);
      std::string token;
      std::vector<std::string> tokens;
      while (getline(iss, token, ',')) {
          tokens.push_back(token);
      }

      if (state == 1) {
        // Read nodes
        if (tokens.size() != 3) {
          local_error("A node definition needs 3 values");
          continue;
        }
        int nodeNr = std::stoi(tokens[0]) - 1; // zero indexed
        double xx = std::stod(tokens[1]);
        double yy = std::stod(tokens[2]);
        _nodes.push_back(Eigen::Vector2d(xx, yy)); // assume the nodes are ordered 1, 2, 3...
      } else if (state == 2) {
        // Read elements
        if (tokens.size() != 5) {
          local_error("An element definition needs 5 values");
          continue;
        }
        int elemNr = std::stoi(tokens[0]); // not used in this context
        std::vector<unsigned int> element;
        for (int i = 1; i <= 4; ++i) {
          element.push_back(std::stoi(tokens[i]) - 1); // zero indexed
        }
        // Assuming the order needs to be adjusted from the Python code
        _elements.push_back({element[0], element[1], element[2], element[3]});
        //_elements.push_back({element[0], element[3], element[2], element[1]});
      } else if (state == 3) {
        // Read displacement boundary conditions
        if (tokens.size() != 4) {
          local_error("A displacement boundary condition needs 4 values");
          continue;
        }
        int nodeNr = std::stoi(tokens[0]) - 1; // zero indexed
        int dof1 = std::stoi(tokens[1]);
        int dof2 = std::stoi(tokens[2]);
        double val = std::stod(tokens[3]);
        if (dof1 == 1) {
          boundary.push_back({static_cast<double>(nodeNr), 1.0, val});
        }
        if (dof2 == 2) {
          boundary.push_back({static_cast<double>(nodeNr), 2.0, val});
        }
      }
    }

    inpFile.close();

    std::cout << "Nodes: " << _nodes.size() << std::endl;
    std::cout << "Elements: " << _elements.size() << std::endl;
    std::cout << "Boundary: " << boundary.size() << std::endl;

  }


  void ApplyDisplacements(std::vector<Eigen::Vector2d> &displacements, double scale = 1.0) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      _nodes_deformed.push_back(_nodes[n] + scale * displacements[n]);
    }
  }
    
  std::vector<Eigen::Vector2d> _nodes, _nodes_deformed;
  std::vector<std::vector<unsigned int>> _elements;
  std::vector<std::vector<double>> boundary;

  double _x0 = 1.0;
  double _y0 = 1.0;
  double _z0 = 3.0;
  double _w = 1.0;

  std::string Name() { return "abaqus_c2d4_1"; }
  double LoadLocation() { return _z0; }
  std::string ElementType() { return "C2D4"; }
};




class AbaqusC3D6_1 {
public:
  AbaqusC3D6_1(std::string filename) {
    read_inp_file(filename);
  }

  // Assuming `local_error` is a function that handles errors. You'll need to define it.
  void local_error(const std::string &message) {
      std::cerr << message << std::endl;
      // Possibly throw an exception or handle the error appropriately
  }

  void read_inp_file(const std::string& inpFileName) {
    std::cout << "\n** Read input file" << std::endl;

    std::ifstream inpFile(inpFileName);
    if (!inpFile.is_open()) {
        local_error("Unable to open file");
        return;
    }

    std::string line;

    while (getline(inpFile, line)) {
      // Trim leading and trailing whitespace
      line.erase(0, line.find_first_not_of(" \t"));
      line.erase(line.find_last_not_of(" \t") + 1);

      if (line.empty()) continue;

      std::istringstream iss(line);
      std::string token;
      std::vector<std::string> tokens;
      while (getline(iss, token, ',')) {
          tokens.push_back(token);
      }

      // Read nodes
      if (tokens.size() != 4) {
        local_error("A node definition needs 4 values");
        continue;
      }
      int nodeNr = std::stoi(tokens[0]); // zero indexed
      double xx = std::stod(tokens[1]);
      double yy = std::stod(tokens[2]);
      double zz = std::stod(tokens[3]);
      _nodes.push_back(Eigen::Vector3d(xx, yy, zz)); // assume the nodes are ordered 1, 2, 3...
    }

    inpFile.close();

    std::cout << "Nodes: " << _nodes.size() << std::endl;

    std::vector<unsigned int> element1 = {1,3,13,0,2,12};
    std::vector<unsigned int> element2 = {15,13,3,14,12,2};
    _elements.push_back(element1);
    _elements.push_back(element2);


    for (unsigned int i=2; i<_num_elements; i++) {
      std::vector<unsigned int> element;
      for (unsigned int j=0; j<_elements[i-2].size(); j++) {
        element.push_back(_elements[i-2][j] + 2);
      }
      _elements.push_back(element);
    }



    std::cout << "Elements: " << _elements.size() << std::endl;
    for (auto elt: _elements) {
      std::cout << "Element: ";
      for (auto n: elt) {
        std::cout << n << " ";
      }
      std::cout << std::endl;
    }
  }


  void ApplyDisplacements(std::vector<Eigen::Vector3d> &displacements, double scale = 1.0) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      _nodes_deformed.push_back(_nodes[n] + scale * displacements[n]);
    }
  }

  void ApplyDisplacements(std::vector<Eigen::Vector3d> &displacements, std::vector<double> scale) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      Eigen::Vector3d node = _nodes[n];
      for (unsigned int i=0; i<3; i++) {
        node[i] += scale[i] * displacements[n][i];
      }
      _nodes_deformed.push_back(node);
    }
  }
    
  std::vector<Eigen::Vector3d> _nodes, _nodes_deformed;
  std::vector<std::vector<unsigned int>> _elements;

  int _num_elements = 10;
  double _x0 = 1.0;
  double _y0 = 1.0;
  double _z0 = 5.0;
  double _w = 1.0;

  std::string Name() { return "c3d6_1"; }
  double LoadLocation() { return _z0; }
  std::string LoadDirection() { return "z"; }
  int LoadDirectionKey() { return 3; }
  double LoadMagnitude() { return 1.0; }
  std::string ElementType() { return "C3D6"; }
};








class AbaqusC3D6_2 {
public:
  AbaqusC3D6_2(std::string filename) {
    read_inp_file(filename);
  }

  // Assuming `local_error` is a function that handles errors. You'll need to define it.
  void local_error(const std::string &message) {
      std::cerr << message << std::endl;
      // Possibly throw an exception or handle the error appropriately
  }

  void read_inp_file(const std::string& inpFileName) {
    std::cout << "\n** Read input file" << std::endl;

    std::ifstream inpFile(inpFileName);
    if (!inpFile.is_open()) {
        local_error("Unable to open file");
        return;
    }

    std::string line;
    int state = 0;
    while (getline(inpFile, line)) {
      // Trim leading and trailing whitespace
      line.erase(0, line.find_first_not_of(" \t"));
      line.erase(line.find_last_not_of(" \t") + 1);

      if (line.empty()) continue;
      if (line[0] == '*') {
        state = 0;
      }
      if (line.length() >= 5 && line.substr(0, 5).compare("*Node") == 0) {
          state = 1;
          continue;
      }
      if (line.length() >= 8 && line.substr(0, 8).compare("*Element") == 0) {
          state = 2;
          continue;
      }
      if (line.length() >= 9 && line.substr(0, 9).compare("*Boundary") == 0) {
          state = 3;
          continue;
      }

      //std::cout << "state: " << state << " : " << line << std::endl;

      if (state == 0) {
          continue;
      }

      std::istringstream iss(line);
      std::string token;
      std::vector<std::string> tokens;
      while (getline(iss, token, ',')) {
          tokens.push_back(token);
      }

      if (state == 1) {
        // Read nodes
        if (tokens.size() != 4) {
          local_error("A node definition needs 4 values");
          continue;
        }
        int nodeNr = std::stoi(tokens[0]) - 1; // zero indexed
        double xx = std::stod(tokens[1]) - 4.0;
        double yy = std::stod(tokens[2]) - 4.0;
        double zz = std::stod(tokens[3]);
        _nodes.push_back(Eigen::Vector3d(xx, yy, zz)); // assume the nodes are ordered 1, 2, 3...
      } else if (state == 2) {
        // Read elements
        if (tokens.size() != 7) {
          local_error("An element definition needs 7 values");
          continue;
        }
        int elemNr = std::stoi(tokens[0]); // not used in this context
        std::vector<unsigned int> element;
        for (int i = 1; i <= 6; ++i) {
          element.push_back(std::stoi(tokens[i]) - 1); // zero indexed
        }
        // Assuming the order needs to be adjusted from the Python code
        _elements.push_back({element[0], element[1], element[2], element[3], element[4], element[5]});
        //_elements.push_back({element[0], element[3], element[2], element[1]});
      } else if (state == 3) {
        std::cout << "BCs not implemented in file read" << std::endl;
        // Read displacement boundary conditions
        //if (tokens.size() != 4) {
        //  local_error("A displacement boundary condition needs 4 values");
        //  continue;
        //}
        //int nodeNr = std::stoi(tokens[0]) - 1; // zero indexed
        //int dof1 = std::stoi(tokens[1]);
        //int dof2 = std::stoi(tokens[2]);
        //double val = std::stod(tokens[3]);
        //if (dof1 == 1) {
        //  boundary.push_back({static_cast<double>(nodeNr), 1.0, val});
        //}
        //if (dof2 == 2) {
        //  boundary.push_back({static_cast<double>(nodeNr), 2.0, val});
        //}
      }
    }

    inpFile.close();

    std::cout << "Nodes: " << _nodes.size() << std::endl;
    std::cout << "Elements: " << _elements.size() << std::endl;
    //std::cout << "Boundary: " << boundary.size() << std::endl;

  }


  void ApplyDisplacements(std::vector<Eigen::Vector3d> &displacements, double scale = 1.0) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      _nodes_deformed.push_back(_nodes[n] + scale * displacements[n]);
      std::cout << "Node " << n << ": " << _nodes[n].transpose() << " -> " << _nodes_deformed[n].transpose() << std::endl;
    }
  }
    
  std::vector<Eigen::Vector3d> _nodes, _nodes_deformed;
  std::vector<std::vector<unsigned int>> _elements;

  int _num_elements = 10;
  double _x0 = 1.0;
  double _y0 = 1.0;
  double _z0 = 5.0;
  double _w = 1.0;

  std::string Name() { return "c3d6_2"; }
  double LoadLocation() { return _z0; }
  std::string LoadDirection() { return "z"; }
  int LoadDirectionKey() { return 3; }
  double LoadMagnitude() { return 1.0; }
  std::string ElementType() { return "C3D6"; }
};





class AbaqusC3D6_3 {
public:
  AbaqusC3D6_3(std::string filename) {
    read_inp_file(filename);
  }

  // Assuming `local_error` is a function that handles errors. You'll need to define it.
  void local_error(const std::string &message) {
      std::cerr << message << std::endl;
      // Possibly throw an exception or handle the error appropriately
  }

  void read_inp_file(const std::string& inpFileName) {
    std::cout << "\n** Read input file" << std::endl;

    std::ifstream inpFile(inpFileName);
    if (!inpFile.is_open()) {
        local_error("Unable to open file");
        return;
    }

    std::string line;
    int state = 0;
    while (getline(inpFile, line)) {
      // Trim leading and trailing whitespace
      line.erase(0, line.find_first_not_of(" \t"));
      line.erase(line.find_last_not_of(" \t") + 1);

      if (line.empty()) continue;
      if (line[0] == '*') {
        state = 0;
      }
      if (line.length() >= 5 && line.substr(0, 5).compare("*Node") == 0) {
          state = 1;
          continue;
      }
      if (line.length() >= 8 && line.substr(0, 8).compare("*Element") == 0) {
          state = 2;
          continue;
      }
      if (line.length() >= 9 && line.substr(0, 9).compare("*Boundary") == 0) {
          state = 3;
          continue;
      }

      //std::cout << "state: " << state << " : " << line << std::endl;

      if (state == 0) {
          continue;
      }

      std::istringstream iss(line);
      std::string token;
      std::vector<std::string> tokens;
      while (getline(iss, token, ',')) {
          tokens.push_back(token);
      }

      if (state == 1) {
        // Read nodes
        if (tokens.size() != 4) {
          local_error("A node definition needs 4 values");
          continue;
        }
        int nodeNr = std::stoi(tokens[0]) - 1; // zero indexed
        double xx = std::stod(tokens[1]) - 4.0;
        double yy = std::stod(tokens[2]) - 4.0;
        double zz = std::stod(tokens[3]);
        _nodes.push_back(Eigen::Vector3d(xx, yy, zz)); // assume the nodes are ordered 1, 2, 3...
      } else if (state == 2) {
        // Read elements
        if (tokens.size() != 7) {
          local_error("An element definition needs 7 values");
          continue;
        }
        int elemNr = std::stoi(tokens[0]); // not used in this context
        std::vector<unsigned int> element;
        for (int i = 1; i <= 6; ++i) {
          element.push_back(std::stoi(tokens[i]) - 1); // zero indexed
        }
        // Assuming the order needs to be adjusted from the Python code
        _elements.push_back({element[0], element[1], element[2], element[3], element[4], element[5]});
        //_elements.push_back({element[0], element[3], element[2], element[1]});
      } else if (state == 3) {
        std::cout << "BCs not implemented in file read" << std::endl;
        // Read displacement boundary conditions
        //if (tokens.size() != 4) {
        //  local_error("A displacement boundary condition needs 4 values");
        //  continue;
        //}
        //int nodeNr = std::stoi(tokens[0]) - 1; // zero indexed
        //int dof1 = std::stoi(tokens[1]);
        //int dof2 = std::stoi(tokens[2]);
        //double val = std::stod(tokens[3]);
        //if (dof1 == 1) {
        //  boundary.push_back({static_cast<double>(nodeNr), 1.0, val});
        //}
        //if (dof2 == 2) {
        //  boundary.push_back({static_cast<double>(nodeNr), 2.0, val});
        //}
      }
    }

    inpFile.close();

    std::cout << "Nodes: " << _nodes.size() << std::endl;
    std::cout << "Elements: " << _elements.size() << std::endl;
    //std::cout << "Boundary: " << boundary.size() << std::endl;

  }


  void ApplyDisplacements(std::vector<Eigen::Vector3d> &displacements, double scale = 1.0) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      _nodes_deformed.push_back(_nodes[n] + scale * displacements[n]);
      std::cout << "Node " << n << ": " << _nodes[n].transpose() << " -> " << _nodes_deformed[n].transpose() << std::endl;
    }
  }
    
  std::vector<Eigen::Vector3d> _nodes, _nodes_deformed;
  std::vector<std::vector<unsigned int>> _elements;

  int _num_elements = 10;
  double _x0 = 1.0;
  double _y0 = 1.0;
  double _z0 = 1.0;
  double _w = 1.0;

  std::string Name() { return "c3d6_3"; }
  double LoadLocation() { return _z0; }
  std::string LoadDirection() { return "z"; }
  int LoadDirectionKey() { return 3; }
  double LoadMagnitude() { return 1.0; }
  std::string ElementType() { return "C3D6"; }
};





class AbaqusC3D6_4 {
public:
  AbaqusC3D6_4(std::string filename) {
    read_inp_file(filename);
  }

  // Assuming `local_error` is a function that handles errors. You'll need to define it.
  void local_error(const std::string &message) {
      std::cerr << message << std::endl;
      // Possibly throw an exception or handle the error appropriately
  }

  void read_inp_file(const std::string& inpFileName) {
    std::cout << "\n** Read input file" << std::endl;

    std::ifstream inpFile(inpFileName);
    if (!inpFile.is_open()) {
        local_error("Unable to open file");
        return;
    }

    std::string line;
    int state = 0;
    while (getline(inpFile, line)) {
      // Trim leading and trailing whitespace
      line.erase(0, line.find_first_not_of(" \t"));
      line.erase(line.find_last_not_of(" \t") + 1);

      if (line.empty()) continue;
      if (line[0] == '*') {
        state = 0;
      }
      if (line.length() >= 5 && line.substr(0, 5).compare("*Node") == 0) {
          state = 1;
          continue;
      }
      if (line.length() >= 8 && line.substr(0, 8).compare("*Element") == 0) {
          state = 2;
          continue;
      }
      if (line.length() >= 9 && line.substr(0, 9).compare("*Boundary") == 0) {
          state = 3;
          continue;
      }

      //std::cout << "state: " << state << " : " << line << std::endl;

      if (state == 0) {
          continue;
      }

      std::istringstream iss(line);
      std::string token;
      std::vector<std::string> tokens;
      while (getline(iss, token, ',')) {
          tokens.push_back(token);
      }

      if (state == 1) {
        // Read nodes
        if (tokens.size() != 4) {
          local_error("A node definition needs 4 values");
          continue;
        }
        int nodeNr = std::stoi(tokens[0]) - 1; // zero indexed
        double xx = std::stod(tokens[1]);
        double yy = std::stod(tokens[2]);
        double zz = std::stod(tokens[3]);
        _nodes.push_back(Eigen::Vector3d(xx, yy, zz)); // assume the nodes are ordered 1, 2, 3...
      } else if (state == 2) {
        // Read elements
        if (tokens.size() != 7) {
          local_error("An element definition needs 7 values");
          continue;
        }
        int elemNr = std::stoi(tokens[0]); // not used in this context
        std::vector<unsigned int> element;
        for (int i = 1; i <= 6; ++i) {
          element.push_back(std::stoi(tokens[i]) - 1); // zero indexed
        }
        // Assuming the order needs to be adjusted from the Python code
        _elements.push_back({element[0], element[1], element[2], element[3], element[4], element[5]});
        //_elements.push_back({element[0], element[3], element[2], element[1]});
      } else if (state == 3) {
        std::cout << "BCs not implemented in file read" << std::endl;
        // Read displacement boundary conditions
        //if (tokens.size() != 4) {
        //  local_error("A displacement boundary condition needs 4 values");
        //  continue;
        //}
        //int nodeNr = std::stoi(tokens[0]) - 1; // zero indexed
        //int dof1 = std::stoi(tokens[1]);
        //int dof2 = std::stoi(tokens[2]);
        //double val = std::stod(tokens[3]);
        //if (dof1 == 1) {
        //  boundary.push_back({static_cast<double>(nodeNr), 1.0, val});
        //}
        //if (dof2 == 2) {
        //  boundary.push_back({static_cast<double>(nodeNr), 2.0, val});
        //}
      }
    }

    inpFile.close();

    std::cout << "Nodes: " << _nodes.size() << std::endl;
    std::cout << "Elements: " << _elements.size() << std::endl;
    //std::cout << "Boundary: " << boundary.size() << std::endl;

  }


  void ApplyDisplacements(std::vector<Eigen::Vector3d> &displacements, double scale = 1.0) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      _nodes_deformed.push_back(_nodes[n] + scale * displacements[n]);
      std::cout << "Node " << n << ": " << _nodes[n].transpose() << " -> " << _nodes_deformed[n].transpose() << std::endl;
    }
  }
    
  std::vector<Eigen::Vector3d> _nodes, _nodes_deformed;
  std::vector<std::vector<unsigned int>> _elements;

  int _num_elements = 10;
  double _x0 = 1.0;
  double _y0 = 1.0;
  double _z0 = 2.0;
  double _w = 1.0;

  std::string Name() { return "c3d6_4"; }
  double LoadLocation() { return _z0; }
  std::string LoadDirection() { return "z"; }
  int LoadDirectionKey() { return 3; }
  double LoadMagnitude() { return 1.0; }
  std::string ElementType() { return "C3D6"; }
};




/*
Eigen::MatrixXd read_stiffness_matrix_abaqus(std::string inpFileName) {
    std::cout << "\n** Read input file" << std::endl;

    std::ifstream inpFile(inpFileName);
    if (!inpFile.is_open()) {
        std::cout << "ERROR: Unable to open file" << std::endl;
        return Eigen::MatrixXd(0,0);
    }

    std::string line;
    int state = 0;
    while (getline(inpFile, line)) {
      // Trim leading and trailing whitespace
      line.erase(0, line.find_first_not_of(" \t"));
      line.erase(line.find_last_not_of(" \t") + 1);

      if (line.empty()) continue;
      if (line[0] == '*') {
        state = 0;
      }
      if (line.length() >= 5 && line.substr(0, 5).compare("*Node") == 0) {
          state = 1;
          continue;
      }
      if (line.length() >= 8 && line.substr(0, 8).compare("*Element") == 0) {
          state = 2;
          continue;
      }
      if (line.length() >= 9 && line.substr(0, 9).compare("*Boundary") == 0) {
          state = 3;
          continue;
      }

      //std::cout << "state: " << state << " : " << line << std::endl;

      if (state == 0) {
          continue;
      }

      std::istringstream iss(line);
      std::string token;
      std::vector<std::string> tokens;
      while (getline(iss, token, ',')) {
          tokens.push_back(token);
      }

      if (state == 1) {
        // Read nodes
        if (tokens.size() != 4) {
          local_error("A node definition needs 4 values");
          continue;
        }
        int nodeNr = std::stoi(tokens[0]) - 1; // zero indexed
        double xx = std::stod(tokens[1]);
        double yy = std::stod(tokens[2]);
        double zz = std::stod(tokens[3]);
        _nodes.push_back(Eigen::Vector3d(xx, yy, zz)); // assume the nodes are ordered 1, 2, 3...
      } else if (state == 2) {
        // Read elements
        if (tokens.size() != 7) {
          local_error("An element definition needs 7 values");
          continue;
        }
        int elemNr = std::stoi(tokens[0]); // not used in this context
        std::vector<unsigned int> element;
        for (int i = 1; i <= 6; ++i) {
          element.push_back(std::stoi(tokens[i]) - 1); // zero indexed
        }
        // Assuming the order needs to be adjusted from the Python code
        _elements.push_back({element[0], element[1], element[2], element[3], element[4], element[5]});
        //_elements.push_back({element[0], element[3], element[2], element[1]});
      } else if (state == 3) {
        std::cout << "BCs not implemented in file read" << std::endl;
        // Read displacement boundary conditions
        //if (tokens.size() != 4) {
        //  local_error("A displacement boundary condition needs 4 values");
        //  continue;
        //}
        //int nodeNr = std::stoi(tokens[0]) - 1; // zero indexed
        //int dof1 = std::stoi(tokens[1]);
        //int dof2 = std::stoi(tokens[2]);
        //double val = std::stod(tokens[3]);
        //if (dof1 == 1) {
        //  boundary.push_back({static_cast<double>(nodeNr), 1.0, val});
        //}
        //if (dof2 == 2) {
        //  boundary.push_back({static_cast<double>(nodeNr), 2.0, val});
        //}
      }
    }

    inpFile.close();

    std::cout << "Nodes: " << _nodes.size() << std::endl;
    std::cout << "Elements: " << _elements.size() << std::endl;
    //std::cout << "Boundary: " << boundary.size() << std::endl;

  }
}
*/