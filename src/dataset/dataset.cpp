#include "dataset.hpp"

// Load data from csv file with custom delimiter.
std::vector<std::vector<float> > get_from_file_vvf (std::string inputPath, std::string delimiter) {
    std::ifstream file(inputPath);
    std::vector<std::vector<std::string> > datalist;
    std::string line = "";

    while (getline(file,line)) {
        std::vector<std::string> vec;
        boost::algorithm::split(vec,line,boost::is_any_of(delimiter));
        datalist.push_back(vec);
    }
    file.close();

    std::vector<std::vector<float> > vPoints;
    for (unsigned int i=0; i<datalist.size(); i++) {
        std::vector<float> vPoint;
        for (unsigned int j=0; j<datalist[i].size(); j++) {
            std::string sPoint = datalist[i][j];
            float fPoint = strtof((sPoint).c_str(),0);
            vPoint.push_back(fPoint);
        }
        vPoints.push_back(vPoint);
    }

    return vPoints;
}


std::vector<std::vector<int> > get_from_file_vvn (std::string inputPath, std::string delimiter) {
    std::ifstream file(inputPath);
    std::vector<std::vector<std::string> > datalist;
    std::string line = "";

    while (getline(file,line)) {
        std::vector<std::string> vec;
        boost::algorithm::split(vec,line,boost::is_any_of(delimiter));
        datalist.push_back(vec);
    }
    file.close();

    std::vector<std::vector<int> > vPoints;
    for (unsigned int i=0; i<datalist.size(); i++) {
        std::vector<int> vPoint;
        for (unsigned int j=0; j<datalist[i].size(); j++) {
            std::string sPoint = datalist[i][j];
            int nPoint = std::stoi((sPoint).c_str(),0);
            vPoint.push_back(nPoint);
        }
        vPoints.push_back(vPoint);
    }

    return vPoints;
}


void put_to_file_vvf (std::string outputPath, std::string delimiter, std::vector<std::vector<float> > vvoutput, bool append){
    
    std::ofstream os;
    if (append)
        os.open(outputPath.c_str(), std::ios::out | std::ios::app );
    else
        os.open(outputPath.c_str(), std::ios::out);
    
    for (unsigned int i=0; i<vvoutput.size(); i++) {
        for (unsigned int j=0; j<vvoutput[i].size(); j++) {
            os << vvoutput[i][j] << " ";
        }
        os << std::endl;
    }
    os << std::endl;
    
    os.close();
}


std::vector<std::vector<float> > vector_resize_cols(std::vector<std::vector<float> > v1, int n){
    std::vector<std::vector<float> > v2;
    std::vector<float> v2i;
    for (unsigned int i=0; i<v1.size(); i++){
        for (unsigned int j=0; j<v1[i].size(); j++){
            v2i.push_back(v1[i][j]);
            if (v2i.size() == n){
                v2.push_back(v2i);
                v2i.clear();
            }
        }
    }
    return v2;
}



Dataset::Dataset(std::string data_path, std::string element_type) {
  path_ = data_path;
  element_type_ = element_type;

  vpts_ = get_from_file_vvf(path_ + "/input_" + element_type_ + "/input_points.csv",",");
  vElts_ = get_from_file_vvn(path_ + "/input_" + element_type_ + "/input_elements.csv",",");
  vvF_ = get_from_file_vvf(path_ + "/input_" + element_type_ + "/input_forces.csv",",");
  vDir_ = get_from_file_vvn(path_ + "/input_" + element_type_ + "/input_Dirichlet.csv",",");

  // Node numbers in Abaqus start at 1, change to 0 for cpp
  for (unsigned int i=0; i<vElts_.size(); i++){
      for (unsigned int j=0; j<vElts_[i].size(); j++){
          vElts_[i][j] -= 1;
      }
  }
}

std::vector<std::vector<float>> Dataset::points() {
  return vpts_;
}

std::vector<std::vector<int>> Dataset::elements() {
  return vElts_;
}

std::vector<std::vector<float>> Dataset::forces() {
  return vvF_;
}

std::vector<std::vector<int>> Dataset::dirichlet() {
  return vDir_;
}