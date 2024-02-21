#ifndef DATA_HPP
#define DATA_HPP

#include <iostream>
#include <vector>
#include <string>

#include <fstream>
#include <iterator>
#include <algorithm>
#include <boost/algorithm/string.hpp>


// Load data from csv file with custom delimiter.
std::vector<std::vector<float> > get_from_file_vvf (std::string inputPath, std::string delimiter);

std::vector<std::vector<int> > get_from_file_vvn (std::string inputPath, std::string delimiter);

void put_to_file_vvf (std::string outputPath, std::string delimiter, std::vector<std::vector<float> > vvoutput, bool append);

std::vector<std::vector<float> > vector_resize_cols(std::vector<std::vector<float> > v1, int n);


class Dataset {

public:
  Dataset(std::string data_path, std::string element_type);

  std::vector<std::vector<float>> points();
  std::vector<std::vector<int>> elements();
  std::vector<std::vector<float>> forces();
  std::vector<std::vector<int>> dirichlet();



private:
  std::vector<std::vector<float>> vpts_;
  std::vector<std::vector<int>> vElts_;
  std::vector<std::vector<float>> vvF_;
  std::vector<std::vector<int>> vDir_;

  std::string path_;
  std::string element_type_;

};


#endif