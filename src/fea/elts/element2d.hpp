#ifndef ELEMENT2D_HPP
#define ELEMENT2D_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>

#include "element.hpp"


class Element2D : public Element {
public:
  Element2D(double E, double nu) : Element(E, nu) {
    computeElasticityMatrix();
  }

  void computeElasticityMatrix() override;

  // Function to compute the Jacobian matrix for a triangular prism
  Eigen::MatrixXd computeJacobian(const std::vector<Eigen::Vector3d>& nodes, Eigen::MatrixXd &dN);

  // Function to compute the inverse of the Jacobian matrix and its determinant
  std::pair<Eigen::MatrixXd, double> computeInverseJacobianAndDet(const Eigen::MatrixXd& J);

  // Function to compute the stiffness matrix for a triangular prism
  Eigen::MatrixXd computeStiffnessMatrix(const std::vector<Eigen::Vector3d>& nodes);

  // Function to assemble the global stiffness matrix
  Eigen::MatrixXd matAssembly(std::vector<Eigen::Vector3d> &vpts, 
                              std::vector<std::vector<unsigned int>> &velts);

};


#endif // ELEMENT2D_HPP