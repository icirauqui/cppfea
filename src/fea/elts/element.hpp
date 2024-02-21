#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>


class FEAData {

public:
  std::vector<std::vector<double>> node_stress;
  std::vector<std::vector<double>> node_strain;
  std::vector<double> emin;
  std::vector<double> emax;
  std::vector<double> smin;
  std::vector<double> smax;
  std::vector<double> umin;
  std::vector<double> umax;
  double strain_energy;
  int number_of_nodes;
  int number_of_elements;
  std::string element_name;
  //std::string analysis_type;
  //std::string material_model;
  //std::string element_type;
  //std::string element_order;
}; // FEAData



class Element {
public:
  Element(double E, double nu) {
    _E = E;
    _nu = nu;
  }

  // Function to compute shape functions for a triangular prism
  virtual Eigen::VectorXd computeShapeFunctions(double xi, double eta, double zeta) = 0;

  // Function to compute the derivative of shape functions for a triangular prism
  virtual Eigen::MatrixXd computeShapeFunctionDerivatives(double xi, double eta, double zeta) = 0;

  // Function to compute the Jacobian matrix for a triangular prism
  virtual Eigen::MatrixXd computeJacobian(const std::vector<Eigen::Vector3d>& nodes, Eigen::MatrixXd &dN) = 0;

  // Function to compute the inverse of the Jacobian matrix and its determinant
  virtual std::pair<Eigen::MatrixXd, double> computeInverseJacobianAndDet(const Eigen::MatrixXd& J) = 0;

  // Function to compute the Strain-Displacement Matrix (B)
  virtual Eigen::MatrixXd computeStrainDisplacementMatrix(Eigen::MatrixXd &dN, Eigen::MatrixXd &invJ, double &detJ) = 0;

  // Function to compute the elasticity matrix
  virtual void computeElasticityMatrix() = 0;

  // Function to compute the stiffness matrix for a triangular prism
  virtual Eigen::MatrixXd computeStiffnessMatrix(const std::vector<Eigen::Vector3d>& nodes) = 0;

  // Function to assemble the global stiffness matrix
  virtual Eigen::MatrixXd matAssembly(std::vector<Eigen::Vector3d> &vpts, 
                                      std::vector<std::vector<unsigned int>> &velts) = 0;

  virtual void postProcess(std::vector<Eigen::Vector3d> &vpts, 
                           std::vector<std::vector<unsigned int>> &velts,
                           Eigen::MatrixXd U,
                           FEAData &data);


  // Accessors
  std::string getElementName() const { return _element_name; }
  int getNumNodes() const { return _num_nodes; }
  int getDofPerNode() const { return _dof_per_node; }
  double getE() const { return _E; }
  double getNu() const { return _nu; }
  Eigen::MatrixXd getD() const { return _D; }

protected:
  std::string _element_name;
  int _num_nodes;
  int _dof_per_node;

  double _E;
  double _nu;
  Eigen::MatrixXd _D;
};


#endif // ELEMENT_HPP