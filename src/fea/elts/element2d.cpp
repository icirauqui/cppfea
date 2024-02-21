#include "element2d.hpp"


void Element2D::computeElasticityMatrix() {

  //double mult = _E/(1.0+_nu)/(1.0-2.0*_nu);
  //_D = Eigen::MatrixXd::Zero(3, 3);
  //_D(0, 0) = _D(1, 1) = (1-_nu);
  //_D(0, 1) = _D(1, 0) = _nu;
  //_D(2, 2) = 0.5 - _nu;
  //_D *= mult;

  //std::cout << "D: " << std::endl << _D << std::endl;

  double L = _E * _nu / ((1 + _nu) * (1 - 2 * _nu));
  double G = _E / (2 * (1 + _nu));
    
  double Ls = 2*L*G / (L + 2*G);
  double Gs = G;

  Ls = L;
    
  _D = Eigen::MatrixXd::Zero(3, 3);
  _D(0, 0) = _D(1, 1) = Ls + 2 * Gs;
  _D(2, 2) = Gs;
  _D(0, 1) = _D(1, 0) = Ls;

  std::cout << "D: " << std::endl << _D << std::endl;
}

Eigen::MatrixXd Element2D::computeJacobian(const std::vector<Eigen::Vector3d>& nodes, Eigen::MatrixXd &dN) {
  // Initialize the Jacobian matrix
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(2, 2);

  // Compute the Jacobian matrix
  for (int i = 0; i < 4; ++i) {
    Eigen::Vector2d node = nodes[i].head(2);
    J.col(0) += node * dN(i, 0); // Contribution to the first column of J
    J.col(1) += node * dN(i, 1); // Contribution to the second column of J
  }

  return J;
}

// Function to compute the inverse of the Jacobian matrix and its determinant
std::pair<Eigen::MatrixXd, double> Element2D::computeInverseJacobianAndDet(const Eigen::MatrixXd& J) {
  double detJ = J.determinant();
  Eigen::MatrixXd invJ = J.inverse();
  return {invJ, detJ};
}


Eigen::MatrixXd Element2D::computeStiffnessMatrix(const std::vector<Eigen::Vector3d>& nodes) {
  // Initialize the stiffness matrix: 24x24 for 4 nodes, 2 DOF each
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(_num_nodes*_dof_per_node, _num_nodes*_dof_per_node);

  // Gauss quadrature points and weights (2-point quadrature)
  std::array<double, 2> gaussPoints = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
  std::array<double, 2> gaussWeights = {1.0, 1.0};

  // Integration (simplified example, replace with appropriate Gauss quadrature in real applications)
  for (int i = 0; i < gaussPoints.size(); ++i) {
    for (int j = 0; j < gaussPoints.size(); ++j) {
      double xi = gaussPoints[i];
      double eta = gaussPoints[j];
      double zeta = 0.0;

      Eigen::MatrixXd dN = computeShapeFunctionDerivatives(xi, eta, zeta);
                
      Eigen::MatrixXd J = computeJacobian(nodes, dN);
      auto [invJ, detJ] = computeInverseJacobianAndDet(J);
      Eigen::MatrixXd B = computeStrainDisplacementMatrix(dN, invJ, detJ);

      // Weight calculation considering different weights
      double weight = gaussWeights[i] * gaussWeights[j];

      // K += B^T * D * B * detJ
      K += B.transpose() * _D * B * detJ;
    }
  }

  return K;
}


Eigen::MatrixXd Element2D::matAssembly(std::vector<Eigen::Vector3d> &vpts, 
                                  std::vector<std::vector<unsigned int>> &velts) {
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(_dof_per_node*vpts.size(), _dof_per_node*vpts.size());

  for (auto elt : velts) {

    std::vector<Eigen::Vector3d> xyzi(_num_nodes);
    std::vector<int> mn(_num_nodes);

    for (unsigned int i=0; i<elt.size(); i++) {
      xyzi[i] = vpts[elt[i]];
      mn[i] = elt[i]*_dof_per_node;
    }

    Eigen::MatrixXd Kei = computeStiffnessMatrix(xyzi);

    for (unsigned int ni = 0; ni < mn.size(); ni++) {
      for (unsigned int nj = 0; nj < mn.size(); nj++) {
        for (unsigned int m = 0; m < _dof_per_node; m++) {
          for (unsigned int n = 0; n < _dof_per_node; n++) {
            K(mn[ni]+m, mn[nj]+n) += Kei(ni*_dof_per_node+m, nj*_dof_per_node+n);
          }
        }
      }
    }
  }

  return K;
}