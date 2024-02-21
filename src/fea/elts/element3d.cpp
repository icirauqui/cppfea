#include "element3d.hpp"



void Element3D::computeElasticityMatrix() {
  double lambda = _E * _nu / ((1 + _nu) * (1 - 2 * _nu));
  double G = _E / (2 * (1 + _nu));

  _D = Eigen::MatrixXd::Zero(6, 6);
  _D(0, 0) = _D(1, 1) = _D(2, 2) = lambda + 2 * G;
  _D(0, 1) = _D(0, 2) = _D(1, 0) = _D(1, 2) = _D(2, 0) = _D(2, 1) = lambda;
  _D(3, 3) = _D(4, 4) = _D(5, 5) = G;
}

Eigen::MatrixXd Element3D::computeJacobian(const std::vector<Eigen::Vector3d>& nodes, Eigen::MatrixXd &dN) {
    // Initialize the Jacobian matrix
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, 3);

    // Compute the Jacobian matrix
    for (int i = 0; i < _num_nodes; ++i) {
        J(0, 0) += nodes[i](0) * dN(i, 0);
        J(1, 0) += nodes[i](1) * dN(i, 0);
        J(2, 0) += nodes[i](2) * dN(i, 0);
        J(0, 1) += nodes[i](0) * dN(i, 1);
        J(1, 1) += nodes[i](1) * dN(i, 1);
        J(2, 1) += nodes[i](2) * dN(i, 1);
        J(0, 2) += nodes[i](0) * dN(i, 2);
        J(1, 2) += nodes[i](1) * dN(i, 2);
        J(2, 2) += nodes[i](2) * dN(i, 2);


        //J.col(0) += nodes[i] * dN(i, 0); // Contribution to the first column of J
        //J.col(1) += nodes[i] * dN(i, 1); // Contribution to the second column of J
        //J.col(2) += nodes[i] * dN(i, 2); // Contribution to the third column of J
    }

    return J;
}

std::pair<Eigen::MatrixXd, double> Element3D::computeInverseJacobianAndDet(const Eigen::MatrixXd& J) {
  double detJ = J.determinant();
  Eigen::MatrixXd invJ = J.inverse();
  return {invJ, detJ};
}

Eigen::MatrixXd Element3D::computeStiffnessMatrix(const std::vector<Eigen::Vector3d>& nodes) {
    // Initialize the stiffness matrix: 18x18 for 6 nodes, 3 DOF each
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(_num_nodes*_dof_per_node, _num_nodes*_dof_per_node);

    // Gauss quadrature points and weights (2-point quadrature)
    std::array<double, 2> gaussPoints = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    std::array<double, 2> gaussWeights = {1.0, 1.0};

    // Integration (simplified example, replace with appropriate Gauss quadrature in real applications)
    for (int i = 0; i < gaussPoints.size(); ++i) {
        for (int j = 0; j < gaussPoints.size(); ++j) {
            for (int k = 0; k < gaussPoints.size(); ++k) {
                double xi = gaussPoints[i];
                double eta = gaussPoints[j];
                double zeta = gaussPoints[k];

                Eigen::MatrixXd dN = computeShapeFunctionDerivatives(xi, eta, zeta);

                Eigen::MatrixXd J = computeJacobian(nodes, dN);
                auto [invJ, detJ] = computeInverseJacobianAndDet(J);
                Eigen::MatrixXd B = computeStrainDisplacementMatrix(dN, invJ, detJ);

                std::cout << std::endl << std::endl;

                // Weight calculation considering different weights
                double weight = gaussWeights[i] * gaussWeights[j] * gaussWeights[k];

                // K += B^T * D * B * detJ
                K += B.transpose() * _D * B * detJ;
            }
        }
    }

    return K;
}

Eigen::MatrixXd Element3D::matAssembly(std::vector<Eigen::Vector3d> &vpts, 
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