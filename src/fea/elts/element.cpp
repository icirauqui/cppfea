#include "element.hpp"

void Element::postProcess(std::vector<Eigen::Vector3d> &vpts, 
                           std::vector<std::vector<unsigned int>> &velts,
                           Eigen::MatrixXd U,
                           FEAData &data) {

  data.node_stress = std::vector<std::vector<double>>(vpts.size(), std::vector<double>(3, 0.0));
  data.node_strain = std::vector<std::vector<double>>(vpts.size(), std::vector<double>(3, 0.0));
  data.emin = std::vector<double>(_D.rows(), 0.0);
  data.emax = std::vector<double>(_D.rows(), 0.0);
  data.smin = std::vector<double>(_D.rows(), 0.0);
  data.smax = std::vector<double>(_D.rows(), 0.0);
  data.umin = std::vector<double>(_dof_per_node, 0.0);
  data.umax = std::vector<double>(_dof_per_node, 0.0);

  // Gauss quadrature points and weights (2-point quadrature)
  std::array<double, 2> gaussPoints = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
  std::array<double, 2> gaussWeights = {1.0, 1.0};

  for (auto elt : velts) {

    std::vector<Eigen::Vector3d> xyzi(_num_nodes);
    Eigen::MatrixXd xyzi_u(_num_nodes*_dof_per_node, 1);
    std::vector<int> mn(_num_nodes);

    for (unsigned int i=0; i<elt.size(); i++) {
      xyzi[i] = vpts[elt[i]];
      mn[i] = elt[i]*_dof_per_node;
      for (unsigned int j=0; j<_dof_per_node; j++) {
        xyzi_u(i*_dof_per_node+j, 0) = U(mn[i]+j, 0);
      }
    }

    for (int i = 0; i < gaussPoints.size(); ++i) {
      for (int j = 0; j < gaussPoints.size(); ++j) {
        double xi = gaussPoints[i];
        double eta = gaussPoints[j];
        double zeta = 0.0;

        Eigen::MatrixXd dN = computeShapeFunctionDerivatives(xi, eta, zeta);
        Eigen::MatrixXd J = computeJacobian(xyzi, dN);
        auto [invJ, detJ] = computeInverseJacobianAndDet(J);
        Eigen::MatrixXd B = computeStrainDisplacementMatrix(dN, invJ, detJ);

        Eigen::MatrixXd strain = B * xyzi_u;
        Eigen::MatrixXd stress = _D * strain;

        for (unsigned int ni = 0; ni < _num_nodes; ni++) {
          for (unsigned int m = 0; m < _dof_per_node; m++) {
            data.node_stress[elt[ni]][m] += stress(m, 0);
            data.node_strain[elt[ni]][m] += strain(m, 0);
          }
        }

        for (unsigned int m = 0; m < _D.rows(); m++) {
          data.emin[m] = std::min(data.emin[m], strain(m, 0));
          data.emax[m] = std::max(data.emax[m], strain(m, 0));
          data.smin[m] = std::min(data.smin[m], stress(m, 0));
          data.smax[m] = std::max(data.smax[m], stress(m, 0));
        }
      }
    }
  }

  for (unsigned int i = 0; i < vpts.size(); i++) {
    for (unsigned int j = 0; j < _dof_per_node; j++) {
      data.umin[j] = std::min(data.umin[j], U(i*_dof_per_node+j, 0));
      data.umax[j] = std::max(data.umax[j], U(i*_dof_per_node+j, 0));
    }
  }

}

