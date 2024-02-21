#include "c2d4.hpp"


C2D4::C2D4(double E, double nu) : Element2D(E, nu) {
  _element_name = "C2D4";
  _num_nodes = 4;
  _dof_per_node = 2;
}

Eigen::VectorXd C2D4::computeShapeFunctions(double xi, double eta, double zeta) {
    // Initialize a vector to store the shape functions
    Eigen::VectorXd N(_num_nodes);

    // Compute the shape functions for each node
    N(0) = (1 - xi) * (1 - eta) / 4;
    N(1) = (1 + xi) * (1 - eta) / 4;
    N(2) = (1 + xi) * (1 + eta) / 4;
    N(3) = (1 - xi) * (1 + eta) / 4; 

    return N;
}

Eigen::MatrixXd C2D4::computeShapeFunctionDerivatives(double xi, double eta, double zeta) {
    Eigen::MatrixXd dN(4, 2); // 4 nodes, 2 natural coordinates (xi, eta)

    // Derivatives of shape functions with respect to xi, eta, zeta
    dN(0, 0) = -(1 - eta) / 4; // dN1/dxi
    dN(0, 1) = -(1 - xi)  / 4; // dN1/deta

    dN(1, 0) =  (1 - eta) / 4; // dN2/dxi
    dN(1, 1) = -(1 + xi)  / 4; // dN2/deta

    dN(2, 0) =  (1 + eta) / 4; // dN3/dxi
    dN(2, 1) =  (1 + xi)  / 4; // dN3/deta

    dN(3, 0) = -(1 + eta) / 4; // dN4/dxi
    dN(3, 1) =  (1 - xi)  / 4; // dN4/deta

    return dN;
}


// Function to compute the Strain-Displacement Matrix (B)
Eigen::MatrixXd C2D4::computeStrainDisplacementMatrix(Eigen::MatrixXd &dN, Eigen::MatrixXd &invJ, double &detJ) {
    // Compute the derivatives of shape functions w.r.t. global coordinates
    Eigen::MatrixXd dNdXYZ = invJ * dN.transpose();

    // Initialize the B matrix
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8); // 3 strain components, 8 displacement components (2 per node)

    // Fill the B matrix
    for (int i = 0; i < _num_nodes; ++i) {
        B(0, _dof_per_node * i)     = dNdXYZ(0, i);
        B(1, _dof_per_node * i + 1) = dNdXYZ(1, i);

        B(2, _dof_per_node * i)     = dNdXYZ(1, i);
        B(2, _dof_per_node * i + 1) = dNdXYZ(0, i);
    }

    return B;
}


