#include "c3d8.hpp"


C3D8::C3D8(double E, double nu) : Element3D(E, nu) {
  _element_name = "C3D8";
  _num_nodes = 8;
  _dof_per_node = 3;
}


Eigen::VectorXd C3D8::computeShapeFunctions(double xi, double eta, double zeta) {
    // Initialize a vector to store the shape functions
    Eigen::VectorXd N(_num_nodes);

    // Compute the shape functions for each node
    N(0) = (1 - xi) * (1 - eta) * (1 - zeta) / 8;
    N(1) = (1 + xi) * (1 - eta) * (1 - zeta) / 8;
    N(2) = (1 + xi) * (1 + eta) * (1 - zeta) / 8;
    N(3) = (1 - xi) * (1 + eta) * (1 - zeta) / 8;
    N(4) = (1 - xi) * (1 - eta) * (1 + zeta) / 8;
    N(5) = (1 + xi) * (1 - eta) * (1 + zeta) / 8;
    N(6) = (1 + xi) * (1 + eta) * (1 + zeta) / 8;
    N(7) = (1 - xi) * (1 + eta) * (1 + zeta) / 8;    

    return N;
}


Eigen::MatrixXd C3D8::computeShapeFunctionDerivatives(double xi, double eta, double zeta) {
    Eigen::MatrixXd dN(8, 3); // 6 nodes, 3 natural coordinates (xi, eta, zeta)

    // Derivatives of shape functions with respect to xi, eta, zeta
    dN(0, 0) = -(1 - eta) * (1 - zeta) / 8; // dN1/dxi
    dN(0, 1) = -(1 - xi)  * (1 - zeta) / 8; // dN1/deta
    dN(0, 2) = -(1 - xi)  * (1 - eta)  / 8; // dN1/dzeta

    dN(1, 0) =  (1 - eta) * (1 - zeta) / 8; // dN2/dxi
    dN(1, 1) = -(1 + xi)  * (1 - zeta) / 8; // dN2/deta
    dN(1, 2) = -(1 + xi)  * (1 - eta)  / 8; // dN2/dzeta

    dN(2, 0) =  (1 + eta) * (1 - zeta) / 8; // dN3/dxi
    dN(2, 1) =  (1 + xi)  * (1 - zeta) / 8; // dN3/deta
    dN(2, 2) = -(1 + xi)  * (1 + eta)  / 8; // dN3/dzeta

    dN(3, 0) = -(1 + eta) * (1 - zeta) / 8; // dN4/dxi
    dN(3, 1) =  (1 - xi)  * (1 - zeta) / 8; // dN4/deta
    dN(3, 2) = -(1 - xi)  * (1 + eta)  / 8; // dN4/dzeta

    dN(4, 0) = -(1 - eta) * (1 + zeta) / 8; // dN5/dxi
    dN(4, 1) = -(1 - xi)  * (1 + zeta) / 8; // dN5/deta
    dN(4, 2) =  (1 - xi)  * (1 - eta)  / 8; // dN5/dzeta

    dN(5, 0) =  (1 - eta) * (1 + zeta) / 8; // dN6/dxi
    dN(5, 1) = -(1 + xi)  * (1 + zeta) / 8; // dN6/deta
    dN(5, 2) =  (1 + xi)  * (1 - eta)  / 8; // dN6/dzeta

    dN(6, 0) =  (1 + eta) * (1 + zeta) / 8; // dN7/dxi
    dN(6, 1) =  (1 + xi)  * (1 + zeta) / 8; // dN7/deta
    dN(6, 2) =  (1 + xi)  * (1 + eta)  / 8; // dN7/dzeta

    dN(7, 0) = -(1 + eta) * (1 + zeta) / 8; // dN8/dxi
    dN(7, 1) =  (1 - xi)  * (1 + zeta) / 8; // dN8/deta
    dN(7, 2) =  (1 - xi)  * (1 + eta)  / 8; // dN8/dzeta

    return dN;
}


// Function to compute the Strain-Displacement Matrix (B)
Eigen::MatrixXd C3D8::computeStrainDisplacementMatrix(Eigen::MatrixXd &dN, Eigen::MatrixXd &invJ, double &detJ) {
    // Compute the derivatives of shape functions w.r.t. global coordinates
    Eigen::MatrixXd dNdXYZ = invJ * dN.transpose();

    // Initialize the B matrix
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 24); // 6 strain components, 24 displacement components (3 per node)

    // Fill the B matrix
    for (int i = 0; i < _num_nodes; ++i) {
        B(0, 3 * i) = dNdXYZ(0, i);
        B(1, 3 * i + 1) = dNdXYZ(1, i);
        B(2, 3 * i + 2) = dNdXYZ(2, i);

        B(3, 3 * i) = dNdXYZ(1, i);
        B(3, 3 * i + 1) = dNdXYZ(0, i);

        B(4, 3 * i + 1) = dNdXYZ(2, i);
        B(4, 3 * i + 2) = dNdXYZ(1, i);

        B(5, 3 * i) = dNdXYZ(2, i);
        B(5, 3 * i + 2) = dNdXYZ(0, i);
    }

    return B;
}
