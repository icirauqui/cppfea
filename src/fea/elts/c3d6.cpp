#include "c3d6.hpp"


C3D6::C3D6(double E, double nu) : Element3D(E, nu) {
  _element_name = "C3D6";
  _num_nodes = 6;
  _dof_per_node = 3;
}


Eigen::VectorXd C3D6::computeShapeFunctions(double xi, double eta, double zeta) {
    // Initialize a vector to store the shape functions
    Eigen::VectorXd N(_num_nodes);

    // Compute the shape functions for each node
    N(0) = (1 - xi - eta) * (1 - zeta) / 2;
    N(1) =             xi * (1 - zeta) / 2;
    N(2) =            eta * (1 - zeta) / 2;
    N(3) = (1 - xi - eta) * (1 + zeta) / 2;
    N(4) =             xi * (1 + zeta) / 2;
    N(5) =            eta * (1 + zeta) / 2;

    return N;
}


Eigen::MatrixXd C3D6::computeShapeFunctionDerivatives(double xi, double eta, double zeta) {
    Eigen::MatrixXd dN(6, 3); // 6 nodes, 3 natural coordinates (xi, eta, zeta)

    // Derivatives of shape functions with respect to xi, eta, zeta
    dN(0, 0) = -(1 - zeta) / 2; // dN1/dxi
    dN(0, 1) = -(1 - zeta) / 2; // dN1/deta
    dN(0, 2) = -(1 - xi - eta) / 2; // dN1/dzeta

    dN(1, 0) = (1 - zeta) / 2;  // dN2/dxi
    dN(1, 1) = 0;               // dN2/deta
    dN(1, 2) = -xi / 2;         // dN2/dzeta

    dN(2, 0) = 0;               // dN3/dxi
    dN(2, 1) = (1 - zeta) / 2;  // dN3/deta
    dN(2, 2) = -eta / 2;        // dN3/dzeta

    dN(3, 0) = -(1 + zeta) / 2; // dN4/dxi
    dN(3, 1) = -(1 + zeta) / 2; // dN4/deta
    dN(3, 2) = (1 - xi - eta) / 2; // dN4/dzeta

    dN(4, 0) = (1 + zeta) / 2;  // dN5/dxi
    dN(4, 1) = 0;               // dN5/deta
    dN(4, 2) = xi / 2;          // dN5/dzeta

    dN(5, 0) = 0;               // dN6/dxi
    dN(5, 1) = (1 + zeta) / 2;  // dN6/deta
    dN(5, 2) = eta / 2;         // dN6/dzeta

    return dN;
}


// Function to compute the Strain-Displacement Matrix (B)
Eigen::MatrixXd C3D6::computeStrainDisplacementMatrix(Eigen::MatrixXd &dN, Eigen::MatrixXd &invJ, double &detJ) {
    // Compute the derivatives of shape functions w.r.t. global coordinates
    Eigen::MatrixXd dNdXYZ = invJ * dN.transpose();

    //std::cout << " dNdXYZ " << std::endl << dNdXYZ << std::endl;

    // Initialize the B matrix
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 18); // 6 strain components, 18 displacement components (3 per node)

    // Fill the B matrix
    for (int i = 0; i < 6; ++i) {
        B(0, i * 3)     = dNdXYZ(0, i); // Strain εxx
        B(1, i * 3 + 1) = dNdXYZ(1, i); // Strain εyy
        B(2, i * 3 + 2) = dNdXYZ(2, i); // Strain εzz

        //B(3, i * 3 + 0) = dNdXYZ(1, i); // Shear γyx
        //B(3, i * 3 + 1) = dNdXYZ(0, i);
        //
        //B(4, i * 3 + 1) = dNdXYZ(2, i); // Shear γzy
        //B(4, i * 3 + 2) = dNdXYZ(1, i);
        //
        //B(5, i * 3 + 0) = dNdXYZ(2, i); // Shear γzx
        //B(5, i * 3 + 2) = dNdXYZ(0, i);

        B(3, i * 3 + 0) = dNdXYZ(1, i); // Shear γyx
        B(3, i * 3 + 1) = dNdXYZ(0, i);
        
        B(4, i * 3 + 0) = dNdXYZ(2, i); // Shear γzx
        B(4, i * 3 + 2) = dNdXYZ(0, i);
        
        B(5, i * 3 + 1) = dNdXYZ(2, i); // Shear γzy
        B(5, i * 3 + 2) = dNdXYZ(1, i);
    }

    return B;
}