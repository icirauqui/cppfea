#ifndef C3D8_HPP
#define C3D8_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>

#include "element3d.hpp"

class C3D8 : public Element3D {

public:

    C3D8(double E, double nu);

    // Function to compute shape functions for a triangular prism
    Eigen::VectorXd computeShapeFunctions(double xi, double eta, double zeta);

    // Function to compute the derivative of shape functions for a triangular prism
    Eigen::MatrixXd computeShapeFunctionDerivatives(double xi, double eta, double zeta);

    // Function to compute the Strain-Displacement Matrix (B)
    Eigen::MatrixXd computeStrainDisplacementMatrix(Eigen::MatrixXd &dN, Eigen::MatrixXd &invJ, double &detJ);
}; // class c3d8

#endif // C3D8_HPP