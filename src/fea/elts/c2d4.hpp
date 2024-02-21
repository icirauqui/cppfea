#ifndef C2D4_HPP
#define C2D4_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>

#include "element2d.hpp"

class C2D4 : public Element2D {

public:

    C2D4(double E, double nu);

    // Function to compute shape functions for a triangular prism
    Eigen::VectorXd computeShapeFunctions(double xi, double eta, double zeta);

    // Function to compute the derivative of shape functions for a triangular prism
    Eigen::MatrixXd computeShapeFunctionDerivatives(double xi, double eta, double zeta);

    // Function to compute the Strain-Displacement Matrix (B)
    Eigen::MatrixXd computeStrainDisplacementMatrix(Eigen::MatrixXd &dN, Eigen::MatrixXd &invJ, double &detJ);

}; // class C2D4

#endif // C2D4_HPP