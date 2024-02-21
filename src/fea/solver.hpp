#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/LU>
#include <limits>
#include <cmath>

namespace solver {

bool IsSingularOrIllConditioned1(Eigen::MatrixXd &A);
int IsSingularOrIllConditioned2(Eigen::MatrixXd &A);

Eigen::VectorXd SolveSystemWithPreconditioning(const Eigen::MatrixXd &A, const Eigen::MatrixXd &b, std::string method);

Eigen::VectorXd SolveSystemWithPreconditioningCG(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);
Eigen::VectorXd SolveSystemWithPreconditioningBiCGSTAB(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);

Eigen::VectorXd SolveSystemWithLU(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);
Eigen::VectorXd SolveSystemWithFullPivLU(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);

Eigen::VectorXd MatrixToVector(const Eigen::MatrixXd &columnMatrix);

}

#endif