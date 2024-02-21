#include "solver.hpp"


bool solver::IsSingularOrIllConditioned1(Eigen::MatrixXd &A) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
  double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
  if (cond > 1e10) {
    std::cout << "Condition number: " << cond << std::endl;
    return true;
  }
  return false;
}

int solver::IsSingularOrIllConditioned2(Eigen::MatrixXd &A) {
    // Threshold for determining if the determinant is effectively zero
    const double detThreshold = 1e-10;
    // Threshold for numerical instability in solutions, indicative of ill-conditioning
    const double solutionThreshold = 1e-10;

    // Check if the matrix is singular by determinant
    double det = A.determinant();
    if (std::fabs(det) < detThreshold) {
      std::cout << "   ERROR: Matrix is singular or nearly singular, determinant: " << det << std::endl;
      return 1; // Matrix is singular or near-singular
    }

    // Attempt to solve A*x = b for a unit vector b to check for ill-conditioning
    Eigen::VectorXd b = Eigen::VectorXd::Ones(A.rows());
    Eigen::VectorXd x = A.fullPivLu().solve(b);

    // Check if the solution is significantly different from what is expected
    Eigen::VectorXd b_check = A * x; // Recompute b to check solution
    if ((b - b_check).norm() / b.norm() > solutionThreshold) {
      std::cout << "   WARNING: Significant solution instability detected, indicating ill-conditioning." << std::endl;
      std::cout << "            Relative error: " << (b - b_check).norm() / b.norm() << std::endl;
      return 2; // Solution instability suggests ill-conditioning
    }

    // If neither singularity nor significant solution instability was detected
    return 0;
}



Eigen::VectorXd solver::SolveSystemWithPreconditioning(const Eigen::MatrixXd &A, const Eigen::MatrixXd &b, std::string method) {

  Eigen::VectorXd b_vec = MatrixToVector(b);

  if (method == "CG") {
    return SolveSystemWithPreconditioningCG(A, b_vec);
  } else if (method == "BiCGSTAB") {
    return SolveSystemWithPreconditioningBiCGSTAB(A, b_vec);
  } else if (method == "LU") {
    return SolveSystemWithLU(A, b_vec);
  } else if (method == "LUFull") {
    return SolveSystemWithFullPivLU(A, b_vec);
  } else {
    std::cout << "Method not supported" << std::endl;
    return Eigen::VectorXd::Zero(b.size());
  }

}

// Function to solve a system of linear equations A*x = b using an iterative solver with diagonal preconditioning
Eigen::VectorXd solver::SolveSystemWithPreconditioningCG(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
    // Convert the dense matrix A to a sparse matrix for efficiency with iterative solvers
    Eigen::SparseMatrix<double> Asparse = A.sparseView();

    // Initialize the iterative solver with diagonal preconditioning
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
    cg.compute(Asparse);

    // Check if the factorization was successful
    if(cg.info() != Eigen::Success) {
        std::cerr << "Matrix decomposition failed. Returning zero vector." << std::endl;
        return Eigen::VectorXd::Zero(b.size());
    }

    // Solve the system A*x = b
    Eigen::VectorXd x = cg.solve(b);

    // Return the solution vector
    return x;
}

// Function to solve a system of linear equations A*x = b using Eigen::BiCGSTAB with diagonal preconditioning
Eigen::VectorXd solver::SolveSystemWithPreconditioningBiCGSTAB(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
    // Convert the dense matrix A to a sparse matrix format for efficiency
    Eigen::SparseMatrix<double> Asparse = A.sparseView();

    // Initialize the BiCGSTAB solver with diagonal preconditioning for the sparse matrix
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> bicgstab;
    bicgstab.compute(Asparse);

    // Check if the decomposition was successful
    if (bicgstab.info() != Eigen::Success) {
        std::cerr << "Matrix decomposition failed. Returning zero vector." << std::endl;
        return Eigen::VectorXd::Zero(b.size());
    }

    // Solve the system A*x = b
    Eigen::VectorXd x = bicgstab.solve(b);

    return x;
}

Eigen::VectorXd solver::SolveSystemWithLU(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
    // Compute the LU decomposition of A
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);

    // Use the LU decomposition to solve for x
    Eigen::VectorXd x = lu.solve(b);

    return x;
}


Eigen::VectorXd solver::SolveSystemWithFullPivLU(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
    // Compute the LU decomposition of A with full pivoting
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);

    // Before solving, check if the solution is possible and the system is not singular
    if (!lu.isInvertible()) {
        std::cerr << "Matrix is singular or nearly singular, solution may not be accurate or possible." << std::endl;
        return Eigen::VectorXd::Zero(b.size()); // Return a zero vector if not invertible
    }

    // Use the LU decomposition to solve for x
    Eigen::VectorXd x = lu.solve(b);

    return x;
}

// Function to convert a column matrix (Eigen::MatrixXd) to Eigen::VectorXd
Eigen::VectorXd solver::MatrixToVector(const Eigen::MatrixXd &columnMatrix) {
    // Ensure the matrix is indeed a column matrix
    if (columnMatrix.cols() != 1) {
        std::cerr << "Error: The input matrix is not a column matrix." << std::endl;
        return Eigen::VectorXd(); // Return an empty vector in case of error
    }

    // Directly convert the column matrix to a vector
    Eigen::VectorXd vector = columnMatrix.col(0);

    return vector;
}
