#ifndef NLO_HPP
#define NLO_HPP

#include <eigen3/Eigen/Dense>
#include "../fea/fea.hpp"
#include "../fea/pos.hpp"


#include <cmath>
#include <iostream>



class LevenbergMarquardt 
{

public:

  LevenbergMarquardt(
    POS* pos, FEA* fea,
    int maxIterations = 1000, double lambda = 0.001, 
    double damping_factor = 2.0, double tolerance = 0.0001);

  std::pair<double, Eigen::VectorXd> Optimize(
    const std::pair<Eigen::Matrix<double, 4, 1>, Eigen::Matrix<double, 3, 1> > params_in);

  double ComputeResidual(const Eigen::VectorXd& pose);

  double GetResidual();

private:


  Eigen::VectorXd ComputeJacobian(const Eigen::VectorXd& params, double residual_original, double epsilon = 0.0001);
  Eigen::VectorXd ComputeJacobian3d(const Eigen::VectorXd& params, double residual_original,
                                    double delta_q = 0.0001, double delta_t = 0.01, double delta_s = 0.0001);

  double InitializeLambda(const Eigen::VectorXd& jacobian, double initial_residual);

  Eigen::VectorXd PosePairToMatrix(const std::pair<Eigen::Vector4d, Eigen::Vector3d> pose);
  std::pair<Eigen::Vector4d, Eigen::Vector3d> MatrixToPosePair(const Eigen::VectorXd pose_vec);
  double ScaleFromMatrix(const Eigen::VectorXd pose_vec);

  POS* pos_;
  FEA* fea_;

  double residual_;

  int max_iterations_;
  double lambda_;
  double tolerance_;
  double damping_factor_;

  bool verbose_ = false;




};


#endif