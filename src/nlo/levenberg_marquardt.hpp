#ifndef NLO_HPP
#define NLO_HPP

#include <eigen3/Eigen/Dense>
#include "../fea/fea.hpp"
#include "../fea/pos.hpp"

//#include "node.hpp"
//#include "edge.hpp"

#include <mrpt/math/CLevenbergMarquardt.h>
#include <mrpt/system/CTicTac.h>

#include <cmath>
#include <iostream>

//using namespace mrpt;
//using namespace mrpt::math;
//using namespace mrpt::system;


class LevenbergMarquardt 
{

public:

  LevenbergMarquardt(
    POS* pos, FEA* fea,
    int maxIterations = 1000, double lambda = 0.001, 
    double damping_factor = 2.0, double tolerance = 0.000001);

  std::pair<double, Eigen::VectorXd> Optimize(const Eigen::VectorXd params0);

  double GetResidual();

  void ComputeResidual(const mrpt::math::CVectorDouble& x, 
                       const mrpt::math::CVectorDouble& y, 
                       mrpt::math::CVectorDouble& out_f);

private:


  Eigen::MatrixXd ComputeJacobian(const Eigen::VectorXd& params, double residual_original, 
                                  double delta_t = 0.001, double delta_q = 0.0001, double delta_s = 0.00001);

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