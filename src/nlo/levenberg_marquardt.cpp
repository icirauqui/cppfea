#include "levenberg_marquardt.hpp"

LevenbergMarquardt::LevenbergMarquardt(
  POS* pos, FEA* fea,
  int maxIterations, double lambda, double damping_factor, double tolerance)
  : pos_(pos), fea_(fea),
    max_iterations_(maxIterations), lambda_(lambda), damping_factor_(damping_factor), tolerance_(tolerance) {}



void LevenbergMarquardt::ComputeResidual(const mrpt::math::CVectorDouble& x, 
                                         const mrpt::math::CVectorDouble& y, 
                                         mrpt::math::CVectorDouble& out_f) {
  Eigen::Vector3d tvec(x[0], x[1], x[2]);
  Eigen::Vector4d qvec(x[3], x[4], x[5], x[6]);
  double scale = x[7];
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k1 = std::make_pair(qvec, tvec);
  std::pair<std::pair<Eigen::Vector4d, Eigen::Vector3d>, std::vector<Eigen::Vector3d>> sim = pos_->SimulateTransformToPose(pose_k1, scale);

  std::vector<Eigen::Vector3d> nodes_k0 = pos_->GetTarget();

  out_f.resize(1);
  out_f[0] = fea_->ComputeStrainEnergy(nodes_k0, sim.second);

  //std::cout << std::endl;
  //std::cout << " - Levenberg-Marquardt: Strain energy = " << out_f[0] << std::endl;
  //std::cout << " - Levenberg-Marquardt: Pose = " << pose_k1.first.transpose() << " " << pose_k1.second.transpose() << std::endl;
  //std::cout << " - Levenberg-Marquardt: Scale = " << scale << std::endl;
  //std::cout << " - Levenberg-Marquardt: Residual = " << out_f[0] << std::endl;

}



double LevenbergMarquardt::ComputeResidual(const Eigen::VectorXd& pose) {
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_pair = MatrixToPosePair(pose);
  double scale = ScaleFromMatrix(pose);

  std::pair<std::pair<Eigen::Vector4d, Eigen::Vector3d>, std::vector<Eigen::Vector3d>> sim = pos_->SimulateTransformToPose(pose_pair, scale);

  std::vector<Eigen::Vector3d> nodes_k0 = pos_->GetTarget();

  double residual = fea_->ComputeStrainEnergy(nodes_k0, sim.second);

  //std::cout << std::endl;
  //std::cout << " - Levenberg-Marquardt: Strain energy = " << residual << std::endl;
  //std::cout << " - Levenberg-Marquardt: Pose = " << pose_pair.first.transpose() << " | " << pose_pair.second.transpose() << " | " << scale << std::endl;

  return residual;
}





std::pair<double, Eigen::VectorXd> LevenbergMarquardt::Optimize(
  const std::pair<Eigen::Matrix<double, 4, 1>, Eigen::Matrix<double, 3, 1> > params_in) {

  Eigen::VectorXd params0 = PosePairToMatrix(params_in);
  residual_ = ComputeResidual(params0);
  double jacobian_grad = 1 / residual_;

  // Initialize lambda
  Eigen::VectorXd jacobian_init = ComputeJacobian(params0, residual_, jacobian_grad);
  //std::cout << "   Jacobian Init = " << jacobian_init.transpose() << std::endl;
  //std::cout << "   Jacobian Grad = " << jacobian_grad << std::endl;
  if (lambda_ < 0.0) {
    lambda_ = InitializeLambda(jacobian_init, residual_);
  }

  // set lambda to max value of the Jacobian
  //lambda_ = jacobian_init.maxCoeff();

  std::cout << "   Levenberg-Marquardt: Initialization | Residual = " << residual_ << std::endl;
  std::cout << "      Lambda = " << lambda_ << std::endl;

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(params0.size(), params0.size());
  bool converged = false;
  int iterations = 0;

  while (!converged && iterations < max_iterations_) {
    std::cout << "   Levenberg-Marquardt: Iteration " << iterations+1 << std::endl;

    Eigen::VectorXd jacobian = ComputeJacobian(params0, residual_, 0.0000001);
    std::cout << "      Jacobian = " << jacobian.transpose() << std::endl;
    Eigen::MatrixXd H = jacobian * jacobian.transpose();
    //std::cout << "      Hessian = " << H << std::endl;
    Eigen::VectorXd g = jacobian * residual_;
    std::cout << "      Gradient = " << g.transpose() << std::endl;

    Eigen::VectorXd delta = (H + lambda_ * I).ldlt().solve(-g);
    std::cout << "      Delta = " << delta.transpose() << std::endl;
    Eigen::VectorXd params1 = params0 + delta;
    double residual1 = ComputeResidual(params1);
    std::cout << "      residual1 = " << residual1 << std::endl;

    if (residual1 < residual_) {
      lambda_ /= damping_factor_;
      params0 = params1;
      residual_ = residual1;
    } else {
      lambda_ *= damping_factor_;
    }

    if (delta.norm() < tolerance_) {
      converged = true;
    }

    std::cout << "      Residual = " << residual_ << std::endl;

    iterations++;
  }

  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k1 = MatrixToPosePair(params0);
  double scale = ScaleFromMatrix(params0);

  pos_->TransformToPose(pose_k1, scale);

  return std::make_pair(residual_, params0);
}


std::pair<double, Eigen::VectorXd> LevenbergMarquardt::OptimizeCLM(
  const std::pair<Eigen::Matrix<double, 4, 1>, Eigen::Matrix<double, 3, 1> > params_in) {

  Eigen::VectorXd params0(8);
  for (int i = 0; i < 3; i++) {
    params0(i) = params_in.second(i);
  }
  for (int i = 0; i < 4; i++) {
    params0(i+3) = params_in.first(i);
  }
  params0(7) = 1.0;

  return OptimizeCLM(params0);
}


std::pair<double, Eigen::VectorXd> LevenbergMarquardt::OptimizeCLM(const Eigen::VectorXd params0) {

  mrpt::math::CVectorDouble optimal_x;
  mrpt::math::CVectorDouble initial_x;
  mrpt::math::CVectorDouble y;

  initial_x.resize(8);
  for (int i = 0; i < 8; i++) {
    initial_x[i] = params0(i);
  }

  mrpt::math::CVectorDouble increments_x(8);
  increments_x.fill(0.0001);

  mrpt::math::CLevenbergMarquardt::TResultInfo info;
  mrpt::math::CLevenbergMarquardt lm;

  auto callback = [&](const mrpt::math::CVectorDouble& x, 
                      const mrpt::math::CVectorDouble& initial_x, 
                      mrpt::math::CVectorDouble& result) {
      ComputeResidual(x, initial_x, result);
  };

  lm.execute(optimal_x, initial_x, callback, increments_x, y, info);

  //std::cout << "Levenberg-Marquardt: final error: " << info.final_sqr_err << std::endl;
  double levmarq_final_error = info.final_sqr_err;

  Eigen::Vector3d tvec(optimal_x[0], optimal_x[1], optimal_x[2]);
  Eigen::Vector4d qvec(optimal_x[3], optimal_x[4], optimal_x[5], optimal_x[6]);
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k1 = std::make_pair(qvec, tvec);
  
  Eigen::VectorXd params = params0;
  for (int i = 0; i < 8; i++) {
    params(i) = optimal_x[i];
  }
  double scale = params(7);

  pos_->TransformToPose(pose_k1, scale);

  return std::make_pair(levmarq_final_error, params);
}




double LevenbergMarquardt::GetResidual() {
  return residual_;
}




Eigen::VectorXd LevenbergMarquardt::ComputeJacobian(const Eigen::VectorXd& params, double residual_original, double epsilon) {

  Eigen::VectorXd params_epsilon = params;
  Eigen::VectorXd jacobian(params.size());

  for (int i = 0; i < params.size(); i++) {
    params_epsilon = params;
    params_epsilon(i) += epsilon;
    double residual_e = ComputeResidual(params_epsilon);
    jacobian(i) = (residual_e - residual_original) / epsilon;
  }

  return jacobian;
}



Eigen::VectorXd LevenbergMarquardt::ComputeJacobian3d(const Eigen::VectorXd& params, double residual_original,
  double delta_q, double delta_t, double delta_s) {

  Eigen::VectorXd params_epsilon = params;

  // Vector of zeros
  Eigen::VectorXd params_delta = Eigen::VectorXd::Zero(params.size());
  params_delta(0) += delta_q;
  params_delta(1) += delta_q;
  params_delta(2) += delta_q;
  params_delta(3) += delta_q;
  params_delta(4) += delta_t;
  params_delta(5) += delta_t;
  params_delta(6) += delta_t;
  params_delta(7) += delta_s;

  Eigen::VectorXd jacobian(params.size());

  for (int i = 0; i < params.size(); i++) {
    params_epsilon = params;
    params_epsilon(i) += params_delta(i);
    double residual_e = ComputeResidual(params_epsilon);
    jacobian(i) = (residual_e - residual_original) / params_delta(i);
  }

  return jacobian;
}


double LevenbergMarquardt::InitializeLambda(const Eigen::VectorXd& jacobian, double initial_residual) {
    // Ensure the vector is not empty
    if (jacobian.size() == 0) {
        throw std::invalid_argument("Jacobian vector is empty.");
    }

    // Calculate the square of each element (assuming it represents the diagonal of the Jacobian)
    Eigen::VectorXd squaredJacobian = jacobian.array().square();
    
    // Initialize lambda as a small value times the maximum element in the squared Jacobian
    double lambda = (1/initial_residual) * squaredJacobian.maxCoeff();

    return lambda;
}


//double LevenbergMarquardt::InitializeLambda(const Eigen::VectorXd& jacobianVector, double initial_residual) {
//    // Assuming jacobianVector is the gradient of the residual w.r.t. parameters
//    Eigen::MatrixXd jacobianMatrix = Eigen::MatrixXd::Zero(8, 8);
//    
//    // Convert the vector to a diagonal matrix with squared elements
//    for(int i = 0; i < jacobianVector.size(); ++i) {
//        jacobianMatrix(i, i) = std::pow(jacobianVector(i), 2);
//    }
//    
//    // Initialize lambda based on the average of the diagonal elements
//    double avgDiagonal = jacobianMatrix.diagonal().mean();
//    double lambda = (1/initial_residual) * avgDiagonal; // Example scaling factor
//    
//    return lambda;
//}



Eigen::VectorXd LevenbergMarquardt::PosePairToMatrix(const std::pair<Eigen::Vector4d, Eigen::Vector3d> pose) {
  Eigen::VectorXd pose_vec(8);
  for (int i = 0; i < 4; i++) {
    pose_vec(i) = pose.first(i);
  }
  for (int i = 0; i < 3; i++) {
    pose_vec(i+4) = pose.second(i);
  }
  pose_vec(7) = 1.0;
  return pose_vec;
}

std::pair<Eigen::Vector4d, Eigen::Vector3d> LevenbergMarquardt::MatrixToPosePair(const Eigen::VectorXd pose_vec) {
  Eigen::Vector4d qvec(pose_vec(0), pose_vec(1), pose_vec(2), pose_vec(3));
  Eigen::Vector3d tvec(pose_vec(4), pose_vec(5), pose_vec(6));
  return std::make_pair(qvec, tvec);
}

double LevenbergMarquardt::ScaleFromMatrix(const Eigen::VectorXd pose_vec) {
  return pose_vec(7);
}