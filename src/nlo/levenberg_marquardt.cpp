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

  std::vector<Eigen::Vector3d> nodes_k0 = pos_->GetTarget();

  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k1 = std::make_pair(qvec, tvec);
  std::pair<std::pair<Eigen::Vector4d, Eigen::Vector3d>, std::vector<Eigen::Vector3d>> sim = pos_->SimulateTransformToPose(pose_k1, scale);

  std::vector<Eigen::Vector3d> nodes_k_front, nodes_k_back;
  for (unsigned int i=0; i<sim.second.size(); i++) {
    if (i < sim.second.size()/2) {
      nodes_k_front.push_back(sim.second[i]);
    } else {
      nodes_k_back.push_back(sim.second[i]);
    }
  }

  out_f.resize(1);
  out_f[0] = fea_->ComputeStrainEnergy(nodes_k0, sim.second);
}


std::pair<double, Eigen::VectorXd> LevenbergMarquardt::Optimize(const Eigen::VectorXd params0) {

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
