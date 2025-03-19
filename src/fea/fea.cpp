#include "fea.hpp"

FEA::FEA(std::string element_type,
         float young_modulus, float poisson_coefficient, 
         bool debug_mode) : E_(young_modulus), nu_(poisson_coefficient), debug_mode_(debug_mode) {

  if (element_type == "C2D4") {
    element_ = new C2D4(young_modulus, poisson_coefficient);
  } else if (element_type == "C3D6") {
    element_ = new C3D6(young_modulus, poisson_coefficient);
  } else if (element_type == "C3D8") {
    element_ = new C3D8(young_modulus, poisson_coefficient);
  } else {
    std::cout << "Element not supported" << std::endl;
  }

  fea_data_ = new FEAData();
  fea_data_->element_name = element_->getElementName();

  if (debug_mode_) {
    std::cout << std::endl;
    std::cout << "Finite Element Analysis: " << element_->getElementName() << std::endl;
    std::cout << " - Nodes per element: " << element_->getNumNodes() << std::endl;
    std::cout << " - DOF per node: " << element_->getDofPerNode() << std::endl;
    //std::cout << " - Material Matrix D: " << std::endl << element_->getD() << std::endl;
  }
}


void FEA::MatAssembly(std::vector<Eigen::Vector2d> &vpts, 
                      std::vector<std::vector<unsigned int>> &velts) {

  std::vector<Eigen::Vector3d> vpts3d;
  for (auto v : vpts) {
    vpts3d.push_back(Eigen::Vector3d(v[0], v[1], 0.0));
  }
  
  K_ = element_->matAssembly(vpts3d, velts);

  // Reset F_ and U_ to zero
  F_ = Eigen::MatrixXd::Zero(K_.rows(), 1);
  U_ = Eigen::MatrixXd::Zero(K_.rows(), 1);

  fea_data_->number_of_nodes = vpts.size();
  fea_data_->number_of_elements = velts.size();
}


void FEA::MatAssembly(std::vector<Eigen::Vector3d> &vpts, 
                      std::vector<std::vector<unsigned int>> &velts) {
  K_ = element_->matAssembly(vpts, velts);

  // Reset F_ and U_ to zero
  F_ = Eigen::MatrixXd::Zero(K_.rows(), 1);
  U_ = Eigen::MatrixXd::Zero(K_.rows(), 1);

  bc_ = std::vector<bool>(vpts.size(), false);
  loads_ = std::vector<bool>(vpts.size(), false);
}




void FEA::ApplyBoundaryConditions(BoundaryConditions &bc) {

  int bc_encastre = 0;
  int bc_displ = 0;

  int num_dof = bc.NumDof();

  for (unsigned int node = 0; node < bc.NodeIds().size(); node++) {

    if (!bc.NodeIds()[node])
      continue;
    
    std::vector<double> values = bc.Values(node);

    bool encastre = std::all_of(values.begin(), values.end(), [](double val) { return val == 0.0; });
    
    for (unsigned int i=0; i<values.size(); i++) {
      unsigned int m = num_dof*node + i;

      //if (!encastre && values[i] == 0) {
      //  continue;
      //}

      for (unsigned int j=0; j<K_.cols(); j++) {
        K_(m,j) = 0.0;
        if (!encastre) {
          F_(j,0) -= K_(j,m) * values[i];
        }
        K_(j,m) = 0.0;
      }

      if (encastre) {
        bc_encastre++;
      } else {
        bc_displ++;
      }

      K_(m,m) = 1;
      F_(m,0) = values[i];
    }
  }

  std::cout << " - Encastre conditions: " << bc_encastre << std::endl;
  std::cout << " - Displacement conditions: " << bc_displ << std::endl;
}





void FEA::ApplyLoads(Loads &loads) {

  int loads_node = 0;

  int num_dof = loads.NumDof();

  for (unsigned int node = 0; node < loads.NodeIds().size(); node++) {
    if (!loads.NodeIds()[node])
      continue;
    
    std::vector<double> values = loads.Values(node);

    for (unsigned int i=0; i<values.size(); i++) {
      if (values[i] == 0) {
        continue;
      }

      unsigned int m = num_dof*node + i;
      
      F_(m,0) += values[i];
      
      loads_node++;
    }
  }

  std::cout << " - Node loads: " << loads_node << std::endl;
}



void FEA::Solve(std::string method) {
  // Check if K is singular or ill-conditioned
  if (solver::IsSingularOrIllConditioned2(K_) == 1) {
    return;
  }

  std::cout << " - Solving system with " << method << " method" << std::endl;
  U_ = solver::SolveSystemWithPreconditioning(K_, F_, method);

  Fi_ = K_ * U_;
}


void FEA::PostProcess(std::vector<Eigen::Vector2d> &vpts, 
                      std::vector<std::vector<unsigned int>> &velts) {

  std::vector<Eigen::Vector3d> vpts3d;
  for (auto v : vpts) {
    vpts3d.push_back(Eigen::Vector3d(v[0], v[1], 0.0));
  }

  PostProcess(vpts3d, velts);
}

void FEA::PostProcess(std::vector<Eigen::Vector3d> &vpts, 
                      std::vector<std::vector<unsigned int>> &velts) {
  fea_data_->number_of_nodes = vpts.size();
  fea_data_->number_of_elements = velts.size();
  
  element_->postProcess(vpts, velts, U_, *fea_data_);
  fea_data_->strain_energy = ComputeStrainEnergy();
}



double FEA::ComputeStrainEnergy() {
  sE_ = 0.5 * (U_.transpose() * F_)(0,0);
  return sE_;
}


double FEA::ComputeStrainEnergy(std::vector<Eigen::Vector3d> &u0,
                                std::vector<Eigen::Vector3d> &u1) {
  int dim_in = u0.size() * 3;
  int dim_k = K_.rows();

  //std::cout << "K" << std::endl << K_ << std::endl;

  if (dim_in != dim_k) {
    std::cout << "Error1: dim_in != dim_k: " << dim_in << " != " << dim_k << std::endl;
    return -1.0;
  }

  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(dim_in, 1);
  for (unsigned int i = 0; i < u0.size()/2; i++) {
    U(i*3+0, 0) =  u1[i][0] - u0[i][0];
    U(i*3+1, 0) =  u1[i][1] - u0[i][1];
    U(i*3+2, 0) =  u1[i][2] - u0[i][2];
  }

  Eigen::MatrixXd Fi = K_ * U;

  //std::ofstream file;
  //file.open("../data/meshmeld01/F.txt");
  //file << Fi << std::endl;
  //file << std::endl << "F ( " << Fi.rows() << ", " << Fi.cols() << " ) " << std::endl;
  //file.close();
  //file.open("../data/meshmeld01/U.txt");
  //file << U << std::endl;
  //file << std::endl << "U ( " << U.rows() << ", " << U.cols() << " ) " << std::endl;
  //file.close();

  double sE = 0.5 * (U.transpose() * Fi)(0,0);

  return sE;
}



double FEA::ComputeStrainEnergyPoseOptimization(
  std::vector<Eigen::Vector3d> &u0,
  std::vector<Eigen::Vector3d> &u1) {

  int dim_in = u0.size() * 3;
  int dim_k = K_.rows();

  if (dim_in != dim_k) {
    std::cout << "Error2: dim_in != dim_k: " << dim_in << " != " << dim_k << std::endl;
    return -1.0;
  }

  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(dim_in, 1);
  for (unsigned int i = 0; i < u0.size()/2; i++) {
    U(i*3, 0) = u1[i][0] - u0[i][0];
    U(i*3+1, 0) = u1[i][1] - u0[i][1];
    U(i*3+2, 0) = u1[i][2] - u0[i][2];
  }

  Eigen::MatrixXd Fi = K_ * U;

  double sE = 0.5 * (U.transpose() * Fi)(0,0);

  return sE;
}






// 
// LEGACY FEA
// 


void FEA::EncastreBackLayer() {
  std::vector<int> dir;
  for (unsigned int i = K_.rows()/2; i < K_.rows(); i++) {
    // Set row and column to zero
    for (unsigned int j = 0; j < K_.cols(); j++) {
      K_(i,j) = 0.0;
      K_(j,i) = 0.0;
    }
    // Set diagonal to k_large_
    K_(i,i) = k_large_;
  }
}


void FEA::ImposeDirichletEncastre(std::vector<int> &dir) {
  for (auto d : dir) {
    int mp0 = 3*(d - 1);
    int mp1 = mp0 + 1;
    int mp2 = mp0 + 2;

    K_(mp0,mp0) = k_large_;
    K_(mp1,mp1) = k_large_;
    K_(mp2,mp2) = k_large_;

    F_(mp0,0) = 0.0;
    F_(mp1,0) = 0.0;
    F_(mp2,0) = 0.0;
  }
}


void FEA::ImposeDirichletEncastre(std::vector<std::vector<int>> &dir) {
  for (auto d : dir) {
    int mp0 = 3*(d[0] - 1);
    int mp1 = mp0 + 1;
    int mp2 = mp0 + 2;

    K_(mp0,mp0) = k_large_;
    K_(mp1,mp1) = k_large_;
    K_(mp2,mp2) = k_large_;

    F_(mp0,0) = 0.0;
    F_(mp1,0) = 0.0;
    F_(mp2,0) = 0.0;
  }

}


void FEA::ComputeForces() {
  F_ = Eigen::MatrixXd::Zero(K_.rows(), 1);
  F_ = K_ * U_;
}


void FEA::SetForces(std::vector<std::vector<float>> &vF) {
  F_ = Eigen::MatrixXd::Zero(3*vF.size(), 1);
  for (unsigned int i = 0; i < vF.size(); i++) {
    F_(i*3, 0) = vF[i][0];
    F_(i*3+1, 0) = vF[i][1];
    F_(i*3+2, 0) = vF[i][2];
  }
}


void FEA::ComputeDisplacements() {
  K1_ = K_.inverse();
  U_ = K1_ * F_;
}


// 
// LEGACY FEA
// 





// 
// REPORT FUNCTIONS
// 

std::string buildHeader() {
  int len_field = 12;
  std::vector<std::string> hdrs = {"n", "f", "f.x", "f.y", "f.z", "fi", "fi.x", "fi.y", "fi.z", "u", "u.x", "u.y", "u.z"};

  std::string header;
  for (auto h : hdrs) {
    for (unsigned int i = 0; i < len_field - h.size(); i++) {
      header += " ";
    }
    header += h + " ";
  }

  return header;
}


void FEA::ReportNodes(std::string filename) {

  // Eigen Matrix of size U_.rows()/3 x 7
  // n, u.x, u.y, u.z, f.x, f.y, f.z
  Eigen::MatrixXd exp = Eigen::MatrixXd::Zero(U_.rows()/3, 13);

  for (unsigned int n=0; n<exp.rows(); n++) {
    exp(n,0) = n;

    exp(n,2) = F_(n*3, 0);
    exp(n,3) = F_(n*3+1, 0);
    exp(n,4) = F_(n*3+2, 0);
    exp(n,1) = sqrt(exp(n,2)*exp(n,2) + exp(n,3)*exp(n,3) + exp(n,4)*exp(n,4));

    exp(n,6) = Fi_(n*3, 0);
    exp(n,7) = Fi_(n*3+1, 0);
    exp(n,8) = Fi_(n*3+2, 0);
    exp(n,5) = sqrt(exp(n,6)*exp(n,6) + exp(n,7)*exp(n,7) + exp(n,8)*exp(n,8));

    exp(n,10) = U_(n*3, 0);
    exp(n,11) = U_(n*3+1, 0);
    exp(n,12) = U_(n*3+2, 0);
    exp(n,9) = sqrt(exp(n,10)*exp(n,10) + exp(n,11)*exp(n,11) + exp(n,12)*exp(n,12));
  }

  std::ofstream file;
  file.open(filename);
  file << buildHeader() << std::endl;
  file << exp << std::endl;

  file.close();

  std::cout << " * ReportNodes [U,F]: " << filename << std::endl;
}

void FEA::ExportAll(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << K_ << std::endl;
  file << std::endl << "K ( " << K_.rows() << ", " << K_.cols() << " ) " << std::endl << std::endl;
  file << F_ << std::endl;
  file << std::endl << "F ( " << F_.rows() << ", " << F_.cols() << " ) " << std::endl << std::endl;
  file << U_ << std::endl;
  file << std::endl << "U ( " << U_.rows() << ", " << U_.cols() << " ) " << std::endl << std::endl;
  file.close();

  std::cout << " * ExportAll [K,F,U]: " << filename << std::endl;
}

void FEA::ExportSystem(std::string filename) {
  
  Eigen::MatrixXd KUF = Eigen::MatrixXd::Zero(K_.rows(), K_.cols() + 2);
  KUF.block(0, 0, K_.rows(), K_.cols()) = K_;
  KUF.block(0, K_.cols(), K_.rows(), 1) = U_;
  KUF.block(0, K_.cols()+1, K_.rows(), 1) = F_;
  std::cout << "KUF" << std::endl << KUF << std::endl;
  
  std::ofstream file;

  file.open(filename);

  for (unsigned int i = 0; i < K_.rows(); i++) {
    for (unsigned int j = 0; j < K_.cols(); j++) {
      file << K_(i,j) << " ";
    }
    file << U_(i,0) << " ";
    file << F_(i,0) << " ";
    file << std::endl;
  }


  file.close();

  std::cout << " * ExportSystem KU=F: " << filename << std::endl;
}

void FEA::ExportK(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << K_ << std::endl;
  file << std::endl << "K ( " << K_.rows() << ", " << K_.cols() << " ) " << std::endl;
  file.close();

  std::cout << " * ExportK [K]: " << filename << std::endl;
}

void FEA::ExportF(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << F_ << std::endl;
  file << std::endl << "F ( " << F_.rows() << ", " << F_.cols() << " ) " << std::endl;
  file.close();

  std::cout << " * ExportF [F]: " << filename << std::endl;
}

void FEA::ExportU(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << U_ << std::endl;
  file << std::endl << "U ( " << U_.rows() << ", " << U_.cols() << " ) " << std::endl;
  file.close();

  std::cout << " * ExportU [U]: " << filename << std::endl;
}


void FEA::PrintK() {
  std::cout << std::endl
            << "K_ = " << std::endl
            << K_ << std::endl;
}


void FEA::PrintEigenvaluesK() {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K_);
  std::cout << std::endl
            << "The eigenvalues of K_ are:" << std::endl
            << es.eigenvalues().transpose() << std::endl;

}


void FEA::ReportFEAData(std::string filename) {
  std::ofstream file;
  file.open(filename);

  file << "FEA Data" << std::endl;
  file << "---------" << std::endl;
  file << "Element: " << fea_data_->element_name << std::endl;
  file << "Number of nodes: " << fea_data_->number_of_nodes << std::endl;
  file << "Number of elements: " << fea_data_->number_of_elements << std::endl;
  file << "Strain energy: " << fea_data_->strain_energy << std::endl;
  file << std::endl;

  file << "\nStrain min: ";
  for (auto e : fea_data_->emin) {
    file << e << " ";
  }
  file << "\nStrain max: ";
  for (auto e : fea_data_->emax) {
    file << e << " ";
  }
  file << "\nStress min: ";
  for (auto s : fea_data_->smin) {
    file << s << " ";
  }
  file << "\nStress max: ";
  for (auto s : fea_data_->smax) {
    file << s << " ";
  }
  file << "\nU min: ";
  for (auto u : fea_data_->umin) {
    file << u << " ";
  }
  file << "\nU max: ";
  for (auto u : fea_data_->umax) {
    file << u << " ";
  }
  file << std::endl << std::endl;

  

  file << "\nNode Dislacements: " << std::endl;
  file << "-----------------" << std::endl;
  file << "n, u, u.x, u.y, u.z" << std::endl;
  for (unsigned int n=0; n<fea_data_->number_of_nodes; n++) {
    file << n << ", ";
    file << sqrt(U_(n*3, 0)*U_(n*3, 0) + U_(n*3+1, 0)*U_(n*3+1, 0) + U_(n*3+2, 0)*U_(n*3+2, 0)) << ", ";
    for (unsigned int i=0; i<3; i++) {
      file << U_(n*3+i, 0) << ", ";
    }
    file << std::endl;
  }


  file << "\nNode Strain: " << std::endl;
  file << "-----------------" << std::endl;
  file << "n, e, e.x, e.y, e.z" << std::endl;
  for (unsigned int n=0; n<fea_data_->node_strain.size(); n++) {
    file << n << ", ";
    file << sqrt(fea_data_->node_strain[n][0]*fea_data_->node_strain[n][0] + 
                 fea_data_->node_strain[n][1]*fea_data_->node_strain[n][1] + 
                 fea_data_->node_strain[n][2]*fea_data_->node_strain[n][2]) << ", ";
    for (auto n : fea_data_->node_strain[n]) {
      file << n << ", ";
    }
    file << std::endl;
  }

  file << "\nNode Stress: " << std::endl;
  file << "-----------------" << std::endl;
  file << "n, s, s.x, s.y, s.z" << std::endl;
  for (unsigned int n=0; n<fea_data_->node_stress.size(); n++) {
    file << n << ", ";
    file << sqrt(fea_data_->node_stress[n][0]*fea_data_->node_stress[n][0] + 
                 fea_data_->node_stress[n][1]*fea_data_->node_stress[n][1] + 
                 fea_data_->node_stress[n][2]*fea_data_->node_stress[n][2]) << ", ";
    for (auto n : fea_data_->node_stress[n]) {
      file << n << ", ";
    }
    file << std::endl;
  }

  file.close();

  std::cout << std::endl;
  std::cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
  std::cout << " - ReportFEAData: " << filename << std::endl;
  std::cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
  std::cout << " - Element: " << fea_data_->element_name << std::endl;
  std::cout << " - Number of nodes: " << fea_data_->number_of_nodes << std::endl;
  std::cout << " - Number of elements: " << fea_data_->number_of_elements << std::endl;
  std::cout << " - Strain energy: " << fea_data_->strain_energy << std::endl;
  std::cout << " - Strain min: ";
  for (auto e : fea_data_->emin) {
    std::cout << e << " ";
  }
  std::cout << std::endl;
  std::cout << " - Strain max: ";
  for (auto e : fea_data_->emax) {
    std::cout << e << " ";
  }
  std::cout << std::endl;
  std::cout << " - Stress min: ";
  for (auto s : fea_data_->smin) {
    std::cout << s << " ";
  }
  std::cout << std::endl;
  std::cout << " - Stress max: ";
  for (auto s : fea_data_->smax) {
    std::cout << s << " ";
  }
  std::cout << std::endl;
  std::cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
  std::cout << std::endl;

}