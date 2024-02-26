#include "fea/fea.hpp"
#include "fea/fem.hpp"
#include "fea/pos.hpp"
#include "vis/vis.hpp"

#include "dataset/dataset.hpp"
#include "dataset/simulation_models.cpp"



void simulation_manual() {
  std::cout << "\nSimulation Manual" << std::endl;
  std::string element = "C3D6";



  std::cout << "\nLoad Data" << std::endl;
  SimulationC3D8_1 ds(1);
  std::vector<Eigen::Vector3d> vpts = ds.GetNodes();


  
  std::cout << "\nCreate FEM objects and add points" << std::endl;
  FEM fem1(element);
  FEM fem2(element);
  for (unsigned int i=0; i<vpts.size(); i++) {
    fem1.AddPoint(Eigen::Vector3d(vpts[i][0], vpts[i][1], vpts[i][2]));
    fem2.AddPoint(Eigen::Vector3d(vpts[i][0], vpts[i][1], vpts[i][2]));
  }



  std::cout << "\nTriangulate and compute poses" << std::endl;
  fem1.Compute(true);
  fem1.ComputeExtrusion();
  fem1.ViewMesh(true, 0);
  fem2.InitCloud();
  fem2.SetExtrusion(fem1.GetExtrusionDelta(), fem1.GetElementHeight());
  //fem1.ViewMesh(true, fem2.GetCloud(), fem2.GetExtrusion(), fem2.GetPose(), 0);



  std::cout << "\nTransform pose 2 for simulation, impose a rotation of x degrees around each axis" << std::endl;
  POS pos(fem2.GetEigenNodes(), fem2.GetPose());
  pos.SetTarget(fem1.GetEigenNodes());

  double ang = 5*M_PI/180;
  Eigen::Vector3d axis(1,1,1);

  Eigen::Vector4d imposed_angle_q = pos.QuaternionFromAngleAxis(axis, ang);
  Eigen::Vector3d model_offset(1.0, 1.0, 1.0); 
  pos.Transform(imposed_angle_q, model_offset, 1.0);
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k0 = pos.GetPose();
  //std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes_k0 = pos.GetPointLayers();
  //fem1.ViewMesh(true, nodes_k0.first, nodes_k0.second, pose_k0, 0);



  std::cout << "\nInitialize FEA object, use conectivity to compute stiffness matrix" << std::endl;
  FEA fea(element, 10000.0, 0.3, true);
  std::vector<Eigen::Vector3d> nodes = fem1.GetEigenNodes();
  std::vector<std::vector<unsigned int>> elements = fem1.GetElements();
  fea.MatAssembly(nodes, elements);



  std::cout << "\nSimulate displacement steps, computing strain energy in each" << std::endl;
  int num_steps = 5;
  Eigen::Vector3d step_t = -model_offset / num_steps;
  Eigen::Vector4d step_r = pos.QuaternionFromAngleAxis(axis, -ang/num_steps);
  std::vector<Eigen::Vector3d> nodes_target = pos.GetTarget();

  for (unsigned int step = 0; step < num_steps + 1; step++) {
    std::cout << " - Step " << step << " / " << num_steps << std::endl;

    if (step > 0) {
      pos.Transform(step_r, step_t, 1.0);
    } 

    std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_i = pos.GetPose();
    std::vector<Eigen::Vector3d> nodes_i = pos.GetPoints();

    double sE = fea.ComputeStrainEnergy(nodes_target, nodes_i);
    std::cout << "   Strain energy = " << sE << std::endl;

    std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes_layers = pos.GetPointLayers();
    fem1.ViewMesh(true, nodes_layers.first, nodes_layers.second, pose_i, 0);
  }
}




/*
void simulation_optimizer() {
  std::cout << "Simulation with Optimizer" << std::endl;

  // 1. Load data
  Dataset ds("../data", element);
  std::vector<std::vector<float>> vpts = ds.points();



  // 2. Create FEM objects and add points
  FEM fem1(element);
  FEM fem2(element);

  for (int i = 0; i < vpts.size(); i++) {
    Eigen::Vector3d pt1(vpts[i][0], vpts[i][1], vpts[i][2]);
    Eigen::Vector3d pt2(vpts[i][0], vpts[i][1], vpts[i][2]);
    //for (unsigned int j=0; j<3; j++) {
    //  //pt1(j) += noise_multiplier * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    //  pt2(j) += noise_multiplier * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    //}
    fem1.AddPoint(pt1);
    fem2.AddPoint(pt2);
  }



  // 3. Triangulate andcompute poses.
  fem1.Compute(true);
  fem2.InitCloud();
  fem1.ComputeExtrusion();
  fem2.SetExtrusion(fem1.GetExtrusionDelta(), fem1.GetElementHeight());

  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose1 = ApproximatePose(fem1.GetEigenNodes());
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2 = ApproximatePose(fem2.GetEigenNodes());

  fem1.ViewMesh(true, fem2.GetCloud(), fem2.GetExtrusion(), pose1, pose2, view_quick);

  POS pos_tmp(fem1.GetNodes(), pose1);
  std::vector<Eigen::Vector3d> nodes_k0 = pos_tmp.GetPoints();




  // 4. Transform Pose 2 for simulation, impose a rotation of x degrees around each axis
  POS pos(fem2.GetNodes(), pose2);
  pos.SetTarget(fem1.GetNodes());

  double ang = 5*M_PI/180;
  Eigen::Vector3d axis(1,1,1);
  Eigen::Vector4d imposed_angle_q = pos.QuaternionFromAngleAxis(axis, ang);
  std::pair<Eigen::Vector4d, Eigen::Vector3d> latest_pose = pos.GetPose();
  pos.Transform(imposed_angle_q, model_offset, 1.0);

  std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes_k0_layers = pos.GetPointLayers();
  std::pair<Eigen::Vector4d, Eigen::Vector3d> new_pose = pos.GetPose();
  fem1.ViewMesh(true, nodes_k0_layers.first, nodes_k0_layers.second, pose1, new_pose, view_quick);



  //pos.TransformToPose(pose1, 1.0);
  //std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes_k1_layers = pos.GetPointLayers();
  //std::pair<Eigen::Vector4d, Eigen::Vector3d> new_pose1 = pos.GetPose();
  //fem1.ViewMesh(true, nodes_k1_layers.first, nodes_k1_layers.second, pose1, new_pose1, view_slow);


  Eigen::Vector4d applied_rotation = pos.ComputeQuaternionRotation(latest_pose.first, new_pose.first);


  // 5. Initialize FEA object, use conectivity to compute stiffness matrix
  FEA fea(0, element, E, nu, depth, fg, false);
  std::vector<std::vector<float>> nodes = fem1.GetNodes();
  std::vector<std::vector<int>> elements = fem1.GetElements();
  fea.MatAssembly(nodes, elements);


  // Initialize optimizer
  LevenbergMarquardt lm(&pos, &fea, 5, 10.0, 2.0, 0.01);
  std::vector<Eigen::Vector3d> nodes_k0 = pos.GetTarget();
  std::vector<Eigen::Vector3d> nodes_k1 = pos.GetPoints();
  double sE0 = fea.ComputeStrainEnergy(nodes_k0, nodes_k1);

  Eigen::VectorXd pose_k0 = pos.GetPoseVector();
  std::pair<double, Eigen::VectorXd> res = lm.Optimize(pose_k0);

  std::cout << "Len History 3 = " << pos.LenHistory() << std::endl;

  //double sE1 = lm.GetResidual();
  double sE1 = res.first;
  std::cout << "  Strain energy = " << sE0 << " -> " << sE1 << std::endl;
  std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes_k2_layers = pos.GetPointLayers();
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k2 = pos.GetPose();
  fem1.ViewMesh(true, nodes_k2_layers.first, nodes_k2_layers.second, pose1, pose_k2, 1);



  Eigen::VectorXd pose_k2v = pos.GetPoseVector();

  std::cout << "Pose in = " << pose1.first.transpose() << " " << pose1.second.transpose() << std::endl;
  std::cout << "Pose k0 = " << pose_k0.transpose() << std::endl;
  std::cout << "Pose k2 = " << pose_k2v.transpose() << std::endl;

  std::vector<Eigen::Vector3d> new_nodes_k2 = pos.GetPoints();
  std::vector<Eigen::Vector3d> new_nodes_k2_front, new_nodes_k2_back;
  for (unsigned int i=0; i<new_nodes_k2.size(); i++) {
    if (i < new_nodes_k2.size()/2) {
      new_nodes_k2_front.push_back(new_nodes_k2[i]);
    } else {
      new_nodes_k2_back.push_back(new_nodes_k2[i]);
    }
  }
  std::pair<Eigen::Vector4d, Eigen::Vector3d> new_pose_k2 = pos.GetPose();

  double sE_man = fea.ComputeStrainEnergy(nodes_k0, new_nodes_k2);
  std::cout << "Strain Energy Manual: " << sE0 << " -> " << sE_man << std::endl;

  //std::cout << "A" << std::endl;
  //for (unsigned int i=0; i<20; i++) {
  //  std::cout << new_nodes_front[i].transpose()
  //            << "\t | \t" 
  //            << new_nodes_k2_front[i].transpose()
  //            << "\t"
  //            << std::endl;
  //}

  fem1.ViewMesh(true, new_nodes_k2_front, new_nodes_k2_back, pose1, new_pose_k2, 10);
}
*/
