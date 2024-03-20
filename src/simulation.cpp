#include "fea/fea.hpp"
#include "fea/fem.hpp"
#include "fea/pos.hpp"
//#include "vis/vis.hpp"

#include "dataset/dataset.hpp"
#include "dataset/simulation_models.cpp"
#include "nlo/levenberg_marquardt.hpp"



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



  std::cout << "\nSet Boundary Conditions" << std::endl;
  std::vector<unsigned int> bc_nodes = fem1.GetExtrusionIndices();
  BoundaryConditions3d bc(fea.NumDof(), &nodes);
  bc.AddEncastreByNodeIds(bc_nodes);
  fea.ApplyBoundaryConditions(bc);






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
    //double sE1 = fea.ComputeStrainEnergyPoseOptimization(nodes_target, nodes_i);
    std::cout << "   Strain energy = " << sE << std::endl;

    std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes_layers = pos.GetPointLayers();
    fem1.ViewMesh(true, nodes_layers.first, nodes_layers.second, pose_i, 0);
  }
}





void simulation_optimizer() {
  std::cout << "Simulation with Optimizer" << std::endl;
  std::string element = "C3D6";



  std::cout << "\nLoad Data" << std::endl;
  SimulationC3D8_1 ds(1);
  std::vector<Eigen::Vector3d> vpts = ds.GetNodes();


  
  std::cout << "\nCreate FEM objects and add points" << std::endl;
  FEM fem1(element);
  for (unsigned int i=0; i<vpts.size(); i++) {
    fem1.AddPoint(Eigen::Vector3d(vpts[i][0], vpts[i][1], vpts[i][2]));
  }  
  
  std::cout << "\nTriangulate and compute poses" << std::endl;
  fem1.Compute(true, true);
  fem1.ComputeExtrusion();
  fem1.ViewMesh(true, 0);


  FEM fem2(fem1);

  //FEM fem2(element);
  //std::vector<Eigen::Vector3d> fem1_alive_points = fem1.GetPoints(true);
  //for (unsigned int i=0; i<fem1_alive_points.size(); i++) {
  //  fem2.AddPoint(Eigen::Vector3d(fem1_alive_points[i][0], fem1_alive_points[i][1], fem1_alive_points[i][2]));
  //}
  //fem2.SetTriangles(fem1.GetTriangles());
  ////bool ok = fem2.ClearNotTriangulated();
  //fem2.InitCloud();
  //fem2.SetExtrusion(fem1.GetExtrusionDelta(), fem1.GetElementHeight());
  fem1.ViewMesh(true, fem2.GetCloud(), fem2.GetExtrusion(), fem2.GetPose(), 0);

  //return;

  std::cout << "\nTransform pose 2 for simulation, impose a rotation of x degrees around each axis" << std::endl;
  POS pos(fem2.GetEigenNodes(true), fem2.GetPose());
  pos.SetTarget(fem1.GetEigenNodes(true));

  double ang = 5*M_PI/180;
  Eigen::Vector3d axis(1,1,1);

  Eigen::Vector4d imposed_angle_q = pos.QuaternionFromAngleAxis(axis, ang);
  Eigen::Vector3d model_offset(1.0, 1.0, 1.0); 
  pos.Transform(imposed_angle_q, model_offset, 1.0);
  //std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes_k0 = pos.GetPointLayers();
  //fem1.ViewMesh(true, nodes_k0.first, nodes_k0.second, pose_k0, 0);



  std::cout << "\nInitialize FEA object, use conectivity to compute stiffness matrix" << std::endl;
  FEA fea(element, 10000.0, 0.3, true);
  std::vector<Eigen::Vector3d> nodes = fem1.GetEigenNodes(true);
  std::vector<std::vector<unsigned int>> elements = fem1.GetElements();
  fea.MatAssembly(nodes, elements);


  //return;

  std::cout << "\nSet Boundary Conditions" << std::endl;
  std::vector<unsigned int> bc_nodes = fem1.GetExtrusionIndices();
  BoundaryConditions3d bc(fea.NumDof(), &nodes);
  bc.AddEncastreByNodeIds(bc_nodes);
  fea.ApplyBoundaryConditions(bc);









  
  std::cout << "\nInitialize Optimizer" << std::endl;
  LevenbergMarquardt lm(&pos, &fea, 10000, -1.0, 2.0, 0.000001);

  std::vector<Eigen::Vector3d> nodes_k0 = pos.GetTarget();
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k0 = fem1.GetPose();

  std::vector<Eigen::Vector3d> nodes_k1 = pos.GetPoints();
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k1 = pos.GetPose();

  double sE1 = fea.ComputeStrainEnergy(nodes_k0, nodes_k1);
  std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes_k1l = pos.GetPointLayers();

  std::cout << " - Initial Strain energy = " << sE1 << std::endl;
  std::cout << " - Initial Pose = " << pose_k1.first.transpose() << " | " << pose_k1.second.transpose() << std::endl;
  std::cout << " - Target Pose  = " << pose_k0.first.transpose() << " | " << pose_k0.second.transpose() << std::endl;
  fem1.ViewMesh(true, nodes_k1l.first, nodes_k1l.second, pose_k1, 0);
  

  std::cout << "\nOptimize" << std::endl;
  std::pair<double, Eigen::VectorXd> res = lm.Optimize(pose_k1);
  //std::pair<double, Eigen::VectorXd> res = lm.OptimizeCLM(pose_k1);


  std::vector<Eigen::Vector3d> nodes_k2 = pos.GetPoints();
  double sE2 = fea.ComputeStrainEnergy(nodes_k0, nodes_k2);
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k2 = pos.GetPose();
  std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes_k2l = pos.GetPointLayers();

  std::cout << std::endl;
  std::cout << " - Final Strain Energy = " << sE2 << std::endl;
  std::cout << " - Final pose = " << pose_k2.first.transpose() << " " << pose_k2.second.transpose() << std::endl;
  fem1.ViewMesh(true, nodes_k2l.first, nodes_k2l.second, pose_k2, 0);



}


