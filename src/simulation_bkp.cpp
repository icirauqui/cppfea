#include <iostream>
#include <vector>
#include <math.h>

#include "fea/fem.hpp"
#include "fea/fea.hpp"
#include "fea/pos.hpp"
#include "dataset/dataset.hpp"
#include "nlo/levenberg_marquardt.hpp"

#include "test_models.cpp"

#include "fea/elts/c3d8.hpp"

#include <chrono>

// Parameters
int image_id = 0;
std::string element = "C3D8";
float E = 10000.0;
float nu = 0.495;
float depth = 1.0;
float fg = 0.577350269;
float Klarge = 100000000.0;

int view_quick = 1;
int view_slow = 3;

float noise_multiplier = 0.1;
Eigen::Vector3d model_offset(0.2, 0.2, 0.2); 

std::pair<Eigen::Vector4d, Eigen::Vector3d> ApproximatePose(std::vector<Eigen::Vector3d> pts) {
  Eigen::Vector3d centroid(0.0, 0.0, 0.0);
  for (int i = 0; i < pts.size(); i++) {
    centroid += pts[i];
  }
  centroid /= pts.size();

  double dist = sqrt(pow(centroid(0), 2) + pow(centroid(1), 2) + pow(centroid(2), 2));
  Eigen::Vector3d tvec = centroid + Eigen::Vector3d(0, 0, dist/2);

  // Compute rotation, from t to centroid, as a quaternion
  Eigen::Vector3d direction = centroid - tvec;
  direction.normalize();
  Eigen::Quaterniond quaternion;
  quaternion.setFromTwoVectors(-Eigen::Vector3d::UnitX(), direction);
  Eigen::Vector4d qvec = quaternion.coeffs();

  return std::make_pair(qvec, tvec);
}

void test_fea() {
  // Load data
  Dataset ds("../data", element);
  std::vector<std::vector<float>> vpts = ds.points();
  std::vector<std::vector<int>>   vElts = ds.elements();
  std::vector<std::vector<float>> vvF = ds.forces();
  std::vector<std::vector<int>>   vDir = ds.dirichlet();

  // Test new class
  FEA fea(0, element, E, nu, depth, fg, false);
  
  auto start1 = std::chrono::high_resolution_clock::now();
  fea.MatAssembly(vpts, vElts);
  fea.SetForces(vvF);
  fea.ImposeDirichletEncastre(vDir, Klarge);
  fea.ComputeDisplacements();
  fea.ComputeStrainEnergy();
  auto stop1 = std::chrono::high_resolution_clock::now();

  std::cout << " - Strain energy = " << fea.StrainEnergy() << std::endl;
  std::cout << "   Time fea1 = " << std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1).count() << std::endl;
}

void test_fe(bool optimizer = false) {

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
  //fem1.ComputeExtrusion();
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


  if (optimizer == false) {
    // 6.1 Simulate 5 steps of translation and rotation, model 2 will move towards model 1.
    int num_steps = 10;
    int div_steps = 5;
    Eigen::Vector3d step = (model_offset) / div_steps;
    Eigen::Vector4d step_rot = pos.QuaternionFromAngleAxis(axis, -ang/div_steps);

    for (unsigned int s=0; s<num_steps; s++) {
      std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k1 = pos.GetPose();
      std::vector<Eigen::Vector3d> nodes_k1 = pos.GetPoints();
      std::vector<Eigen::Vector3d> nodes_k0 = pos.GetTarget();

      pos.Transform(step_rot, -step, 1.0);

      std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k = pos.GetPose();

      std::vector<Eigen::Vector3d> nodes_k = pos.GetPoints();
      std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> nodes_k_layers = pos.GetPointLayers();

      double sE = fea.ComputeStrainEnergy(nodes_k0, nodes_k);

      //std::cout << std::endl;
      //std::cout << "Step " << s << std::endl;
      //std::cout << "  Translation = " << pose_k1.first.transpose() << " -> " << pose_k.first.transpose() << std::endl;
      //std::cout << "  Orientation = " << pose_k1.second.transpose() << " -> " << pose_k.second.transpose() << std::endl;
      std::cout << "  Strain energy = " << sE << std::endl;

      //std::cout << "A" << std::endl;
      //for (unsigned int i=0; i<20; i++) {
      //  std::cout << new_nodes_front[i].transpose()
      //            << "\t | \t" 
      //            << nodes_k_layers.first[i].transpose()
      //            << "\t"
      //            << std::endl;
      //}

      fem1.ViewMesh(true, nodes_k_layers.first, nodes_k_layers.second, pose1, pose_k, view_quick);
    }
  } else {


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




}

void eval_fea() {
  Dataset ds("../data_test_1", element);
  std::vector<std::vector<float>> vpts = ds.points();

  FEM fem1(element);
  FEM fem2(element);

  for (int i = 0; i < vpts.size(); i++) {
    Eigen::Vector3d pt1(vpts[i][0], vpts[i][1], vpts[i][2]);
    Eigen::Vector3d pt2(vpts[i][0], vpts[i][1], vpts[i][2]);
    //for (unsigned int j=0; j<3; j++) {
    //  //pt1(j) += noise_multiplier * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    //  pt2(j) += noise_multiplier * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    //}
    //if (i==vpts.size()-1) {
    //  pt2(0) += 0.5;
    //  pt2(1) += 0.5;
    //}
    pt2 += model_offset;
    fem1.AddPoint(pt1);
    fem2.AddPoint(pt2);
  }

  fem1.Compute(true);
  fem2.InitCloud();
  //fem1.ComputeExtrusion();
  fem2.SetExtrusion(fem1.GetExtrusionDelta(), fem1.GetElementHeight());

  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose1 = ApproximatePose(fem1.GetEigenNodes());
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2 = ApproximatePose(fem2.GetEigenNodes());


  POS pos(fem2.GetNodes(), pose2);
  pos.SetTarget(fem1.GetNodes());

  std::vector<Eigen::Vector3d> nodes_k0 = pos.GetTarget();
  std::vector<Eigen::Vector3d> nodes_k1 = pos.GetPoints();

  FEA fea(0, element, E, nu, depth, fg, false);
  std::vector<std::vector<float>> nodes = fem1.GetNodes();
  std::vector<std::vector<int>> elements = fem1.GetElements();
  fea.MatAssembly(nodes, elements);
  fea.EncastreBackLayer();
  fea.Eigenvalues_K();

  double sE = fea.ComputeStrainEnergy(nodes_k0, nodes_k1);
  std::cout << "Strain energy = " << sE << std::endl;

  fem1.ViewMesh(true, fem2.GetCloud(), fem2.GetExtrusion(), pose1, pose2, 10);

}

void test_c3d6() {
  // Define material properties
  double E = 3500; // Young's modulus (Pa)
  double nu = 0.495;  // Poisson's ratio

    // Define the nodes of the triangular prism
    std::array<Eigen::Vector3d, 6> nodes = {
        Eigen::Vector3d(0, 0, 0),
        Eigen::Vector3d(1, 0, 0),
        Eigen::Vector3d(0, 1, 0),
        Eigen::Vector3d(0, 0, 0.5),
        Eigen::Vector3d(1, 0, 0.5),
        Eigen::Vector3d(0, 1, 0.5)
    };

    // Compute the stiffness matrix
    Eigen::MatrixXd K = c3d6::computeStiffnessMatrix(nodes, E, nu);

    // Output the stiffness matrix
    std::cout << "Stiffness Matrix: \n" << K << std::endl;

}

void test_c3d8() {
  // Define material properties
  double E = 3500; // Young's modulus (Pa)
  double nu = 0.495;  // Poisson's ratio

    // Define the nodes of the triangular prism
    std::array<Eigen::Vector3d, 8> nodes = {
        Eigen::Vector3d(0, 0, 0),
        Eigen::Vector3d(1, 0, 0),
        Eigen::Vector3d(0, 1, 0),
        Eigen::Vector3d(1, 1, 0),
        Eigen::Vector3d(0, 0, 1),
        Eigen::Vector3d(1, 0, 1),
        Eigen::Vector3d(0, 1, 1),
        Eigen::Vector3d(1, 1, 1)
    };

    // Compute the stiffness matrix
    Eigen::MatrixXd K = c3d8::computeStiffnessMatrix(nodes, E, nu);

    // Output the stiffness matrix
    std::cout << "Stiffness Matrix: \n" << K << std::endl;
}



int main(int argc, char** argv) {
  AbaqusC3D8_1 model;
  C3D8 fea = C3D8();

  Eigen::MatrixXd K = c3d8::matAssembly(model._nodes, model._elements, E, nu);
  std::cout << "Stiffness Matrix: \n" << K << std::endl;
  
  return 0;
}
