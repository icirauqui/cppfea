#include "dataset/test_models.cpp"
#include "fea/fea.hpp"
#include "vis/vis.hpp"



void test_c3d6_displacements() {
  double E = 10000.0;
  double nu = 0.3;

  AbaqusC3D6_1 model("../data/c3d6_1/a.input");

  FEA fea("C3D6", E, nu, true);

  std::cout << "\nBuild BoundaryConditions" << std::endl;
  BoundaryConditions3d bc(fea.NumDof(), &model._nodes);
  std::cout << " - Encastre in z = 0" << std::endl;
  bc.EncastreInCoordZ(0.0);
  std::cout << " - Force of magnitude " << model.LoadMagnitude() << " in direction of " << model.LoadDirection() 
            << " on nodes in " << model.LoadDirection() << " = " << model.LoadLocation() << std::endl;
  bc.AddNodalDisplacementByCoordZ(model.LoadLocation(), {0, 0, 0.5});

  std::cout << "\nMatAssembly" << std::endl;
  fea.MatAssembly(model._nodes, model._elements);
  fea.ExportK("../data/" + model.Name() + "/K_a.csv");

  std::cout << "\nApplyBoundaryConditions" << std::endl;
  fea.ApplyBoundaryConditions(bc);
  fea.ExportK("../data/" + model.Name() + "/K_b.csv");
  fea.ExportF("../data/" + model.Name() + "/F.csv");

  std::cout << "\nSolve" << std::endl;
  fea.Solve("CG");
  fea.ExportU("../data/" + model.Name() + "/U.csv");

  std::cout << "\nPostProcess" << std::endl;
  fea.PostProcess(model._nodes, model._elements);
  fea.ReportFEAData("../data/" + model.Name() + "/FEAData.csv");

  std::cout << "Compute deformed positions" << std::endl;
  Eigen::MatrixXd U = fea.U();
  std::vector<Eigen::Vector3d> u;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    u.push_back(U.block(n*3, 0, 3, 1));
  }
  double scale = 100000000000000.0;
  model.ApplyDisplacements(u, {scale, scale, 1.0});

  fea.ExportAll("../data/" + model.Name() + "/KFU.csv");
  fea.ReportNodes("../data/" + model.Name() + "/nodes.csv");

  std::cout << "\nVisualization" << std::endl;
  PCLViewer viewer;
  viewer.AddNodes(model._nodes, "original", Eigen::Vector3d(0.2, 0.2, 1.0));
  viewer.AddNodes(model._nodes_deformed, "deformed", Eigen::Vector3d(0.0, 0.7, 0.0));
  viewer.AddBCs(bc.NodeIds(), bc.Values(), false);
  viewer.AddEdges(model._elements, model.ElementType());
  viewer.Render();
}



void test_c3d6_loads() {
  double E = 10000.0;
  double nu = 0.3;

  AbaqusC3D6_1 model("../data/c3d6_1/a.input");

  FEA fea("C3D6", E, nu, true);

  std::cout << "\nBuild BoundaryConditions" << std::endl;
  BoundaryConditions3d bc(fea.NumDof(), &model._nodes);
  std::cout << " - Encastre in x = 0" << std::endl;
  bc.EncastreInCoordZ(0.0);

  std::cout << "\nBuild Loads" << std::endl;
  Loads3d loads(fea.NumDof(), &model._nodes);
  std::cout << " - Force of magnitude " << model.LoadMagnitude() << " in direction of " << model.LoadDirection() 
            << " on nodes in " << model.LoadDirection() << " = " << model.LoadLocation() << std::endl;
  loads.AddNodalZ(model.LoadLocation(), {0.0, 0.0, -model.LoadMagnitude()});


  std::cout << "\nMatAssembly" << std::endl;
  fea.MatAssembly(model._nodes, model._elements);
  fea.ExportK("../data/" + model.Name() + "/K_a.csv");

  std::cout << "\nApplyBoundaryConditions" << std::endl;
  fea.ApplyBoundaryConditions(bc);
  fea.ApplyLoads(loads);
  fea.ExportK("../data/" + model.Name() + "/K_b.csv");
  fea.ExportF("../data/" + model.Name() + "/F.csv");

  std::cout << "\nSolve" << std::endl;
  fea.Solve("LU");
  fea.ExportU("../data/" + model.Name() + "/U.csv");

  std::cout << "\nPostProcess" << std::endl;
  fea.PostProcess(model._nodes, model._elements);
  fea.ReportFEAData("../data/" + model.Name() + "/FEAData.csv");

  std::cout << "Compute deformed positions" << std::endl;
  Eigen::MatrixXd U = fea.U();
  std::vector<Eigen::Vector3d> u;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    u.push_back(U.block(n*3, 0, 3, 1));
  }
  double scale = 1000.0;
  model.ApplyDisplacements(u, scale);

  fea.ExportAll("../data/" + model.Name() + "/KFU.csv");
  fea.ExportSystem("../data/" + model.Name() + "/K_U_F.csv");
  fea.ReportNodes("../data/" + model.Name() + "/nodes.csv");

  std::cout << "\nVisualization" << std::endl;
  PCLViewer viewer;
  viewer.AddNodes(model._nodes, "original", Eigen::Vector3d(0.2, 0.2, 1.0));
  viewer.AddNodes(model._nodes_deformed, "deformed", Eigen::Vector3d(0.0, 0.7, 0.0));
  viewer.AddBCs(bc.NodeIds(), bc.Values(), false);
  viewer.AddBCs(loads.NodeIds(), loads.Values(), true);
  viewer.AddEdges(model._elements, model.ElementType());
  viewer.Render();
}



void test_c3d8_displacements() {
  double E = 10000.0;
  double nu = 0.3;

  AbaqusC3D8_2 model;

  FEA fea("C3D8", E, nu, true);

  std::cout << "\nBuild BoundaryConditions" << std::endl;
  BoundaryConditions3d bc(fea.NumDof(), &model._nodes);
  std::cout << " - Encastre in z = 0" << std::endl;
  bc.EncastreInCoordZ(0.0);
  std::cout << " - Displacement of magnitude 1 in direction x on nodes in x = 5.0" << std::endl;
  bc.AddNodalDisplacementByCoordZ(model.LoadLocation(), {0, 0, 0.5});

  std::cout << "\nMatAssembly" << std::endl;
  fea.MatAssembly(model._nodes, model._elements);
  fea.ExportK("../data/" + model.Name() + "/K_a.csv");

  std::cout << "\nApplyBoundaryConditions" << std::endl;
  fea.ApplyBoundaryConditions(bc);
  fea.ExportK("../data/" + model.Name() + "/K_b.csv");
  fea.ExportF("../data/" + model.Name() + "/F.csv");

  std::cout << "\nSolve" << std::endl;
  fea.Solve("CG");
  fea.ExportU("../data/" + model.Name() + "/U.csv");

  std::cout << "\nPostProcess" << std::endl;
  fea.PostProcess(model._nodes, model._elements);
  fea.ReportFEAData("../data/" + model.Name() + "/FEAData.csv");

  std::cout << "Compute deformed positions" << std::endl;
  Eigen::MatrixXd U = fea.U();
  std::vector<Eigen::Vector3d> u;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    u.push_back(U.block(n*3, 0, 3, 1));
  }
  double scale = 1.0;
  model.ApplyDisplacements(u, scale);

  fea.ExportAll("../data/" + model.Name() + "/KFU.csv");
  fea.ReportNodes("../data/" + model.Name() + "/nodes.csv");

  std::cout << "\nVisualization" << std::endl;
  PCLViewer viewer;
  viewer.AddNodes(model._nodes, "original", Eigen::Vector3d(0.2, 0.2, 1.0));
  viewer.AddNodes(model._nodes_deformed, "deformed", Eigen::Vector3d(0.0, 0.7, 0.0));
  viewer.AddBCs(bc.NodeIds(), bc.Values(), false);
  viewer.AddEdges(model._elements, model.ElementType());
  viewer.Render();
}



void test_c3d8_loads() {
  double E = 10000.0;
  double nu = 0.3;

  AbaqusC3D8_2 model;

  FEA fea("C3D8", E, nu, true);

  std::cout << "\nBuild BoundaryConditions" << std::endl;
  BoundaryConditions3d bc(fea.NumDof(), &model._nodes);
  std::cout << " - Encastre in z = 0" << std::endl;
  bc.EncastreInCoordZ(0.0);

  std::cout << "\nBuild Loads" << std::endl;
  Loads3d loads(fea.NumDof(), &model._nodes);
  std::cout << " - Force of magnitude 1 in direction of z on nodes in z = " << model.LoadLocation() << std::endl;
  loads.AddNodalZ(model.LoadLocation(), {0.0, 0.0, 0.5});


  std::cout << "\nMatAssembly" << std::endl;
  fea.MatAssembly(model._nodes, model._elements);
  fea.ExportK("../data/" + model.Name() + "/K_a.csv");

  std::cout << "\nApplyBoundaryConditions" << std::endl;
  fea.ApplyBoundaryConditions(bc);
  fea.ApplyLoads(loads);
  fea.ExportK("../data/" + model.Name() + "/K_b.csv");
  fea.ExportF("../data/" + model.Name() + "/F.csv");

  std::cout << "\nSolve" << std::endl;
  fea.Solve("CG");
  fea.ExportU("../data/" + model.Name() + "/U.csv");

  std::cout << "\nPostProcess" << std::endl;
  fea.PostProcess(model._nodes, model._elements);
  fea.ReportFEAData("../data/" + model.Name() + "/FEAData.csv");

  std::cout << "Compute deformed positions" << std::endl;
  Eigen::MatrixXd U = fea.U();
  std::vector<Eigen::Vector3d> u;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    u.push_back(U.block(n*3, 0, 3, 1));
  }
  double scale = 100.0;
  model.ApplyDisplacements(u, scale);

  fea.ExportAll("../data/" + model.Name() + "/KFU.csv");
  fea.ReportNodes("../data/" + model.Name() + "/nodes.csv");

  std::cout << "\nVisualization" << std::endl;
  PCLViewer viewer;
  viewer.AddNodes(model._nodes, "original", Eigen::Vector3d(0.2, 0.2, 1.0));
  viewer.AddNodes(model._nodes_deformed, "deformed", Eigen::Vector3d(0.0, 0.7, 0.0));
  viewer.AddBCs(bc.NodeIds(), bc.Values(), false);
  viewer.AddBCs(loads.NodeIds(), loads.Values(), true);
  viewer.AddEdges(model._elements, model.ElementType());
  viewer.Render();
}



void test_c2d4_displacements() {
  double E = 10000.0;
  double nu = 0.3;

  AbaqusC2D4_1 model("../py/tn5.input");

  FEA fea("C2D4", 100.0, 0.3, true);

  std::cout << "\nBuild BoundaryConditions" << std::endl;
  BoundaryConditions2d bc(fea.NumDof(), &model._nodes);
  std::cout << " - Encastre in x = 0" << std::endl;
  bc.EncastreInCoordX(0.0);
  std::cout << " - Displacement of magnitude {1,1} in direction x on nodes in x = 5.0" << std::endl;
  bc.AddNodalDisplacementByCoordX(5.0, {1.0, 1.0});
 
  std::cout << "\nMatAssembly" << std::endl;
  fea.MatAssembly(model._nodes, model._elements);
  fea.ExportK("../data/" + model.Name() + "/K_a.csv");

  std::cout << "\nApplyBoundaryConditions" << std::endl;
  fea.ApplyBoundaryConditions(bc);
  fea.ExportK("../data/" + model.Name() + "/K_b.csv");
  fea.ExportF("../data/" + model.Name() + "/F.csv");

  std::cout << "\nSolve" << std::endl;
  fea.Solve("LU");
  fea.ExportU("../data/" + model.Name() + "/U.csv");

  std::cout << "\nPostProcess" << std::endl;
  fea.PostProcess(model._nodes, model._elements);
  fea.ReportFEAData("../data/" + model.Name() + "/FEAData.csv");

  std::cout << "Compute deformed positions" << std::endl;
  Eigen::MatrixXd U = fea.U();
  std::vector<Eigen::Vector2d> u;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    u.push_back(U.block(n*2, 0, 2, 1));
  }
  double scale = 1.0;
  model.ApplyDisplacements(u, scale);

  fea.ExportAll("../data/" + model.Name() + "/KFU.csv");
  fea.ReportNodes("../data/" + model.Name() + "/nodes.csv");

  std::vector<Eigen::Vector3d> nodes_vis, nodes_deformed_vis;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    nodes_vis.push_back(Eigen::Vector3d(model._nodes[n][0], model._nodes[n][1], 0.0));
    nodes_deformed_vis.push_back(Eigen::Vector3d(model._nodes_deformed[n][0], model._nodes_deformed[n][1], 0.0));
  }


  std::cout << "\nVisualization" << std::endl;
  PCLViewer viewer;
  viewer.AddNodes(nodes_vis, "original", Eigen::Vector3d(0.2, 0.2, 1.0));
  viewer.AddNodes(nodes_deformed_vis, "deformed", Eigen::Vector3d(0.0, 0.7, 0.0));
  viewer.AddBCs(bc.NodeIds(), bc.Values(), false);
  viewer.AddEdges(model._elements, model.ElementType());
  viewer.Render();
}



void test_c2d4_loads() {
  double E = 10000.0;
  double nu = 0.3;

  AbaqusC2D4_1 model("../py/tn5.input");

  FEA fea("C2D4", 100.0, 0.3, true);

  std::cout << "\nBuild BoundaryConditions" << std::endl;
  BoundaryConditions2d bc(fea.NumDof(), &model._nodes);
  std::cout << " - Encastre in x = 0" << std::endl;
  bc.EncastreInCoordX(0.0);

  std::cout << "\nBuild Loads" << std::endl;
  Loads2d loads(fea.NumDof(), &model._nodes);
  std::cout << " - Force of magnitude 1 in direction x on nodes in x = 115.0" << std::endl;
  loads.AddNodalX(5.0, {1.0, 1.0});
  
  std::cout << "\nMatAssembly" << std::endl;
  fea.MatAssembly(model._nodes, model._elements);
  fea.ExportK("../data/" + model.Name() + "/K_a.csv");

  std::cout << "\nApplyBoundaryConditions" << std::endl;
  fea.ApplyBoundaryConditions(bc);
  fea.ApplyLoads(loads);
  fea.ExportK("../data/" + model.Name() + "/K_b.csv");
  fea.ExportF("../data/" + model.Name() + "/F.csv");

  std::cout << "\nSolve" << std::endl;
  fea.Solve("LU");
  fea.ExportU("../data/" + model.Name() + "/U.csv");

  std::cout << "\nPostProcess" << std::endl;
  fea.PostProcess(model._nodes, model._elements);
  fea.ReportFEAData("../data/" + model.Name() + "/FEAData.csv");

  std::cout << "Compute deformed positions" << std::endl;
  Eigen::MatrixXd U = fea.U();
  std::vector<Eigen::Vector2d> u;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    u.push_back(U.block(n*2, 0, 2, 1));
  }
  double scale = 1.0;
  model.ApplyDisplacements(u, scale);

  fea.ExportAll("../data/" + model.Name() + "/KFU.csv");
  fea.ReportNodes("../data/" + model.Name() + "/nodes.csv");

  std::vector<Eigen::Vector3d> nodes_vis, nodes_deformed_vis;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    nodes_vis.push_back(Eigen::Vector3d(model._nodes[n][0], model._nodes[n][1], 0.0));
    nodes_deformed_vis.push_back(Eigen::Vector3d(model._nodes_deformed[n][0], model._nodes_deformed[n][1], 0.0));
  }

  std::cout << "\nVisualization" << std::endl;
  PCLViewer viewer;
  viewer.AddNodes(nodes_vis, "original", Eigen::Vector3d(0.2, 0.2, 1.0));
  viewer.AddNodes(nodes_deformed_vis, "deformed", Eigen::Vector3d(0.0, 0.7, 0.0));
  viewer.AddBCs(bc.NodeIds(), bc.Values(), false);
  viewer.AddBCs(loads.NodeIds(), loads.Values(), true);
  viewer.AddEdges(model._elements, model.ElementType());
  viewer.Render();
}
