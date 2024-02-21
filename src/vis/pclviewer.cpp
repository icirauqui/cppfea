#include "pclviewer.h"

PCLViewer::PCLViewer(unsigned int nPointClouds) {
  images_.resize(nPointClouds);
  meshes_.resize(nPointClouds);

  clouds_images_.resize(nPointClouds);
  clouds_p3D_.resize(nPointClouds);
  clouds_p3Dd_.resize(nPointClouds);
  clouds_mesh_.resize(nPointClouds);

  clouds_p3D_colors_.resize(nPointClouds);
  clouds_p3Dd_colors_.resize(nPointClouds);

  for (unsigned int i = 0; i < nPointClouds; i++) {
    clouds_images_[i] = nullptr;
    clouds_p3D_[i] = nullptr;
    clouds_p3Dd_[i] = nullptr;
    clouds_mesh_[i] = nullptr;
  }
}

PCLViewer::~PCLViewer() {}

void PCLViewer::initializeViewer() {
  viewer_->setBackgroundColor(1.0, 1.0, 1.0);
  viewer_->addCoordinateSystem(1.0);
  viewer_->initCameraParameters();
  viewer_->setSize(1070, 820);
}

void PCLViewer::registerWindowSize(int width, int height) {
  win_width_1_ = width;
  win_height_1_ = height;
}

void PCLViewer::addImages(std::vector<Image*> images, int modelId) {
  clearImages(modelId);

  for (auto im : images)
    images_[modelId].push_back(new PCLImage(im->Qvec(), im->Tvec(), modelId));

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ptr(
      new pcl::PointCloud<pcl::PointXYZ>);

  float l = 0.1;
  float l1 = 1;  // Director vector length
  std::vector<int> indices = {0, 1, 2, 0, 4, 3, 0, 1, 4, 3, 2, 0};

  // Create image pyramid for visualization (5 points)
  for (unsigned int i = 0; i < images_[modelId].size(); i++) {
    std::vector<pcl::PointXYZ> points(5, pcl::PointXYZ());

    points[0].x = 0.0;
    points[0].y = 0.0;
    points[0].z = 0.0;

    points[1].x = -l;
    points[1].y = +0.75 * l;
    points[1].z = +l;

    points[2].x = +l;
    points[2].y = +0.75 * l;
    points[2].z = +l;

    points[3].x = +l;
    points[3].y = -0.75 * l;
    points[3].z = +l;

    points[4].x = -l;
    points[4].y = -0.75 * l;
    points[4].z = +l;

    // Rotate points of camera frame
    for (unsigned int j = 0; j < points.size(); j++) {
      Eigen::Vector3d pte(points[j].x, points[j].y, points[j].z);
      Eigen::Vector3d pte2 =
          QuaternionRotatePoint(images_[modelId][i]->get_q(), pte);
      points[j].x = pte2(0);
      points[j].y = pte2(1);
      points[j].z = pte2(2);
    }

    // Translate points of rotated camera frame
    for (unsigned int j = 0; j < points.size(); j++) {
      points[j].x += images_[modelId][i]->get_t(0);
      points[j].y += images_[modelId][i]->get_t(1);
      points[j].z += images_[modelId][i]->get_t(2);
    }

    // Draw rotated camera frame
    for (unsigned int j = 0; j < indices.size(); j++)
      images_[modelId][i]->poly_->points.push_back(points[indices[j]]);

    // Set director vector base point in origin
    std::vector<pcl::PointXYZ> points_vdir(2, pcl::PointXYZ());
    points_vdir[0].x = 0.0;
    points_vdir[0].y = 0.0;
    points_vdir[0].z = 0.0;

    // Set end point of director vector
    Eigen::Vector3d vdir =
        QuaternionToDirectorVector(images_[modelId][i]->get_q());
    points_vdir[1].x = l1 * vdir(0);
    points_vdir[1].y = l1 * vdir(1);
    points_vdir[1].z = l1 * vdir(2);

    // Translate director vector
    for (unsigned int j = 0; j < points_vdir.size(); j++) {
      points_vdir[j].x += images_[modelId][i]->get_t(0);
      points_vdir[j].y += images_[modelId][i]->get_t(1);
      points_vdir[j].z += images_[modelId][i]->get_t(2);
    }
  }

  cloud_ptr->width = cloud_ptr->size();
  cloud_ptr->height = 1;

  clouds_images_[modelId] = cloud_ptr;
}

void PCLViewer::addImages2(std::vector<Image*> images, int modelId) {
  clearImages(modelId);

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ptr(
      new pcl::PointCloud<pcl::PointXYZ>);

  std::vector<int> indices = {0, 1, 2, 0, 4, 3, 0, 1, 4, 3, 2, 0};

  float image_width = 0.1;
  float image_height = 0.75 * image_width;
  float focal_length = image_width;

  for (auto im : images) {
    PCLImage* pclim = new PCLImage(im->Qvec(), im->Tvec(), modelId);

    const Eigen::Matrix<float, 3, 4> inv_proj_matrix =
        im->InverseProjectionMatrix().cast<float>();

    // Projection center, top-left, top-right, bottom-right, bottom-left
    // corners.

    const Eigen::Vector3f pc = inv_proj_matrix.rightCols<1>();
    const Eigen::Vector3f tl =
        inv_proj_matrix *
        Eigen::Vector4f(-image_width, image_height, focal_length, 1);
    const Eigen::Vector3f tr =
        inv_proj_matrix *
        Eigen::Vector4f(image_width, image_height, focal_length, 1);
    const Eigen::Vector3f br =
        inv_proj_matrix *
        Eigen::Vector4f(image_width, -image_height, focal_length, 1);
    const Eigen::Vector3f bl =
        inv_proj_matrix *
        Eigen::Vector4f(-image_width, -image_height, focal_length, 1);

    std::vector<pcl::PointXYZ> points(5, pcl::PointXYZ());

    points[0].x = pc(0);
    points[0].y = pc(1);
    points[0].z = pc(2);

    points[1].x = tl(0);
    points[1].y = tl(1);
    points[1].z = tl(2);

    points[2].x = tr(0);
    points[2].y = tr(1);
    points[2].z = tr(2);

    points[3].x = br(0);
    points[3].y = br(1);
    points[3].z = br(2);

    points[4].x = bl(0);
    points[4].y = bl(1);
    points[4].z = bl(2);

    // Draw rotated camera frame
    for (unsigned int j = 0; j < indices.size(); j++)
      pclim->poly_->points.push_back(points[indices[j]]);

    images_[modelId].push_back(pclim);
  }

  cloud_ptr->width = cloud_ptr->size();
  cloud_ptr->height = 1;

  clouds_images_[modelId] = cloud_ptr;
}

void PCLViewer::clearImages(int modelId) {
  for (unsigned int i = 0; i < images_[modelId].size(); i++) {
    delete images_[modelId][i];
    images_[modelId][i] = nullptr;
  }
  images_[modelId].clear();
}

void PCLViewer::addPointCloud(std::vector<Eigen::Vector3d> points,
                              int cloudId) {
  if (clouds_p3D_[cloudId]) delPointCloud(cloudId);

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ptr(
      new pcl::PointCloud<pcl::PointXYZ>);

  for (Eigen::Vector3d point : points) {
    pcl::PointXYZ pt;
    pt.x = point(0);
    pt.y = point(1);
    pt.z = point(2);
    cloud_ptr->points.push_back(pt);
  }

  cloud_ptr->width = cloud_ptr->size();
  cloud_ptr->height = 1;

  clouds_p3D_[cloudId] = cloud_ptr;
}

void PCLViewer::addFusedPointCloud(std::vector<Eigen::Vector3d> points,
                                   std::vector<Eigen::Vector3d> colors,
                                   std::vector<std::pair<double, double>> points_range,
                                   std::vector<std::vector<Eigen::Vector3d>> xyzs,
                                   int cloudId) {
  if (clouds_p3D_[cloudId]) delPointCloud(cloudId);

  clouds_p3D_range_.clear();
  clouds_p3D_xyzs_.clear();

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ptr(
      new pcl::PointCloud<pcl::PointXYZ>);

  for (unsigned int i=0; i<points.size(); i++) {
    pcl::PointXYZ pt;
    pt.x = points[i](0);
    pt.y = points[i](1);
    pt.z = points[i](2);
    cloud_ptr->points.push_back(pt);

    clouds_p3D_range_.push_back(points_range[i]);
    
    //std::cout << "Size in pclviewer = " << xyzs[i].size() << std::endl;
    std::vector<Eigen::Vector3d> xyzs_i;
    for (unsigned int j=0; j<xyzs[i].size(); j++) {
      xyzs_i.push_back(xyzs[i][j]);
    }
    clouds_p3D_xyzs_.push_back(xyzs_i);
  }

  cloud_ptr->width = cloud_ptr->size();
  cloud_ptr->height = 1;

  clouds_p3D_[cloudId] = cloud_ptr;

  clouds_p3D_colors_[cloudId] = colors;
}

void PCLViewer::delPointCloud(int cloudId) { 
  clouds_p3D_[cloudId] = nullptr; 
  clouds_p3D_colors_[cloudId].clear();
}

void PCLViewer::addPointCloudDef(std::vector<Eigen::Vector3d> points,
                                 int cloudId) {
  if (clouds_p3Dd_[cloudId]) delPointCloudDef(cloudId);

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ptr(
      new pcl::PointCloud<pcl::PointXYZ>);

  for (Eigen::Vector3d point : points) {
    pcl::PointXYZ pt;
    pt.x = point(0);
    pt.y = point(1);
    pt.z = point(2);
    cloud_ptr->points.push_back(pt);
  }

  cloud_ptr->width = cloud_ptr->size();
  cloud_ptr->height = 1;

  clouds_p3Dd_[cloudId] = cloud_ptr;
}

void PCLViewer::addFusedPointCloudDef(std::vector<Eigen::Vector3d> points,
                                      std::vector<std::pair<double, double>> points_range,
                                      std::vector<Eigen::Vector3d> xyzs,
                                      int cloudId) {
  if (clouds_p3Dd_[cloudId]) delPointCloudDef(cloudId);

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ptr(
      new pcl::PointCloud<pcl::PointXYZ>);

  for (unsigned int i=0; i<points.size(); i++) {
    pcl::PointXYZ pt;
    pt.x = points[i](0);
    pt.y = points[i](1);
    pt.z = points[i](2);
    cloud_ptr->points.push_back(pt);

    clouds_p3Dd_range_.push_back(points_range[i]);
  }

  cloud_ptr->width = cloud_ptr->size();
  cloud_ptr->height = 1;

  clouds_p3Dd_[cloudId] = cloud_ptr;
}

void PCLViewer::delPointCloudDef(int cloudId) {
  clouds_p3Dd_[cloudId] = nullptr;
}

void PCLViewer::AddMesh(std::vector<std::vector<int>> mesh, int cloudId) {
  if (clouds_mesh_[cloudId]) DelMesh(cloudId);

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ptr(
      new pcl::PointCloud<pcl::PointXYZ>);

  std::vector<int> indices = {0, 1, 2, 0};

  for (auto tri : mesh) {
    PCLTriangle* pcltri = new PCLTriangle(cloudId);

    std::vector<pcl::PointXYZ> points(3, pcl::PointXYZ());

    for (unsigned int i = 0; i < tri.size(); i++) {
      points[i].x = clouds_p3D_[cloudId]->points[tri[i]].x;
      points[i].y = clouds_p3D_[cloudId]->points[tri[i]].y;
      points[i].z = clouds_p3D_[cloudId]->points[tri[i]].z;

      pcl::PointXYZ pt;
      pt.x = clouds_p3D_[cloudId]->points[tri[i]].x;
      pt.y = clouds_p3D_[cloudId]->points[tri[i]].y;
      pt.z = clouds_p3D_[cloudId]->points[tri[i]].z;
      cloud_ptr->points.push_back(pt);
    }

    for (unsigned int i = 0; i < indices.size(); i++)
      pcltri->poly_->points.push_back(points[indices[i]]);

    meshes_[cloudId].push_back(pcltri);
  }

  cloud_ptr->width = cloud_ptr->size();
  cloud_ptr->height = 1;

  clouds_mesh_[cloudId] = cloud_ptr;
}

void PCLViewer::DelMesh(int cloudId) {
  for (unsigned int i = 0; i < meshes_[cloudId].size(); i++) {
    delete meshes_[cloudId][i];
    meshes_[cloudId][i] = nullptr;
  }
  meshes_[cloudId].clear();
}

void PCLViewer::render() {
  viewer_->removeAllPointClouds();
  viewer_->removeAllShapes();

  if (display_meshes_) {
    render_meshes(viewer_);
  } else if (display_compare_) {
    render_meshes(viewer_);
  } else {
    if (display_images_1_ || display_images_2_ || display_images_agg_) {
      render_images(viewer_);
    } 
    if (display_points3D_1_ || display_points3D_2_ || display_points3D_agg_) {
      render_clouds_p3D(viewer_);
    }
    if (display_points3Dd_1_ || display_points3Dd_2_ || display_points3Dd_agg_) {
      render_clouds_p3Dd(viewer_);
    }
  }

  viewer_->spin();
}

void PCLViewer::render_images(pcl::visualization::PCLVisualizer::Ptr viewer) {
  for (unsigned int i = 0; i < clouds_images_.size(); i++) {
    if (!clouds_images_[i] || 
        (i == 0 && (!display_images_1_ || display_agg_)) ||
        (i == 1 && (!display_images_2_ || display_agg_)) ||
        (i == 2 && (!display_images_agg_ || !display_agg_)))
      continue;

    std::string id = "im" + std::to_string(i);

    pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud = clouds_images_[i];
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZ> rgb(
        cloud);
    viewer->addPointCloud<pcl::PointXYZ>(cloud, id);
    // if (i == 0){
    //     viewer->setPointCloudRenderingProperties
    //     (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, id);
    //     viewer->setPointCloudRenderingProperties
    //     (pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, id);
    // }
    // else{
    // if (i != 0){
    //     viewer->setPointCloudRenderingProperties
    //     (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, id);
    //     viewer->setPointCloudRenderingProperties
    //     (pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 0.14, id);
    // }

    for (unsigned int j = 0; j < images_[i].size(); j++) {
      std::string id_poly =
          "poly_" + std::to_string(i) + "_" + std::to_string(j);
      pcl::PointCloud<pcl::PointXYZ>::ConstPtr polygon = images_[i][j]->poly_;
      if (i == 0 || i == 2)
        viewer->addPolygon<pcl::PointXYZ>(polygon, 0.0, 1.0, 0.0, id_poly, 0);
      else
        viewer->addPolygon<pcl::PointXYZ>(polygon, 1.0, 1.0, 0.14, id_poly, 0);
    }
  }
}

void PCLViewer::render_clouds_p3D(
    pcl::visualization::PCLVisualizer::Ptr viewer) {
  for (unsigned int i = 0; i < clouds_p3D_.size(); i++) {
    if (!clouds_p3D_[i] || 
        (i == 0 && (!display_points3D_1_ || display_agg_)) ||
        (i == 1 && (!display_points3D_2_ || display_agg_)) ||
        (i == 2 && (!display_points3D_agg_ || !display_agg_)))
      continue;

    std::string id = "p3D" + std::to_string(i);

    pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud = clouds_p3D_[i];
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZ> rgb(cloud);

    if (i==2) {
      std::vector<Eigen::Vector3d> colors = clouds_p3D_colors_[i];

      pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_rgb(new pcl::PointCloud<pcl::PointXYZRGB>);
      cloud_rgb->width = cloud->width;
      cloud_rgb->height = cloud->height;
      cloud_rgb->points.resize(cloud_rgb->width * cloud_rgb->height);

      for (unsigned int j = 0; j < cloud->size(); j++) {
        cloud_rgb->points[j].x = cloud->points[j].x;
        cloud_rgb->points[j].y = cloud->points[j].y;
        cloud_rgb->points[j].z = cloud->points[j].z;
        cloud_rgb->points[j].r = colors[j](0);
        cloud_rgb->points[j].g = colors[j](1);
        cloud_rgb->points[j].b = colors[j](2);
      }

      viewer->addPointCloud<pcl::PointXYZRGB>(cloud_rgb, id);
    } else {
      viewer->addPointCloud<pcl::PointXYZ>(cloud, id);
    }

    if (i == 0 || i == 2) {
      viewer->setPointCloudRenderingProperties(
          pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, id);

      if (i != 2) {
        viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 1.0, id);
      }

      if (i==2) { // Add spheres with the same color as the points
        if (display_agg_spheres_ || display_agg_trajectories_) {
          for (unsigned int j = 0; j < cloud->size(); j++) {
            if (display_agg_spheres_) {
              std::string id_sphere_avg = "sphere_" + std::to_string(i) + "_" + std::to_string(j) + "_avg";
              std::string id_sphere_max = "sphere_" + std::to_string(i) + "_" + std::to_string(j) + "_max";
              
              // Add sphere with radious clouds_p3D_range_[j].first
              viewer->addSphere(cloud->points[j], clouds_p3D_range_[j].first, 0.25, 0.53, 0.56, id_sphere_avg, 0);
              viewer->addSphere(cloud->points[j], clouds_p3D_range_[j].second, 0.25, 0.53, 0.56, id_sphere_max, 0);

              // Set transparency: avg = 0.3, max = 0.1
              viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, sphere_transparency_avg_, id_sphere_avg);
              viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, sphere_transparency_max_, id_sphere_max);
            }

            // Draw trajectories from clouds_p3D_xyzs_
            if (display_agg_trajectories_) {
              if (clouds_p3D_xyzs_[j].size() < 2) {
                //std::cout << " Display trajectory " << j << " with NA" << std::endl;
                continue;
              }
              //std::cout << " Display trajectory " << j << " with " << clouds_p3D_xyzs_[j].size() << " points" << std::endl;
              for (unsigned int k = 0; k < clouds_p3D_xyzs_[j].size() - 1; k++) {
                std::string id_line = "line_traj_" + std::to_string(j) + "_" + std::to_string(j) + "_" + std::to_string(k);
                pcl::PointXYZ p1(clouds_p3D_xyzs_[j][k](0), clouds_p3D_xyzs_[j][k](1), clouds_p3D_xyzs_[j][k](2));
                pcl::PointXYZ p2(clouds_p3D_xyzs_[j][k+1](0), clouds_p3D_xyzs_[j][k+1](1), clouds_p3D_xyzs_[j][k+1](2));
                viewer->addLine(p1, p2, 0.0, 0.0, 0.99, id_line, 0);
              }
            }
          }
        }
      }

    } else {
      viewer->setPointCloudRenderingProperties(
          pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, id);
      viewer->setPointCloudRenderingProperties(
          pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, id);
    }
  }
}

void PCLViewer::render_clouds_p3Dd(
    pcl::visualization::PCLVisualizer::Ptr viewer) {
  for (unsigned int i = 0; i < clouds_p3Dd_.size(); i++) {
    if (!clouds_p3Dd_[i] || 
        (i == 0 && (!display_points3Dd_1_ || display_agg_)) ||
        (i == 1 && (!display_points3Dd_2_ || display_agg_)))
      continue;

    std::string id = "p3Dd" + std::to_string(i);

    pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud = clouds_p3Dd_[i];
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZ> rgb(
        cloud);
    viewer->addPointCloud<pcl::PointXYZ>(cloud, id);
    if (i == 0 || i == 2) {
      viewer->setPointCloudRenderingProperties(
          pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, id);
      viewer->setPointCloudRenderingProperties(
          pcl::visualization::PCL_VISUALIZER_COLOR, 0.25, 0.53, 0.56, id);

      if (i==2) { // Add spheres with the same color as the points
        for (unsigned int j = 0; j < cloud->size(); j++) {
          std::string id_sphere_avg = "sphere_" + std::to_string(i) + "_" + std::to_string(j) + "_avg";
          std::string id_sphere_max = "sphere_" + std::to_string(i) + "_" + std::to_string(j) + "_max";
          
          // Add sphere with radious clouds_p3D_range_[j].first
          viewer->addSphere(cloud->points[j], clouds_p3D_range_[j].first, 0.25, 0.53, 0.56, id_sphere_avg, 0);
          viewer->addSphere(cloud->points[j], clouds_p3D_range_[j].second, 0.25, 0.53, 0.56, id_sphere_max, 0);

          // Set transparency: avg = 0.3, max = 0.1
          viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, sphere_transparency_avg_, id_sphere_avg);
          viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, sphere_transparency_max_, id_sphere_max);
        }
      }
      
    } else {
      viewer->setPointCloudRenderingProperties(
          pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, id);
      viewer->setPointCloudRenderingProperties(
          pcl::visualization::PCL_VISUALIZER_COLOR, 0.81, 0.23, 0.32, id);
    }
  }
}

void PCLViewer::render_meshes(pcl::visualization::PCLVisualizer::Ptr viewer) {
  for (unsigned int i = 0; i < clouds_mesh_.size(); i++) {
    if (!clouds_mesh_[i]) continue;

    std::string id = "tri" + std::to_string(i);

    if (display_mesh_pts_) {
      pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud = clouds_mesh_[i];
      pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZ> rgb(
          cloud);
      viewer->addPointCloud<pcl::PointXYZ>(cloud, id);
      if (i == 0) {
        viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, id);
        viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, id);
      } else {
        viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, id);
        viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 0.14, id);
      }
    }

    if (display_mesh_tri_) {
      for (unsigned int j = 0; j < meshes_[i].size(); j++) {
        std::string id_poly =
            "tri_" + std::to_string(i) + "_" + std::to_string(j);
        pcl::PointCloud<pcl::PointXYZ>::ConstPtr polygon = meshes_[i][j]->poly_;
        if (i == 0 || i == 2)
          viewer->addPolygon<pcl::PointXYZ>(polygon, 0.0, 1.0, 0.0, id_poly, 0);
        else
          viewer->addPolygon<pcl::PointXYZ>(polygon, 0.1, 0.1, 1.0, id_poly, 0);
      }
    }
  }
}

void PCLViewer::DisplayModeMesh() {
  display_meshes_ = true;
  display_analyze_ = false;
  display_compare_ = false;
  display_agg_ = false;
}

void PCLViewer::DisplayModeAnalyze() {
  display_meshes_ = false;
  display_analyze_ = true;
  display_compare_ = false;
  display_agg_ = false;
}

void PCLViewer::DisplayModeCompare() {
  display_meshes_ = false;
  display_analyze_ = false;
  display_compare_ = true;
  display_agg_ = false;
}

void PCLViewer::DisplayModeAgg() {
  display_meshes_ = false;
  display_analyze_ = false;
  display_compare_ = false;
  display_agg_ = true;
}
