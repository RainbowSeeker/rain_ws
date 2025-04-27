#include <pcl/io/pcd_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/features/moment_of_inertia_estimation.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/cvfh.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <chrono>

typedef pcl::PointXYZ PointT;

struct ObjectTemplate {
    pcl::PointCloud<PointT>::Ptr cloud;  // 模板点云
    pcl::VFHSignature308 signature;     // 模板特征
    std::string label;                  // 类别标签
};

pcl::VFHSignature308 extract_cvfh_feature(pcl::PointCloud<PointT>::Ptr cloud) {
    // 计算法线
    pcl::NormalEstimation<PointT, pcl::Normal> ne;
    ne.setInputCloud(cloud);
    pcl::search::KdTree<PointT>::Ptr tree(new pcl::search::KdTree<PointT>());
    ne.setSearchMethod(tree);
    pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
    ne.setRadiusSearch(0.3);
    ne.compute(*normals);

    // CVFH特征提取
    pcl::CVFHEstimation<PointT, pcl::Normal, pcl::VFHSignature308> cvfh;
    cvfh.setInputCloud(cloud);
    cvfh.setInputNormals(normals);
    cvfh.setSearchMethod(tree);
    
    pcl::PointCloud<pcl::VFHSignature308>::Ptr features(new pcl::PointCloud<pcl::VFHSignature308>);
    cvfh.compute(*features);
    if (features->size() > 0) {
        return features->at(0);
    } else {
        throw std::runtime_error("No features extracted");
    }
}

bool flann_template_matching(pcl::PointCloud<PointT>::Ptr cluster,
    const std::vector<ObjectTemplate>& templates,
    pcl::KdTreeFLANN<pcl::VFHSignature308>& feature_tree) 
{
    // 提取查询特征
    pcl::VFHSignature308 query_feature = extract_cvfh_feature(cluster);

    // K近邻搜索
    int K = 1; // 找最相似的一个模板
    std::vector<int> indices(K);
    std::vector<float> distances(K);

    if(feature_tree.nearestKSearch(query_feature, K, indices, distances) > 0 
        && distances[0] < 1000) {
        
        std::cout << "Matched template: " << templates[indices[0]].label
                    << " | Distance: " << distances[0] << std::endl;

        return true;
    }

    std::cout << "No match found." << std::endl;
    return false;
}

int main() {
    // 1. 读取点云
    pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
    pcl::io::loadPCDFile("/home/yangyu/rain_ws/install/pcds/cloud_test.pcd", *cloud);

    std::cout << "PointCloud before filtering has: " << cloud->size () << " data points." << std::endl; //*

    // 1.5 读取模板
    std::vector<ObjectTemplate> templates;

    ObjectTemplate uav;
    uav.label = "uav";
    uav.cloud.reset(new pcl::PointCloud<PointT>);
    pcl::io::loadPCDFile("/home/yangyu/rain_ws/install/pcds/cluster_uav.pcd", *uav.cloud);
    uav.signature = extract_cvfh_feature(uav.cloud);
    templates.push_back(uav);

    // 构建特征KD树
    pcl::PointCloud<pcl::VFHSignature308>::Ptr template_features(new pcl::PointCloud<pcl::VFHSignature308>);
    for(const auto& temp : templates) {
        template_features->push_back(temp.signature);
    }

    pcl::KdTreeFLANN<pcl::VFHSignature308> feature_tree;
    feature_tree.setInputCloud(template_features);

    // start to timer
    auto start_time = std::chrono::high_resolution_clock::now();

    // 2. 预处理
    // 降采样
    pcl::VoxelGrid<PointT> vg;
    vg.setInputCloud(cloud);
    vg.setLeafSize(0.01f, 0.01f, 0.01f); // 根据数据调整
    pcl::PointCloud<PointT>::Ptr filtered_cloud(new pcl::PointCloud<PointT>);
    vg.filter(*filtered_cloud);

    // 去噪
    pcl::StatisticalOutlierRemoval<PointT> sor;
    sor.setInputCloud(filtered_cloud);
    sor.setMeanK(50);
    sor.setStddevMulThresh(1.0);
    sor.filter(*filtered_cloud);

    std::cout << "PointCloud after filtering has: " << filtered_cloud->size ()  << " data points." << std::endl; //*


    // 3. 地面分割（RANSAC）
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
    pcl::SACSegmentation<PointT> seg;
    seg.setOptimizeCoefficients(true);
    seg.setModelType(pcl::SACMODEL_PLANE);
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setMaxIterations (100);
    seg.setDistanceThreshold(0.05); // 地面点阈值

    // 多次提取
    int nr_points = (int)filtered_cloud->size();
    while (filtered_cloud->size() > 0.3 * nr_points) {
        // 设置输入点云
        seg.setInputCloud(filtered_cloud);
        // 执行分割，得到内点和模型系数
        seg.segment(*inliers, *coefficients);
        if (inliers->indices.size() == 0) {
            std::cout << "Could not estimate a planar model for the given dataset." << std::endl;
            break;
        }
        // 提取平面点
        pcl::ExtractIndices<PointT> extract;
        extract.setInputCloud(filtered_cloud);
        extract.setIndices(inliers);
        extract.setNegative(false); // 提取地面点
        pcl::PointCloud<PointT>::Ptr cloud_plane(new pcl::PointCloud<PointT>);
        extract.filter(*cloud_plane);
        std::cout << "PointCloud representing the planar component: " << cloud_plane->size() << " data points." << std::endl;
        
        // 提取非平面点
        extract.setNegative(true); // 提取非地面点
        pcl::PointCloud<PointT>::Ptr cloud_f(new pcl::PointCloud<PointT>);
        extract.filter(*cloud_f);
        *filtered_cloud = *cloud_f; // 更新输入点云
        std::cout << "PointCloud representing the non-planar component: " << filtered_cloud->size() << " data points." << std::endl;
    }

    // 4. 欧几里得聚类
    pcl::search::KdTree<PointT>::Ptr tree(new pcl::search::KdTree<PointT>);
    tree->setInputCloud(filtered_cloud);

    std::vector<pcl::PointIndices> cluster_indices;
    pcl::EuclideanClusterExtraction<PointT> ec;
    ec.setClusterTolerance(0.1); // 聚类距离阈值
    ec.setMinClusterSize(20);    // 最小点数
    ec.setMaxClusterSize(25000); // 最大点数
    ec.setSearchMethod(tree);
    ec.setInputCloud(filtered_cloud);
    ec.extract(cluster_indices);

    std::cout << "Number of clusters found: " << cluster_indices.size() << std::endl;

    // 5. 特征提取与分类
    pcl::PCDWriter writer;
    std::vector<pcl::PointCloud<PointT>::Ptr> clusters;
    pcl::visualization::PCLVisualizer viewer("Detection Result");
    viewer.addPointCloud<PointT>(filtered_cloud, "original");

    int cluster_id = 0;
    for (const auto& indices : cluster_indices) {
        pcl::PointCloud<PointT>::Ptr cluster(new pcl::PointCloud<PointT>);
        for (const auto& idx : indices.indices) {
            cluster->push_back((*filtered_cloud)[idx]);
        }
        clusters.push_back(cluster);

        // 计算包围盒
        pcl::MomentOfInertiaEstimation<PointT> feature_extractor;
        feature_extractor.setInputCloud(cluster);
        feature_extractor.compute();

        PointT min_point, max_point, center;
        Eigen::Matrix3f rotation;
        feature_extractor.getOBB(min_point, max_point, center, rotation);

        // 过滤条件：根据尺寸过滤
        float length = max_point.x - min_point.x;
        float width = max_point.y - min_point.y;
        float height = max_point.z - min_point.z;

        if (length > 0.05 && length < 1.5 && 
            width > 0.05 && width < 1.0 &&
            height > 0.05 && height < 1.0) 
        {
            Eigen::Vector3f position(center.x, center.y, center.z);
            Eigen::Quaternionf quat(rotation);
            // 添加包围盒到可视化
            viewer.addCube(position, quat,
                          max_point.x - min_point.x,
                          max_point.y - min_point.y,
                          max_point.z - min_point.z,
                          "cluster_" + std::to_string(cluster_id));
            viewer.setShapeRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_REPRESENTATION, 
                pcl::visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME,
                "cluster_" + std::to_string(cluster_id)
            );

            std::cout << "Cluster " << cluster_id << ": size = " << cluster->size() 
                      << ", length = " << length 
                      << ", width = " << width 
                      << ", height = " << height 
                      << std::endl;

            writer.write<PointT>(
                "/home/yangyu/rain_ws/install/pcds/cluster_" + std::to_string(cluster_id) + ".pcd", 
                *cluster, 
                false);
            
            // 对每个聚类执行快速匹配
            if (flann_template_matching(cluster, templates, feature_tree)) {
                // set color
                viewer.setShapeRenderingProperties(
                    pcl::visualization::PCL_VISUALIZER_COLOR, 
                    1.0, 0.0, 0.0,
                    "cluster_" + std::to_string(cluster_id));
            }
        }
        cluster_id++;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Processing time: " << duration.count() << " ms" << std::endl;

    // 6. 可视化
    while (!viewer.wasStopped()) {
        viewer.spinOnce();
    }

    return 0;
}