#include <rclcpp/rclcpp.hpp>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/features/moment_of_inertia_estimation.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/cvfh.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <iostream>
#include <thread>
#include <string>
#include <filesystem>

using PointT = pcl::PointXYZI;

class PCLDetector : public rclcpp::Node {
public:
    PCLDetector() : Node("pcl_detector") {
        _cloud_sub = this->create_subscription<sensor_msgs::msg::PointCloud2>(
            "/livox/lidar", 10, std::bind(&PCLDetector::cloud_callback, this, std::placeholders::_1));

        // _viewer = std::make_shared<pcl::visualization::PCLVisualizer>("3D Viewer");
        // _viewer->setBackgroundColor(0, 0, 0);
        // _viewer->addCoordinateSystem(1.0);
        // _viewer->initCameraParameters();

        RCLCPP_INFO(this->get_logger(), "PCLDetector node started, waiting for point clouds...");
    }

    ~PCLDetector() {
    }

private:
    void cloud_callback(const sensor_msgs::msg::PointCloud2::SharedPtr msg) {
        pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>), filtered_cloud(new pcl::PointCloud<PointT>);
        pcl::fromROSMsg(*msg, *cloud);

        // exclude points with intensity < 150
        pcl::ConditionAnd<PointT>::Ptr condition(new pcl::ConditionAnd<PointT>());
        condition->addComparison(pcl::FieldComparison<PointT>::Ptr(
            new pcl::FieldComparison<PointT>("intensity", pcl::ComparisonOps::LT, 50)
        ));
        
        pcl::ConditionalRemoval<PointT> condrem;
        condrem.setCondition(condition);
        condrem.setInputCloud(cloud);
        condrem.setKeepOrganized(false);
        condrem.filter(*cloud);

        _cloud_viewer->showCloud(cloud);
        return;

        pcl::search::KdTree<PointT>::Ptr tree(new pcl::search::KdTree<PointT>);
        tree->setInputCloud(filtered_cloud);

        std::vector<pcl::PointIndices> cluster_indices;
        pcl::EuclideanClusterExtraction<PointT> ec;
        ec.setClusterTolerance(0.1); // 聚类距离阈值
        ec.setMinClusterSize(10);    // 最小点数
        ec.setMaxClusterSize(25000); // 最大点数
        ec.setSearchMethod(tree);
        ec.setInputCloud(filtered_cloud);
        ec.extract(cluster_indices);

        std::cout << "Number of clusters found: " << cluster_indices.size() << std::endl;

        // 5. 特征提取与分类
        int cluster_id = 0;
        for (const auto& indices : cluster_indices) {
            pcl::PointCloud<PointT>::Ptr cluster(new pcl::PointCloud<PointT>);
            for (const auto& idx : indices.indices) {
                cluster->push_back((*filtered_cloud)[idx]);
            }

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

            // if (length > 0.05 && length < 1.5 && 
            //     width > 0.05 && width < 1.0 &&
            //     height > 0.05 && height < 1.0) 
            {
                Eigen::Vector3f position(center.x, center.y, center.z);
                Eigen::Quaternionf quat(rotation);
                // 添加包围盒到可视化
                // _viewer->addCube(position, quat,
                //             max_point.x - min_point.x,
                //             max_point.y - min_point.y,
                //             max_point.z - min_point.z,
                //             "cluster_" + std::to_string(cluster_id));
                // _viewer->setShapeRenderingProperties(
                //     pcl::visualization::PCL_VISUALIZER_REPRESENTATION, 
                //     pcl::visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME,
                //     "cluster_" + std::to_string(cluster_id)
                // );

                std::cout << "Cluster " << cluster_id << ": size = " << cluster->size() 
                        << ", length = " << length 
                        << ", width = " << width 
                        << ", height = " << height 
                        << std::endl;

                // writer.write<PointT>(
                //     "/home/yangyu/rain_ws/install/pcds/cluster_" + std::to_string(cluster_id) + ".pcd", 
                //     *cluster, 
                //     false);
            }
            cluster_id++;
        }

        // visualizer point cloud
        if (!_viewer->updatePointCloud<PointT>(filtered_cloud, "cloud")) {
            pcl::visualization::PointCloudColorHandlerGenericField<PointT> 
                intensity_color(filtered_cloud, "intensity");
            _viewer->addPointCloud<PointT>(filtered_cloud, intensity_color, "cloud");
            _viewer->setPointCloudRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "cloud");
        }
        _viewer->spinOnce();
        _viewer->removeAllShapes();

        // _writer.write<PointT>(
        //         "/home/yangyu/rain_ws/install/pcds/cloud_" + std::to_string(1) + ".pcd", 
        //         *cloud,
        //         false);

        // exit(0);
        // Save the point cloud to a PCD file
        // std::string filename = _save_path + "/cloud_" + std::to_string(msg->header.stamp.sec) + ".pcd";
        // if (pcl::io::savePCDFileASCII(filename, *cloud) == 0) {
        //     RCLCPP_INFO(this->get_logger(), "Saved point cloud to %s", filename.c_str());
        // } else {
        //     RCLCPP_ERROR(this->get_logger(), "Failed to save point cloud to %s", filename.c_str());
        // }        
    }

    // Subscriber
    rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr  _cloud_sub;

    // PCL
    pcl::visualization::PCLVisualizer::Ptr _viewer;
    pcl::visualization::CloudViewer::Ptr _cloud_viewer { new pcl::visualization::CloudViewer("PCL Viewer") };

    // detail
    pcl::PCDWriter _writer;
};


int main(int argc, char** argv) 
{
    rclcpp::init(argc, argv);

    rclcpp::spin(std::make_shared<PCLDetector>());

    rclcpp::shutdown();
}