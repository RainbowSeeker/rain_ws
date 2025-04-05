#include <rclcpp/rclcpp.hpp>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/visualization/cloud_viewer.h>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <iostream>
#include <thread>
#include <string>
#include <filesystem>

class PCLDetector : public rclcpp::Node {
public:
    PCLDetector() : Node("pcl_detector") {
        _cloud_sub = this->create_subscription<sensor_msgs::msg::PointCloud2>(
            "/livox/lidar", 10, std::bind(&PCLDetector::cloud_callback, this, std::placeholders::_1));

        _cloud_viewer->runOnVisualizationThread([this](pcl::visualization::PCLVisualizer& viewer) {

        });

        RCLCPP_INFO(this->get_logger(), "PCLDetector node started, waiting for point clouds...");
    }

    ~PCLDetector() {
    }

private:
    void cloud_callback(const sensor_msgs::msg::PointCloud2::SharedPtr msg) {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::fromROSMsg(*msg, *cloud);

        // visualizer point cloud
        _cloud_viewer->showCloud(cloud);

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
    pcl::visualization::CloudViewer::Ptr _cloud_viewer {new pcl::visualization::CloudViewer("PCL Viewer")};

    // detail
};


int main(int argc, char** argv) 
{
    rclcpp::init(argc, argv);

    rclcpp::spin(std::make_shared<PCLDetector>());

    rclcpp::shutdown();
}