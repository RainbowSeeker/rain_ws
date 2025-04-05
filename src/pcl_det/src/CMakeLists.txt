
add_executable(pcl_det_node 
   pcl_det.cpp
)
ament_target_dependencies(pcl_det_node rclcpp std_msgs pcl_ros)
target_link_libraries(pcl_det_node ${PCL_LIBRARIES})
install(TARGETS pcl_det_node DESTINATION lib/${PROJECT_NAME})
