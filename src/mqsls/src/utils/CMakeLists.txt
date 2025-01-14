
add_executable(ros2_mocap_bridge
    ros2_mocap_bridge.cpp
)
target_link_libraries(ros2_mocap_bridge "${cpp_typesupport_target}")
ament_target_dependencies(ros2_mocap_bridge 
        rclcpp px4_msgs Eigen3)
install(TARGETS ros2_mocap_bridge DESTINATION lib/${PROJECT_NAME})