find_package(control_toolbox REQUIRED)

add_executable(mc_mqsls_control_node
  mc_mqsls_control.cpp
)

ament_target_dependencies(mc_mqsls_control_node 
        rclcpp px4_msgs px4_ros_com formation 
        control_toolbox Eigen3)
install(TARGETS mc_mqsls_control_node DESTINATION lib/${PROJECT_NAME})
