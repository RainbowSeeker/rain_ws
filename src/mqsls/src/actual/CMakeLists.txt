
add_executable(actual_mqsls_control_node
    actual_mqsls_control.cpp
)
target_link_libraries(actual_mqsls_control_node model_controller "${cpp_typesupport_target}")
ament_target_dependencies(actual_mqsls_control_node 
        rclcpp px4_msgs 
        control_toolbox Eigen3)
install(TARGETS actual_mqsls_control_node DESTINATION lib/${PROJECT_NAME})