
add_executable(mocap_plugin_node
    mocap_plugin.cpp
)
target_link_libraries(mocap_plugin_node "${cpp_typesupport_target}")
ament_target_dependencies(mocap_plugin_node rclcpp px4_msgs Eigen3)
install(TARGETS mocap_plugin_node DESTINATION lib/${PROJECT_NAME})