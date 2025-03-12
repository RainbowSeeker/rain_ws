
add_executable(rtk_plugin_node
    rtk_plugin.cpp
)
target_link_libraries(rtk_plugin_node "${cpp_typesupport_target}")
ament_target_dependencies(rtk_plugin_node rclcpp)
install(TARGETS rtk_plugin_node DESTINATION lib/${PROJECT_NAME})