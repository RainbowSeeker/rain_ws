
find_package(control_toolbox REQUIRED)

add_executable(mc_formation_control_node
    mc_formation_control.cpp
)
target_link_libraries(mc_formation_control_node "${cpp_typesupport_target}")
ament_target_dependencies(mc_formation_control_node 
        rclcpp px4_msgs px4_ros_com control_toolbox Eigen3)
install(TARGETS mc_formation_control_node DESTINATION lib/${PROJECT_NAME})