add_subdirectory(model)

add_executable(3dof_mqsls_controller_node
  3dof_mqsls_controller.cpp
)
target_link_libraries(3dof_mqsls_controller_node model_controller "${cpp_typesupport_target}")
ament_target_dependencies(3dof_mqsls_controller_node rclcpp px4_msgs control_toolbox Eigen3)
install(TARGETS 3dof_mqsls_controller_node DESTINATION lib/${PROJECT_NAME})

add_executable(px4_actuator_node
  px4_actuator.cpp
)
target_link_libraries(px4_actuator_node "${cpp_typesupport_target}")
ament_target_dependencies(px4_actuator_node rclcpp px4_msgs Eigen3)
install(TARGETS px4_actuator_node DESTINATION lib/${PROJECT_NAME})