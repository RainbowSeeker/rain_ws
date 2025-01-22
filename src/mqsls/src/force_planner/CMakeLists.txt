add_subdirectory(model)

add_executable(force_planner_node
    force_planner.cpp
)
target_link_libraries(force_planner_node
    model_force_opt "${cpp_typesupport_target}"
)
ament_target_dependencies(force_planner_node rclcpp Eigen3)
install(TARGETS force_planner_node DESTINATION lib/${PROJECT_NAME})
