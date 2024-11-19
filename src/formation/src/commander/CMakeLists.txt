
add_executable(commander_node
    commander.cpp
)

ament_target_dependencies(commander_node 
        rclcpp px4_msgs)
install(TARGETS commander_node DESTINATION lib/${PROJECT_NAME})