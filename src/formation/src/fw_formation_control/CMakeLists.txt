
add_subdirectory(model)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/model
)
set(SRC 
    fw_formation_control.cpp
)

add_executable(fw_formation_control_node ${SRC})
target_link_libraries(fw_formation_control_node FMS Eigen3::Eigen "${cpp_typesupport_target}")
ament_target_dependencies(fw_formation_control_node rclcpp px4_msgs px4_ros_com)
install(TARGETS fw_formation_control_node DESTINATION lib/${PROJECT_NAME})

# add_executable(cross_node cross_node.cpp)
# ament_target_dependencies(cross_node rclcpp px4_msgs px4_ros_com)
# install(TARGETS cross_node DESTINATION lib/${PROJECT_NAME})

# set(TEST false)
# if (TEST)
#   add_executable(test_node src/test.cpp)
#   ament_target_dependencies(test_node rclcpp px4_msgs px4_ros_com)
#   install(TARGETS test_node DESTINATION lib/${PROJECT_NAME})
# endif ()
