
add_subdirectory(model)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/model
)
set(SRC 
    fw_formation_control.cpp
)

add_executable(fw_formation_control_node ${SRC})
target_link_libraries(fw_formation_control_node FMS Eigen3::Eigen "${cpp_typesupport_target}")
ament_target_dependencies(fw_formation_control_node rclcpp px4_msgs)
install(TARGETS fw_formation_control_node DESTINATION lib/${PROJECT_NAME})