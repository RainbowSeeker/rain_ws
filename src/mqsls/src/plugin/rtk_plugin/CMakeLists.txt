aux_source_directory(. RTK_SRC)

add_executable(rtk_plugin_node
    ${RTK_SRC}
)

option(RTK_TESTING "Enable testing for RTK plugin" OFF)
if(RTK_TESTING)
    target_compile_options(rtk_plugin_node PRIVATE -g -O0)
    target_link_options(rtk_plugin_node PRIVATE -g)
endif()

target_link_libraries(rtk_plugin_node "${cpp_typesupport_target}")
ament_target_dependencies(rtk_plugin_node rclcpp)
install(TARGETS rtk_plugin_node DESTINATION lib/${PROJECT_NAME})