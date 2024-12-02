find_package(control_toolbox REQUIRED)
find_package(gz-transport NAMES gz-transport12)

if(gz-transport_FOUND)
  message(STATUS "Found Gazebo transport: ${gz-transport_VERSION}")
else()
  message(STATUS "Gazebo transport not found")
  return()
endif()

set(GZ_TRANSPORT_VER ${gz-transport_VERSION_MAJOR})

if(GZ_TRANSPORT_VER GREATER_EQUAL 12)
  set(GZ_TRANSPORT_LIB gz-transport${GZ_TRANSPORT_VER}::core)
else()
  set(GZ_TRANSPORT_LIB ignition-transport${GZ_TRANSPORT_VER}::core)
endif()

add_executable(ideal_mqsls_control_node
  ideal_mqsls_control.cpp
)
target_link_libraries(ideal_mqsls_control_node ${GZ_TRANSPORT_LIB} mqsls_model)

ament_target_dependencies(ideal_mqsls_control_node 
        rclcpp px4_msgs px4_ros_com formation 
        control_toolbox Eigen3)
install(TARGETS ideal_mqsls_control_node DESTINATION lib/${PROJECT_NAME})


# Compile Plugin
gz_find_package(gz-plugin2 REQUIRED COMPONENTS register)
gz_find_package(gz-common5 REQUIRED)
gz_find_package(gz-sim7 REQUIRED)
set(GZ_PLUGIN_VER ${gz-plugin2_VERSION_MAJOR})
set(GZ_COMMON_VER ${gz-common5_VERSION_MAJOR})
set(GZ_SIM_VER ${gz-sim7_VERSION_MAJOR})
 
# Add sources for each plugin to be registered.
add_library(virtual_wrench_plugin SHARED virtual_wrench_plugin.cpp)
target_link_libraries(virtual_wrench_plugin
  gz-common${GZ_COMMON_VER}::gz-common${GZ_COMMON_VER}
  ${GZ_TRANSPORT_LIB}
  gz-plugin${GZ_PLUGIN_VER}::gz-plugin${GZ_PLUGIN_VER}
  gz-sim${GZ_SIM_VER}::gz-sim${GZ_SIM_VER}
)
install(TARGETS virtual_wrench_plugin DESTINATION share/${PROJECT_NAME}/gz/plugins/)