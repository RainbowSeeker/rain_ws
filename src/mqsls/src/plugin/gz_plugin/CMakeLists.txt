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

add_executable(gz_plugin_node gz_plugin.cpp)
target_link_libraries(gz_plugin_node ${GZ_TRANSPORT_LIB} "${cpp_typesupport_target}")
ament_target_dependencies(gz_plugin_node rclcpp Eigen3)
install(TARGETS gz_plugin_node DESTINATION lib/${PROJECT_NAME})