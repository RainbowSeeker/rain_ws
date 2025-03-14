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

# Compile Plugin
# gz_find_package(gz-plugin2 REQUIRED COMPONENTS register)
# gz_find_package(gz-common5 REQUIRED)
# gz_find_package(gz-sim7 REQUIRED)
# set(GZ_PLUGIN_VER ${gz-plugin2_VERSION_MAJOR})
# set(GZ_COMMON_VER ${gz-common5_VERSION_MAJOR})
# set(GZ_SIM_VER ${gz-sim7_VERSION_MAJOR})
 
# # Add sources for each plugin to be registered.
# add_library(auto_run_sim_plugin SHARED auto_run_sim_plugin.cpp)
# target_link_libraries(auto_run_sim_plugin
#   gz-common${GZ_COMMON_VER}::gz-common${GZ_COMMON_VER}
#   ${GZ_TRANSPORT_LIB}
#   gz-plugin${GZ_PLUGIN_VER}::gz-plugin${GZ_PLUGIN_VER}
#   gz-sim${GZ_SIM_VER}::gz-sim${GZ_SIM_VER}
# )
# install(TARGETS auto_run_sim_plugin DESTINATION share/${PROJECT_NAME}/gz/plugins/)