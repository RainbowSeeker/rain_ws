cmake_minimum_required(VERSION 3.8)
project(mqsls)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# cmake flag
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(
  -Wno-unused-function
  -Wno-unused-parameter
  -Wno-unused-variable
)

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rosidl_default_generators REQUIRED)

# ##############################################################################
# Generate ROS messages, ROS2 interfaces and IDL files #
# ##############################################################################

# get all msg files
set(MSGS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/msg")
file(GLOB ALL_MSGS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${MSGS_DIR}/*.*")

# Generate introspection typesupport for C and C++ and IDL files
rosidl_generate_interfaces(${PROJECT_NAME} 
  ${ALL_MSGS}
)

ament_export_dependencies(rosidl_default_runtime)
rosidl_get_typesupport_target(cpp_typesupport_target ${PROJECT_NAME} rosidl_typesupport_cpp)

# add include for library
# install(DIRECTORY include/ DESTINATION include/${PROJECT_NAME})
# ament_export_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# ##############################################################################
# Build Packages #
# ##############################################################################

find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(px4_msgs REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(formation REQUIRED)
find_package(control_toolbox REQUIRED)

include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
  ${formation_INCLUDE_DIRS}
)

# src
add_subdirectory(src)

# Install config files.
install(DIRECTORY config DESTINATION share/${PROJECT_NAME})

# Install gz files.
install(DIRECTORY gz DESTINATION share/${PROJECT_NAME})

# Install launch files.
install(DIRECTORY launch DESTINATION share/${PROJECT_NAME})

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
