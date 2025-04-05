cmake_minimum_required(VERSION 3.8)
project(pcl_det)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PCL_CUSTOM_VERSION "1.15")
if(${PCL_CUSTOM_VERSION} VERSION_EQUAL 1.15)
  set(PCL_DIR "/usr/local/share/pcl-1.15")
elseif(${PCL_CUSTOM_VERSION} VERSION_EQUAL 1.12)
  set(PCL_DIR "/usr/lib/x86_64-linux-gnu/cmake/pcl")
else()
  message(FATAL_ERROR "PCL version not supported")
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(Boost REQUIRED COMPONENTS system date_time)
find_package(PCL ${PCL_CUSTOM_VERSION} REQUIRED)
find_package(pcl_ros REQUIRED)


message(STATUS "PCL version: ${PCL_VERSION}")

include_directories(
  include
  ${PCL_INCLUDE_DIRS}
)

add_subdirectory(src)

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
