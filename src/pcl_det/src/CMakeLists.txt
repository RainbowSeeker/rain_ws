
add_executable(pcl_det_node 
   pcl_det.cpp
)

option(PCL_TESTING "Enable testing for PCL" OFF)
if(PCL_TESTING)
    target_compile_options(pcl_det_node PRIVATE -g -O0)
    target_link_options(pcl_det_node PRIVATE -g)
endif()

ament_target_dependencies(pcl_det_node rclcpp std_msgs pcl_ros)
target_link_libraries(pcl_det_node ${PCL_LIBRARIES})
install(TARGETS pcl_det_node DESTINATION lib/${PROJECT_NAME})
