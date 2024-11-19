
add_executable(serial_udp_bridge
    serial_udp_bridge.cpp
)
target_link_libraries(serial_udp_bridge Boost::system)
install(TARGETS serial_udp_bridge DESTINATION lib/${PROJECT_NAME}) 