
set(MAVLINK_GIT_DIR "${CMAKE_CURRENT_LIST_DIR}")
set(MAVLINK_LIBRARY_DIR "${CMAKE_BINARY_DIR}/mavlink")
set(MAVLINK_DIALECT "gcs_transfer")

# generate the dialect
add_custom_command(
    OUTPUT ${MAVLINK_LIBRARY_DIR}/${MAVLINK_DIALECT}/${MAVLINK_DIALECT}.h
    COMMAND
    ${PYTHON_EXECUTABLE} ${MAVLINK_GIT_DIR}/pymavlink/tools/mavgen.py
    --lang C --wire-protocol 2.0
        #--no-validate
        #--strict-units
    --output ${MAVLINK_LIBRARY_DIR}
    ${CMAKE_CURRENT_LIST_DIR}/message_definitions/${MAVLINK_DIALECT}.xml > ${CMAKE_CURRENT_BINARY_DIR}/mavgen_${MAVLINK_DIALECT}.log
    DEPENDS
    ${MAVLINK_GIT_DIR}/pymavlink/tools/mavgen.py
    ${CMAKE_CURRENT_LIST_DIR}/message_definitions/${MAVLINK_DIALECT}.xml


    COMMENT "Generating Mavlink ${MAVLINK_DIALECT}: ${CMAKE_CURRENT_LIST_DIR}/message_definitions/${MAVLINK_DIALECT}.xml"
)
add_custom_target(mavlink_c_generate DEPENDS ${MAVLINK_LIBRARY_DIR}/${MAVLINK_DIALECT}/${MAVLINK_DIALECT}.h)
set_source_files_properties(${MAVLINK_LIBRARY_DIR}/${MAVLINK_DIALECT}/${MAVLINK_DIALECT}.h PROPERTIES GENERATED true)

# mavlink header only library
add_library(mavlink_c INTERFACE)
target_compile_options(mavlink_c INTERFACE -Wno-address-of-packed-member -Wno-cast-align)
target_sources(mavlink_c
    INTERFACE
    ${MAVLINK_LIBRARY_DIR}/${MAVLINK_DIALECT}/${MAVLINK_DIALECT}.h
)
target_include_directories(mavlink_c
    INTERFACE
    ${MAVLINK_LIBRARY_DIR}
    ${MAVLINK_LIBRARY_DIR}/${MAVLINK_DIALECT}
)

add_executable(gcs_transfer_node
  src/gcs_transfer.cpp
)
target_link_libraries(gcs_transfer_node mavlink_c Boost::system "${cpp_typesupport_target}")
ament_target_dependencies(gcs_transfer_node rclcpp)
target_include_directories(gcs_transfer_node
    INTERFACE
    ${MAVLINK_LIBRARY_DIR}
    ${MAVLINK_LIBRARY_DIR}/${MAVLINK_DIALECT}
    ${Boost_INCLUDE_DIR}
)
target_compile_options(gcs_transfer_node PRIVATE -Wno-pedantic)
install(TARGETS
  gcs_transfer_node
  DESTINATION lib/${PROJECT_NAME}
)