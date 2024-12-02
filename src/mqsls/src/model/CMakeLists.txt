

file(GLOB SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cpp
)

add_library(mqsls_model
    ${SRC}
)
target_compile_options(mqsls_model PRIVATE -Wno-unused-function)