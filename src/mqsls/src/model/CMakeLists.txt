

file(GLOB SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.c*
)

add_library(mqsls_model STATIC
    ${SRC}
)
target_compile_options(mqsls_model PRIVATE -Wno-unused-function)