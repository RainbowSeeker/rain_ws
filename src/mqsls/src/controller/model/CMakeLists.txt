aux_source_directory(controller SRC)

add_library(model_controller STATIC
    ${SRC}
)
target_compile_options(model_controller PRIVATE -Wno-unused-function)