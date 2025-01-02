
# controller
aux_source_directory(controller CTL_SRC)

add_library(model_controller STATIC
    ${CTL_SRC}
)
target_compile_options(model_controller PRIVATE -Wno-unused-function)

# force optimization
aux_source_directory(force_opt OPT_SRC)

add_library(model_force_opt STATIC
    ${OPT_SRC}
)