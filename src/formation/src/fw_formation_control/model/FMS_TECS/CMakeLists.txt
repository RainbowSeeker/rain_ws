
file(GLOB_RECURSE FMS_SRC
    src/*.*
)

add_library(FMS
    ${FMS_SRC}
    FMS.h
)
