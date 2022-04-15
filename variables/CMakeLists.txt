project(variables)

add_library(variables INTERFACE)

add_library(sub::variables ALIAS variables)

target_include_directories(variables
    INTERFACE ${PROJECT_SOURCE_DIR}/inc
)