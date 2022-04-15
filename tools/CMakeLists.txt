# Set the project name
project(tools)

# Add a library with the above sources
add_library(tools STATIC src/stringtrim.cpp)
add_library(sub::tools ALIAS tools)

target_include_directories(tools
    PUBLIC ${PROJECT_SOURCE_DIR}/inc
)