# Set the project name
project(tools)

# find inmost package
find_package(inmost 0.3.0 REQUIRED)

# include inmost variables
link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})


# Add a library with the above sources
add_library(tools STATIC 
    src/stringtrim.cpp
    src/timer.cpp
    src/vecmath.cpp
    )

# link with external libraries
target_link_libraries(tools PUBLIC ${INMOST_LIBRARIES})

# link with internal libraries
target_link_libraries(tools PRIVATE sub::variables)

# make internal alias
add_library(sub::tools ALIAS tools)

# include headers
target_include_directories(tools
    PUBLIC ${PROJECT_SOURCE_DIR}/inc
)