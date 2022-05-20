# Set the project name
project(forcing)

# find inmost package
find_package(inmost 0.3.0 REQUIRED)

# include inmost variables
link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

# Make library from sources
add_library(forcing STATIC 
  src/forcing.cpp
  )

# link with external libraries
target_link_libraries(forcing PUBLIC ${INMOST_LIBRARIES})

# link with internal libraries
target_link_libraries(forcing PRIVATE 
  sub::variables
  sub::tools
  sub::mesh
)

# make internal alias
add_library(sub::forcing ALIAS forcing)

# include headers
target_include_directories(forcing
    PUBLIC ${PROJECT_SOURCE_DIR}/inc
)