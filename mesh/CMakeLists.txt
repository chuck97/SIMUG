# Set the project name
project(mesh)

# find inmost package
find_package(inmost 0.3.0 REQUIRED)

# include inmost variables
link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

# Make library from sources
add_library(mesh STATIC 
  src/data.cpp
  src/mesh.cpp
  src/gridinfo.cpp
  src/basis.cpp
  )

# link with external libraries
target_link_libraries(mesh PUBLIC ${INMOST_LIBRARIES})

# link with internal libraries
target_link_libraries(mesh PRIVATE 
  sub::variables
  sub::tools
)

# make internal alias
add_library(sub::mesh ALIAS mesh)

# include headers
target_include_directories(mesh
    PUBLIC ${PROJECT_SOURCE_DIR}/inc
)