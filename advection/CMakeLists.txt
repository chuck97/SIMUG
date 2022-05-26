# Set the project name
project(advection)

# find inmost package
find_package(inmost 0.3.0 REQUIRED)

# include inmost variables
link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

# Make library from sources
add_library(advection STATIC 
  src/cgrid_advection.cpp
  )

# link with external libraries
target_link_libraries(advection PUBLIC ${INMOST_LIBRARIES})

# link with internal libraries
target_link_libraries(advection PRIVATE 
  sub::variables
  sub::tools
  sub::mesh
)

# make internal alias
add_library(sub::advection ALIAS advection)

# include headers
target_include_directories(advection
    PUBLIC ${PROJECT_SOURCE_DIR}/inc
)