project(cgrid_plane_advection)

# find inmost package
find_package(inmost 0.3.0 REQUIRED)

# include inmost variables
link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

# create the executable
add_executable(${PROJECT_NAME} 
    main.cpp 
    )

# link with internal libraries
target_link_libraries(${PROJECT_NAME}
    sub::variables
    sub::mesh
    sub::tools
    sub::forcing
    sub::advection
)

# link with external libraries
target_link_libraries(${PROJECT_NAME} ${INMOST_LIBRARIES})