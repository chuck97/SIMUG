project(test_vec_transition)

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
    simug
)

# link with external libraries
target_link_libraries(${PROJECT_NAME} ${INMOST_LIBRARIES})

# to run tests
add_test(${PROJECT_NAME} ${PROJECT_NAME})
add_dependencies(check ${PROJECT_NAME})
