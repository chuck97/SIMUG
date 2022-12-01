project(test_advection)

# find inmost package
find_package(inmost 0.3.0 REQUIRED)

# create the executable
add_executable(${PROJECT_NAME} 
    main.cpp 
    )

# link with internal libraries
target_link_libraries(${PROJECT_NAME}
    simug
)

# include headers
target_include_directories(${PROJECT_NAME}
    PRIVATE ${PROJECT_SOURCE_DIR}/inc
)

# to run tests
add_test(${PROJECT_NAME} ${PROJECT_NAME})
add_dependencies(check ${PROJECT_NAME})