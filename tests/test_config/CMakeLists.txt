project(test_config)

# Create the executable
add_executable(${PROJECT_NAME} main.cpp)

# Link the static library from subprojects
target_link_libraries(${PROJECT_NAME}
    sub::variables
    sub::tools
    sub::configuration
)

# to run tests
add_test(${PROJECT_NAME} ${PROJECT_NAME})
add_dependencies(check ${PROJECT_NAME})
