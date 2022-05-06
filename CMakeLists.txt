cmake_minimum_required (VERSION 3.17.0)

project(SIMUG_v1)

# C/C++ compilers
enable_language (CXX C)
#set(CMAKE_CXX_COMPILER icpc)
#set(CMAKE_C_COMPILER icc)

# default build type is debug
set(CMAKE_BUILD_TYPE Debug)

# set C++ standart
if(NOT "${CMAKE_CXX_STANDARD}")
  set(CMAKE_CXX_STANDARD 20)
endif()

# options
option(USE_MPI "Compile with MPI support" OFF)
option(USE_PARTITIONER "Compile with partitioner" OFF)
option(USE_SOLVER "Compile with solver" OFF)
option(USE_FP64 "Use 64-bit floating-point types (usually on)" ON)
option(COMPILE_TESTS "Compile tests" ON)

if (USE_MPI)
  MESSAGE(STATUS "MPI enabled")
  add_definitions(-DUSE_MPI)
endif()

if (USE_PARTITIONER)
  MESSAGE(STATUS "Partitioner ENABLED")
  add_definitions(-DUSE_PARTITIONER)
endif()

if (USE_PARTITIONER_PARMETIS)
  MESSAGE(STATUS "Partitioner ParMetis ENABLED")
  add_definitions(-DUSE_PARTITIONER)
  add_definitions(-DUSE_PARTITIONER_PARMETIS)
endif()


if (USE_SOLVER)
  MESSAGE(STATUS "Linear system solver ENABLED")
  add_definitions(-DUSE_SOLVER)
endif()

if (USE_FP64)
  MESSAGE(STATUS "Real type = double")
  add_definitions(-DUSE_FP64)
endif()

if (COMPILE_TESTS)
  MESSAGE(STATUS "Test compiling enabled")
endif()



# Add subdirectories
add_subdirectory(variables)
add_subdirectory(tools)
#add_subdirectory(configuration)
add_subdirectory(mesh)


if(COMPILE_TESTS)
  enable_testing()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}) 
  add_subdirectory(tests)
endif(COMPILE_TESTS)

