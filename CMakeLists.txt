cmake_minimum_required (VERSION 3.17.0)

project(SIMUG)

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
option(USE_OMP "Compile with OpenMP support" OFF)
option(USE_PARTITIONER "Compile with partitioner" OFF)
option(USE_SOLVER "Compile with solver" ON)
option(USE_FP64 "Use 64-bit floating-point types (usually on)" ON)
option(COMPILE_TESTS "Compile tests" OFF)
option(COMPILE_EXAMPLES "Compile examples" ON)
option(INTEL_COMPILER "Intel compiler" OFF)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  MESSAGE(STATUS "Intel compiler")
  add_definitions(-DINTEL_COMPILER)
endif()

if (INTEL_COMPILER)
  link_libraries(stdc++fs)
endif()

if (USE_MPI)
  MESSAGE(STATUS "MPI ENABLED")
  add_definitions(-DUSE_MPI)
endif()

if (USE_OMP)
  MESSAGE(STATUS "OpenMP ENABLED")
  add_definitions(-DUSE_OMP)
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
  MESSAGE(STATUS "Test compiling ENABLED")
endif()


# include headers directories
include_directories("source/headers")
include_directories("source/headers/variables")

# Add subdirectories
add_subdirectory(source)

# make library
add_library(simug STATIC ${SOURCE} ${HEADER})

# find inmost package
find_package(inmost 0.3.0 REQUIRED)

# include inmost variables
link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

# link library with inmost
# link with external libraries
target_link_libraries(simug PUBLIC ${INMOST_LIBRARIES})

if(COMPILE_TESTS)
  enable_testing()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}) 
  add_subdirectory(tests)
endif(COMPILE_TESTS)

if (COMPILE_EXAMPLES)
  MESSAGE(STATUS "Examples compiling ENABLED")
  add_subdirectory(examples)
endif()