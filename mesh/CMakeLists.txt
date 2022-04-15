# Set the project name
project(mesh)

# libraries path
set(INMOST_SRC "/home/users/spetrov/INMOST_BUILD_PVTU/INMOST")
set(INMOST_BUILD "/home/users/spetrov/INMOST_BUILD_PVTU/INMOST-intel-build")
set(PETSC_DIR "/home/users/spetrov/INMOST_BUILD_PVTU/petsc-3.5.2")
set(PETSC_ARCH "linux-gnu-debug")

# include libraries headers
include_directories(${INMOST_SRC}/Source/Headers)
include_directories(${INMOST_BUILD})
include_directories(${PETSC_DIR}/include)
include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)

# Add a library with the above sources
add_library(mesh STATIC 
  src/basis.cpp  
  src/mesh.cpp)

# link with variables and tools libraries
target_link_libraries(mesh PRIVATE 
  sub::variables)

  target_link_libraries(mesh PRIVATE ${INMOST_BUILD}/libinmost.a)
  target_link_libraries(mesh PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.a)
  target_link_libraries(mesh PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/lib/libparmetis.a)
  target_link_libraries(mesh PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.a)
  target_link_libraries(mesh PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2clapack.a)
  target_link_libraries(mesh PRIVATE ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2cblas.a)

add_library(sub::mesh ALIAS mesh)

target_include_directories( mesh
    PUBLIC ${PROJECT_SOURCE_DIR}/inc
)