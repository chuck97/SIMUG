project(variables)

#add_library(variables INTERFACE)

add_library(variables STATIC 
  src/variable.cpp
#  src/model_config.cpp
#  src/phys_config.cpp
#  src/mesh_config.cpp
#  src/adv_config.cpp
#  src/dyn_config.cpp
#  src/forc_config.cpp
)
add_library(sub::variables ALIAS variables)

target_include_directories(variables
    PUBLIC ${PROJECT_SOURCE_DIR}/inc
)