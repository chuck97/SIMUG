# Set the project name
project(configuration)

# Find nlohmann::json
find_package(nlohmann_json 3.10.5 REQUIRED)

# Add a library with the above sources
add_library(configuration STATIC 
  src/config.cpp
  src/model_config.cpp
  src/phys_config.cpp
  src/mesh_config.cpp
  src/adv_config.cpp
  src/dyn_config.cpp
  src/forc_config.cpp
)

# link with variables and tools libraries
target_link_libraries(configuration PRIVATE 
  sub::variables 
  sub::tools)

target_link_libraries(configuration PUBLIC 
  nlohmann_json::nlohmann_json)

# Make alias for library
add_library(sub::configuration ALIAS configuration)

# include "inc" directory
target_include_directories( configuration
    PUBLIC ${PROJECT_SOURCE_DIR}/inc
)