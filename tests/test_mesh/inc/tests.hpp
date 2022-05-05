#pragma once

#include "defines.hpp"
#include "inmost.h"
#include "model_var.hpp"
#include "mesh_info.hpp"
#include "data.hpp"
#include "mesh.hpp"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <memory>
#include <vector>
#include <string>
#include <iostream>

#define MESH_PATH "/home/users/spetrov/SIMUG/SIMUG_v0/MESHES/pmf/square8km.pmf"

bool test_data();
bool test_mesh_load();
bool test_bnd_selection();
bool test_id();
bool test_mute();