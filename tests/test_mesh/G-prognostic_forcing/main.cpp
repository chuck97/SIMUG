#include "simug.hpp"
#include "inmost.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <memory>
#include <vector>
#include <string>
#include <iostream>

#define MESH_PATH "../../../../SIMUG_v0/MESHES/pmf/Box_low_res.pmf"

using namespace INMOST;
using namespace SIMUG;
using namespace std;

bool test_prognostic_forcing()
{
    IceMesh mesh_plane(MESH_PATH, mesh::surfType::plane, mesh::gridType::Agrid, 5);
    mesh_plane.SaveVTU("mult_layers");

    return true;
}

int main()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (test_prognostic_forcing())
    {
        if (rank == 0)
            std::cout << "Prognostic_forcing test: OK!\n";
    }
    else
        SIMUG_ERR("Prognostic_forcing test: FAILED!\n");

    BARRIER

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}