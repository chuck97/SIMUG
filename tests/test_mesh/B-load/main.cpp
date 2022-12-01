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

bool test_mesh_load()
{
    IceMesh(MESH_PATH,
            mesh::surfType::plane,
            mesh::gridType::Agrid);
    
    return true;
}

int main()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (test_mesh_load())
    {
        if (rank == 0)
            std::cout << "Mesh load test: OK!\n";
    }
    else
        SIMUG_ERR("Mesh load test: FAILED!\n");

    BARRIER

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}