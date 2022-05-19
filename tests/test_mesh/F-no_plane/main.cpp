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

#define SPHERE_PATH "../../../../SIMUG_v0/MESHES/pmf/Sphere.pmf"

using namespace INMOST;
using namespace SIMUG::mesh;
using namespace std;

bool test_not_plane()
{
    //IceMesh mesh_arctic(ARCTIC_PATH, surfType::basin, gridType::Agrid);
    IceMesh mesh_sphere(SPHERE_PATH, surfType::sphere, gridType::Agrid);

    //if (mesh_arctic.GetMesh()->GetProcessorsNumber() > 1)
    //    mesh_arctic.SaveVTU("./arctic.pvtu");
    //else
     //   mesh_arctic.SaveVTU("./arctic.vtu");

    mesh_sphere.SaveVTU("sphere");
    return true;
}

int main()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (test_not_plane())
    {
        if (rank == 0)
            std::cout << "Not plane test: OK!\n";
    }
    else
        SIMUG_ERR("Not plane test: FAILED!\n");

    BARRIER

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}