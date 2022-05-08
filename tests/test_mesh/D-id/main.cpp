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
#define ARCTIC_PATH "/home/users/spetrov/SIMUG/SIMUG_v0/MESHES/pmf/Arctic.pmf"
#define SPHERE_PATH "/home/users/spetrov/SIMUG/SIMUG_v0/MESHES/pmf/Sphere.pmf"

using namespace INMOST;
using namespace SIMUG::mesh;
using namespace std;

bool test_id()
{
    SIMUG::Logger log(cout);

    IceMesh imesh(MESH_PATH, surfType::plane, gridType::Agrid);


    for (int pr = 0; pr < imesh.GetMesh()->GetProcessorsNumber(); ++pr)
    {
        if (pr == imesh.GetMesh()->GetProcessorRank())
        {
            log.Log("==============================================\n");
            log.Log("Processor " + to_string(imesh.GetMesh()->GetProcessorRank()) + ": \n");
            log.Log("Num nodes = " + to_string(imesh.GetMeshInfo().num_nodes) + ";\n");
            log.Log("Num edges = " + to_string(imesh.GetMeshInfo().num_edges) + ";\n");
            log.Log("Num trians = " + to_string(imesh.GetMeshInfo().num_trians) + ";\n");
            log.Log("Num bnd nodes = " + to_string(imesh.GetMeshInfo().num_bnd_nodes) + ";\n");
            log.Log("Num bnd edges = " + to_string(imesh.GetMeshInfo().num_bnd_edges) + ";\n");
            log.Log("Num bnd trians = " + to_string(imesh.GetMeshInfo().num_bnd_trians) + ";\n");
            log.Log("Id interval nodes: " + to_string(imesh.GetMeshInfo().id_interval_nodes.id_min) + " -- " + to_string(imesh.GetMeshInfo().id_interval_nodes.id_max) + ";\n");
            log.Log("Id interval edges: " + to_string(imesh.GetMeshInfo().id_interval_edges.id_min) + " -- " + to_string(imesh.GetMeshInfo().id_interval_edges.id_max) + ";\n");
            log.Log("Id interval trians: " + to_string(imesh.GetMeshInfo().id_interval_trians.id_min) + " -- " + to_string(imesh.GetMeshInfo().id_interval_trians.id_max) + ";\n");
            log.Log("Id interval nodes no_bnd: " + to_string(imesh.GetMeshInfo().id_interval_nodes_no_bnd.id_min) + " -- " + to_string(imesh.GetMeshInfo().id_interval_nodes_no_bnd.id_max) + ";\n");
            log.Log("Id interval edges no_bnd: " + to_string(imesh.GetMeshInfo().id_interval_edges_no_bnd.id_min) + " -- " + to_string(imesh.GetMeshInfo().id_interval_edges_no_bnd.id_max) + ";\n");
            log.Log("Id interval trians no_bnd: " + to_string(imesh.GetMeshInfo().id_interval_trians_no_bnd.id_min) + " -- " + to_string(imesh.GetMeshInfo().id_interval_trians_no_bnd.id_max) + ";\n");
            log.Log("==============================================\n");
        }
        BARRIER
    }

    if (imesh.GetMesh()->GetProcessorsNumber() > 1)
        imesh.SaveVTU("./bndtest.pvtu");
    else
        imesh.SaveVTU("./bndtest.vtu");

    return true;
}


int main()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (test_id())
    {
        if (rank == 0)
            std::cout << "Id computation test: OK!\n";
    }
    else
        SIMUG_ERR("Id computation test: FAILED!\n");

    BARRIER

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}