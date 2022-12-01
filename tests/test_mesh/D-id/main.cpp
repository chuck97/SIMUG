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

bool test_id()
{
    SIMUG::Logger log(cout);

    IceMesh imesh(MESH_PATH, mesh::surfType::plane, mesh::gridType::Agrid);


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

    imesh.SaveVTU("bndtest");
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