#include "tests.hpp"

using namespace std;
using namespace SIMUG::mesh;
using namespace INMOST;

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