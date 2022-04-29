#include "mesh.hpp"

using namespace std;
using namespace INMOST;
using namespace SIMUG::mesh;

IceMesh::IceMesh(const std::string& path_to_file_,
                 const surfType& surf_type_,
                 const gridType& grid_type_) 
{
    // logger and timer initialization
    SIMUG::Logger mesh_log(cout);
    SIMUG::Timer mesh_timer;
    double duration;

    // update mesh information
    mesh_info.surface_type = surf_type_;
    mesh_info.grid_type = grid_type_;
    num_ice_layers = 1;

    // INMOST::mesh initialization and creating new mesh
    Mesh::Initialize(NULL, NULL);
    ice_mesh = make_shared<INMOST::Mesh>();
    
#if defined(USE_MPI)
    ice_mesh->SetCommunicator(INMOST_MPI_COMM_WORLD);
#endif

    // setting partitioner
#ifdef USE_PARTITIONER
    Partitioner::Initialize(NULL, NULL);
#endif

    // check the existance of mesh .pmf file
    if (!std::filesystem::exists(path_to_file_))
    {
        if (ice_mesh->GetProcessorRank()==0)
            SIMUG_ERR("ERROR! Mesh file \'" + path_to_file_ + "\' doesn't exist");

    }

    mesh_timer.Launch();

    if(ice_mesh->isParallelFileFormat(path_to_file_))
        ice_mesh->Load(path_to_file_); 
    else
    {
        if (ice_mesh->GetProcessorRank() == 0)
            ice_mesh->Load(path_to_file_);
    }

    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();

    if (ice_mesh->GetProcessorRank()==0)
    {
        mesh_log.Log("================== Mesh initialization ==================\n");
        mesh_log.Log("Mesh: \'" + path_to_file_+ "\' loaded successfully!\n");
        mesh_log.Log("Mesh initialization duration: " + to_string(duration) + " ms;\n");
    }
    BARRIER
}