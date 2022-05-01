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

    // read mesh from file
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
        mesh_log.Log("Mesh: \'" + path_to_file_+ "\' loaded successfully! (" + to_string(duration) + " ms)\n");
    }
    BARRIER

    // make mesh partition
#ifdef USE_PARTITIONER
    mesh_timer.Launch();
    Partition();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Repartition finished successfully! (" + to_string(duration)+ " ms)\n");
#endif

    // find boundary elements and setup id intrvals
    ComputeMeshInfo();
}

void IceMesh::Partition()
{
	if (ice_mesh->GetProcessorsNumber() > 1) // need repartition
	{ 
	    unique_ptr<Partitioner> p = make_unique<INMOST::Partitioner>(ice_mesh.get());
#ifdef USE_PARTITIONER_PARMETIS
        p->SetMethod(Partitioner::Parmetis, Partitioner::Partition);
#elif USE_PARTITIONER_ZOLTAN
        p->SetMethod(Partitioner::Zoltan, Partitioner::Partition);
#else
        p->SetMethod(Partitioner::INNER_KMEANS, Partitioner::Partition);
#endif
        BARRIER
		p->Evaluate();

		ice_mesh->Redistribute();
		ice_mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
        ice_mesh->ExchangeGhost(1, NODE); 
        BARRIER

        ice_mesh->AssignGlobalID(CELL|EDGE|FACE|NODE);
        BARRIER
        ice_mesh->ExchangeGhost(1, NODE);
        BARRIER
	}
};

void IceMesh::SelectBndEdges()
{
    ElementArray<Face> all_bnd_edges = ice_mesh->GatherBoundaryFaces();
    for (size_t i = 0; i < all_bnd_edges.size(); ++i)
    {
        if (all_bnd_edges[i].GetStatus() != Element::Ghost)
        {
            bnd_edges.push_back(all_bnd_edges[i]);
        }
    }
    BARRIER
}

void IceMesh::SelectBndNodes()
{
    for (auto edgeit = bnd_edges.begin(); edgeit != bnd_edges.end(); ++edgeit)
    {
        ElementArray<INMOST::Node> ed_nodes = edgeit->getNodes();
        if (ed_nodes[0].GetStatus() != Element::Ghost)
        {
            bnd_nodes.push_back(ed_nodes[0]);
        }
    }
    BARRIER
}

void IceMesh::SelectBndTrians()
{
    for (auto edgeit = bnd_edges.begin(); edgeit != bnd_edges.end(); ++edgeit)
    {
        ElementArray<INMOST::Cell> ed_trians = edgeit->getCells();
        if (ed_trians[0].GetStatus() != Element::Ghost)
            bnd_trians.push_back(ed_trians[0]);
    }
    BARRIER
}


void IceMesh::ComputeMeshInfo()
{
    SIMUG::Logger mesh_log(cout);
    SIMUG::Timer mesh_timer;
    double duration;

    // select boundary edges
    mesh_timer.Launch();
    SelectBndEdges();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Boundary edges selected succeessfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // select boundary nodes
    mesh_timer.Launch();
    SelectBndNodes();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Boundary nodes selected succeessfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // select boundary triangles
    mesh_timer.Launch();
    SelectBndTrians();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Boundary triangles selected succeessfully! (" + to_string(duration) + " ms)\n");
    BARRIER
}

void IceMesh::SaveVTU(const std::string& filename) const
{
    SIMUG::Logger mesh_log(cout);
    SIMUG::Timer mesh_timer;
    double duration;

    std::string no_spaces_filename = filename;
    SIMUG::tools::rtrim(no_spaces_filename);
#ifdef USE_MPI
    if (ice_mesh->GetProcessorsNumber() > 1)
    {
        if(no_spaces_filename.substr(no_spaces_filename.size()-5, 5) != ".pvtu")
            SIMUG_ERR("Can't write mesh to this file. Filename should ended by .pvtu");
    }
    else
    {
        if(no_spaces_filename.substr(no_spaces_filename.size()-4, 4) != ".vtu")
            SIMUG_ERR("Can't write mesh to this file. Filename should ended by .vtu");
    }
#else
    if(no_spaces_filename.substr(no_spaces_filename.size()-4, 4) != ".vtu")
            SIMUG_ERR("Can't write mesh to this file. Filename should ended by .vtu");
#endif

    mesh_timer.Launch();
	ice_mesh->Save(filename);
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();

    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Mesh saved to \'" + filename + "\'! (" + to_string(duration) + " ms)\n");
    BARRIER
}

/*
void IceMesh::SelectNodes()
{
    INMOST::Tag is_bnd = ice_mesh->CreateTag("is node bnd", DATA_INTEGER, NODE, NONE, 1);

    // find all boundary nodes 
    ElementArray<Face> bnd_edges = ice_mesh->GatherBoundaryFaces();
    for (size_t i = 0; i < bnd_edges.size(); ++i)
    {
        ElementArray<Cell> adj_tr = bnd_edges[i].getCells();
        if ((adj_tr.size() == 1) and
            (adj_tr[0].GetStatus() != Element::Ghost))
        {
            ElementArray<Node> edg_nodes = bnd_edges[i].getNodes();
            edg_nodes[0]->Integer(data.IsNodeBnd) = 1;
            edg_nodes[1]->Integer(data.IsNodeBnd) = 1;
        }
    }
    BARRIER
    ice_mesh->ExchangeData(data.IsNodeBnd, NODE, 0);
    BARRIER

    // assign global id
    data.NodeIdGlobal = ice_mesh->CreateTag("node global id", DATA_INTEGER, NODE, NONE, 1);
    for (Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
         nodeit != ice_mesh->EndNode();
         ++nodeit)
    {
        if ((nodeit->GetStatus() != Element::Ghost))
        {
            nodeit->Integer(data.NodeIdGlobal) = nodeit->GlobalID();
        }
    }
    ice_mesh->ExchangeData(data.NodeIdGlobal, NODE, 0);
    BARRIER

    // assign global id without boundary nodes
    data.NodeIdNoBnd = ice_mesh->CreateTag("node global id no bnd", DATA_INTEGER, NODE, NONE, 1);
    int global_counter = 0;
    for (int curr_proc = 0;
         curr_proc < ice_mesh->GetProcessorsNumber();
         ++curr_proc)
    {
        if (ice_mesh->GetProcessorRank() == curr_proc)
        {
            for (Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
                nodeit != ice_mesh->EndNode();
                ++nodeit)
            {
                if ((nodeit->GetStatus() != Element::Ghost) and
                    (nodeit->Integer(data.IsNodeBnd) == 0))
                {
                    nodeit->Integer(data.NodeIdNoBnd) = global_counter;
                    ++global_counter;
                }
            }
        }
        BARRIER
#if defined(USE_MPI)
        MPI_Bcast(&global_counter, 1, MPI_INT, curr_proc, MPI_COMM_WORLD);
#endif
    }
    BARRIER
    ice_mesh->ExchangeData(data.NodeIdNoBnd, NODE, 0);

    if (!display)
    {
        ice_mesh->SetFileOption("Tag:is node bnd", "nosave");
        ice_mesh->SetFileOption("Tag:node global id", "nosave");
        ice_mesh->SetFileOption("Tag:node global id no bnd", "nosave");
    }
}
*/