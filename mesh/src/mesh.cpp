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

    // find boundary elements, setup ids and id intervals
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
            bnd_edges.push_back(all_bnd_edges[i]);
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


void IceMesh::AssignIds()
{
    // assign global ids for nodes, edges and trians
    ice_mesh->AssignGlobalID(CELL|EDGE|FACE|NODE);


    // assign no_bnd ids for nodes
    grid_info[gridElemType::Node]->Create("node id no_bnd", 1, INMOST::DATA_INTEGER);
    int node_id = 0;
    for (int procn = 0; procn < ice_mesh->GetProcessorsNumber(); procn++)
    {
        if (ice_mesh->GetProcessorRank() == procn)
        {
            for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
            {
                if ((nodeit->GetStatus() != Element::Ghost) and
                    (nodeit->Integer(grid_info[gridElemType::Node]->Get("is node bnd")) == 0))
                {
                    nodeit->Integer(grid_info[gridElemType::Node]->Get("node id no_bnd")) = node_id;
                    ++node_id;
                }
            }
        }
        BARRIER
#ifdef USE_MPI
        if (ice_mesh->GetProcessorsNumber() > 1)
            MPI_Bcast(&node_id, 1, MPI_INT, procn, MPI_COMM_WORLD);        
#endif
    }

    // assign no_bnd ids for edges
    grid_info[gridElemType::Edge]->Create("edge id no_bnd", 1, INMOST::DATA_INTEGER);
    int edge_id = 0;
    for (int procn = 0; procn < ice_mesh->GetProcessorsNumber(); procn++)
    {
        if (ice_mesh->GetProcessorRank() == procn)
        {
            for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
            {
                if ((edgeit->GetStatus() != Element::Ghost) and
                    (edgeit->Integer(grid_info[gridElemType::Edge]->Get("is edge bnd")) == 0))
                {
                    edgeit->Integer(grid_info[gridElemType::Edge]->Get("edge id no_bnd")) = edge_id;
                    ++edge_id;
                }
            }
        }
        BARRIER
#ifdef USE_MPI
        if (ice_mesh->GetProcessorsNumber() > 1)
            MPI_Bcast(&edge_id, 1, MPI_INT, procn, MPI_COMM_WORLD);        
#endif
    }

    // assign no_bnd ids for trians
    grid_info[gridElemType::Trian]->Create("trian id no_bnd", 1, INMOST::DATA_INTEGER);
    int trian_id = 0;
    for (int procn = 0; procn < ice_mesh->GetProcessorsNumber(); procn++)
    {
        if (ice_mesh->GetProcessorRank() == procn)
        {
            for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
            {
                if ((trianit->GetStatus() != Element::Ghost) and
                    (trianit->Integer(grid_info[gridElemType::Trian]->Get("is trian bnd")) == 0))
                {
                    trianit->Integer(grid_info[gridElemType::Trian]->Get("trian id no_bnd")) = trian_id;
                    ++trian_id;
                }
            }
        }
        BARRIER
#ifdef USE_MPI
        if (ice_mesh->GetProcessorsNumber() > 1)
            MPI_Bcast(&trian_id, 1, MPI_INT, procn, MPI_COMM_WORLD);        
#endif
    }
    BARRIER
}

void IceMesh::AssignIdIntervals()
{
    // calculate global and no_bnd interval for nodes
    int nodeidmin = std::numeric_limits<int>::max();
    int nobnd_nodeidmin = std::numeric_limits<int>::max();
    int nodeidmax = std::numeric_limits<int>::min();
    int nobnd_nodeidmax = std::numeric_limits<int>::min();

    for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            int pid = nodeit->GlobalID();

            if(pid < nodeidmin)
                nodeidmin = pid;
            if((pid + 1) > nodeidmax)
                nodeidmax = pid + 1;

            if (nodeit->Integer(grid_info[gridElemType::Node]->Get("is node bnd")) == 0)
            {
                int nobnd_pid = nodeit->Integer(grid_info[gridElemType::Node]->Get("node id no_bnd"));
                if(nobnd_pid < nobnd_nodeidmin)
                    nobnd_nodeidmin = nobnd_pid;
                if((nobnd_pid + 1) > nobnd_nodeidmax)
                    nobnd_nodeidmax = nobnd_pid + 1;
            }
        }
    }
    mesh_info.id_interval_nodes = {nodeidmin, nodeidmax};
    mesh_info.id_interval_nodes_no_bnd = {nobnd_nodeidmin, nobnd_nodeidmax};
    mesh_info.num_nodes = nodeidmax - nodeidmin;
    mesh_info.num_bnd_nodes = bnd_nodes.size();
    BARRIER

    // calculate global and no_bnd interval for edges
    int edgeidmin = std::numeric_limits<int>::max();
    int nobnd_edgeidmin = std::numeric_limits<int>::max();
    int edgeidmax = std::numeric_limits<int>::min();
    int nobnd_edgeidmax = std::numeric_limits<int>::min();

    for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        if(edgeit->GetStatus() != Element::Ghost)
        {
            int pid = edgeit->GlobalID();;
            if(pid < edgeidmin)
                edgeidmin = pid;
            if((pid + 1) > edgeidmax)
                edgeidmax = pid + 1;
            
            if (edgeit->Integer(grid_info[gridElemType::Edge]->Get("is edge bnd")) == 0)
            {
                int nobnd_pid = edgeit->Integer(grid_info[gridElemType::Edge]->Get("edge id no_bnd"));
                if(nobnd_pid < nobnd_edgeidmin)
                    nobnd_edgeidmin = nobnd_pid;
                if((nobnd_pid + 1) > nobnd_edgeidmax)
                    nobnd_edgeidmax = nobnd_pid + 1;
            }
        }
    }
    mesh_info.id_interval_edges = {edgeidmin, edgeidmax};
    mesh_info.id_interval_edges_no_bnd = {nobnd_edgeidmin, nobnd_edgeidmax};
    mesh_info.num_edges = edgeidmax - edgeidmin;
    mesh_info.num_bnd_edges = bnd_edges.size();
    BARRIER

    // calculate global and no bnd interval for trians
    int trianidmin = std::numeric_limits<int>::max();
    int nobnd_trianidmin = std::numeric_limits<int>::max();
    int trianidmax = std::numeric_limits<int>::min();
    int nobnd_trianidmax = std::numeric_limits<int>::min();

    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        if(trianit->GetStatus() != Element::Ghost)
        {
            int pid = trianit->GlobalID();;
            if(pid < trianidmin)
            {
                trianidmin = pid;
            } 
            if((pid + 1) > trianidmax)
            {
                trianidmax = pid + 1;
            } 

            if (trianit->Integer(grid_info[gridElemType::Trian]->Get("is trian bnd")) == 0)
            {
                int nobnd_pid = trianit->Integer(grid_info[gridElemType::Trian]->Get("trian id no_bnd"));
                if(nobnd_pid < nobnd_trianidmin)
                    nobnd_trianidmin = nobnd_pid;
                if((nobnd_pid + 1) > nobnd_trianidmax)
                    nobnd_trianidmax = nobnd_pid + 1;
            }
        }
    }
    mesh_info.id_interval_trians = {trianidmin, trianidmax};
    mesh_info.id_interval_trians_no_bnd = {nobnd_trianidmin, nobnd_trianidmax};
    mesh_info.num_trians = trianidmax - trianidmin;
    mesh_info.num_bnd_trians = bnd_trians.size();
    BARRIER
};

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

    // make boundary identification tag for nodes, edges and triangles
    grid_info[gridElemType::Node] = make_shared<NodeData>(ice_mesh.get());
    grid_info[gridElemType::Node]->Create("is node bnd", 1, INMOST::DATA_INTEGER);
    for (size_t i = 0; i < bnd_nodes.size(); ++i)
        bnd_nodes[i].Integer(grid_info[gridElemType::Node]->Get("is node bnd")) = 1;

    grid_info[gridElemType::Edge] = make_shared<EdgeData>(ice_mesh.get());
    grid_info[gridElemType::Edge]->Create("is edge bnd", 1, INMOST::DATA_INTEGER);
    for (size_t i = 0; i < bnd_edges.size(); ++i)
        bnd_edges[i].Integer(grid_info[gridElemType::Edge]->Get("is edge bnd")) = 1;

    grid_info[gridElemType::Trian] = make_shared<TrianData>(ice_mesh.get());
    grid_info[gridElemType::Trian]->Create("is trian bnd", 1, INMOST::DATA_INTEGER);
    for (size_t i = 0; i < bnd_trians.size(); ++i)
        bnd_trians[i].Integer(grid_info[gridElemType::Trian]->Get("is trian bnd")) = 1; 
    BARRIER

    // assign ids 
    mesh_timer.Launch();
    AssignIds();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Assigning IDS succeessfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // assign id intervals
    mesh_timer.Launch();
    AssignIdIntervals();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Assigning ID intervals succeessfully! (" + to_string(duration) + " ms)\n");
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
