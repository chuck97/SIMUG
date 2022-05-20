#include "mesh.hpp"

using namespace std;
using namespace INMOST;
using namespace SIMUG;

IceMesh::IceMesh(const std::string& path_to_file_,
                 const std::string& output_folder_,
                 const mesh::surfType& surf_type_,
                 const mesh::gridType& grid_type_,
                 const int& n_ice_layers_) 
{
    // logger and timer initialization
    SIMUG::Logger mesh_log(cout);
    SIMUG::Timer mesh_timer;
    double duration;

    // update mesh information
    mesh_info.surface_type = surf_type_;
    mesh_info.grid_type = grid_type_;
    mesh_info.num_ice_layers = n_ice_layers_;
    mesh_info.output_folder = output_folder_;  

    if (mesh_info.grid_type == mesh::gridType::Agrid)
    {
        mesh_info.prog_elems = mesh::gridA_progElems;
        mesh_info.forc_elems = mesh::gridA_forcElems;
    }
    else if (mesh_info.grid_type == mesh::gridType::Bgrid)
    {
        SIMUG_ERR("B grid is not implemented yet");
    }
    else if (mesh_info.grid_type == mesh::gridType::Cgrid)
    {
        SIMUG_ERR("C grid is not implemented yet");
    }
    else
        SIMUG_ERR("Unknown type of grid! Possible options: Agrid, Bgrid, Cgrid");

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
    grid_info[mesh::gridElemType::Node] = make_shared<NodeInfo>(ice_mesh.get());
    grid_info[mesh::gridElemType::Edge] = make_shared<EdgeInfo>(ice_mesh.get());
    grid_info[mesh::gridElemType::Trian] = make_shared<TrianInfo>(ice_mesh.get());

    mesh_timer.Launch();
    ComputeMeshInfo();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();

    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Computing mesh information successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // assign coords for grid elements
    mesh_timer.Launch();
    AssignCoords();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();

    // exchange and mute all grid info
    for (auto [key, val]: grid_info)
        val->Mute();
        
    BARRIER

    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Assigning coordinates successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // assign prognostic model grid variables
    for (int lay = 0; lay < mesh_info.num_ice_layers; ++lay)
    {
        if (lay == 0)
        {
            prognostic_data[mesh::gridElemType::Node][lay] = make_shared<NodeData>(ice_mesh.get());
            prognostic_data[mesh::gridElemType::Edge][lay] = make_shared<EdgeData>(ice_mesh.get());
            prognostic_data[mesh::gridElemType::Trian][lay] = make_shared<TrianData>(ice_mesh.get());
            prognostic_data[mesh::gridElemType::bndNode][lay] = make_shared<BndNodeData>(ice_mesh.get());
            prognostic_data[mesh::gridElemType::bndEdge][lay] = make_shared<BndEdgeData>(ice_mesh.get());
            prognostic_data[mesh::gridElemType::bndTrian][lay] = make_shared<BndTrianData>(ice_mesh.get());
        }
        else
        {
            prognostic_data[mesh::gridElemType::Node][lay] = make_shared<NodeData>(ice_mesh.get(), lay);
            prognostic_data[mesh::gridElemType::Edge][lay] = make_shared<EdgeData>(ice_mesh.get(), lay);
            prognostic_data[mesh::gridElemType::Trian][lay] = make_shared<TrianData>(ice_mesh.get(), lay);
            prognostic_data[mesh::gridElemType::bndNode][lay] = make_shared<BndNodeData>(ice_mesh.get(), lay);
            prognostic_data[mesh::gridElemType::bndEdge][lay] = make_shared<BndEdgeData>(ice_mesh.get(), lay);
            prognostic_data[mesh::gridElemType::bndTrian][lay] = make_shared<BndTrianData>(ice_mesh.get(), lay);
        }
    }
    BARRIER

    mesh_timer.Launch();
    AssignPrognosticVariables();
    MuteProgData();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Assigning prognostic variables successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // assign forcing grid variables
    forcing_data[mesh::gridElemType::Node] = make_shared<NodeData>(ice_mesh.get());
    forcing_data[mesh::gridElemType::Edge] = make_shared<EdgeData>(ice_mesh.get());
    forcing_data[mesh::gridElemType::Trian] = make_shared<TrianData>(ice_mesh.get());
    forcing_data[mesh::gridElemType::bndNode] = make_shared<BndNodeData>(ice_mesh.get());
    forcing_data[mesh::gridElemType::bndEdge] = make_shared<BndEdgeData>(ice_mesh.get());
    forcing_data[mesh::gridElemType::bndTrian] = make_shared<BndTrianData>(ice_mesh.get());
    
    BARRIER

    mesh_timer.Launch();
    AssignForcingVariables();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Assigning forcing variables successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    mesh_timer.Launch();
    AssembleBasisData();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Assembling basis data successfully! (" + to_string(duration) + " ms)\n");

    BARRIER
}

IceMesh::IceMesh(const std::string& path_to_file_,
                 const mesh::surfType& surf_type_,
                 const mesh::gridType& grid_type_,
                 const int& n_ice_layers_):
        IceMesh(path_to_file_, "./", surf_type_, grid_type_, n_ice_layers_)
{}

IceMesh::IceMesh(const std::string& path_to_file_,
                 const std::string& output_folder_,
                 const mesh::surfType& surf_type_,
                 const mesh::gridType& grid_type_):
        IceMesh(path_to_file_, output_folder_, surf_type_, grid_type_, 1)
{}

IceMesh::IceMesh(const std::string& path_to_file_,
                 const mesh::surfType& surf_type_,
                 const mesh::gridType& grid_type_):
        IceMesh(path_to_file_, "./", surf_type_, grid_type_, 1)
{}



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
	}
};

void IceMesh::SelectBndEdges()
{
    ElementArray<Face> all_bnd_edges = ice_mesh->GatherBoundaryFaces();
    for (size_t i = 0; i < all_bnd_edges.size(); ++i)
    {
        bnd_edges.push_back(all_bnd_edges[i]);
    }
    BARRIER
}

void IceMesh::SelectBndEdgesNoGhost()
{
    ElementArray<Face> all_bnd_edges = ice_mesh->GatherBoundaryFaces();
    bnd_edges.clear();
    for (size_t i = 0; i < all_bnd_edges.size(); ++i)
    {
        if (all_bnd_edges[i].GetStatus() != Element::Ghost)
            bnd_edges.push_back(all_bnd_edges[i]);
    }
    BARRIER
}

void IceMesh::SelectBndNodes()
{
    for (size_t i = 0; i < bnd_edges.size(); ++i)
    {
        ElementArray<INMOST::Node> ed_nodes = bnd_edges[i].getNodes();

        if (ed_nodes[0].GetStatus() != Element::Ghost)
        {
            bnd_nodes.push_back(ed_nodes[0]);
        }
    }
    BARRIER
}

void IceMesh::SelectBndTrians()
{
    for (size_t i = 0; i < bnd_edges.size(); ++i)
    {
        ElementArray<INMOST::Cell> ed_trians = bnd_edges[i].getCells();
        if (ed_trians[0].GetStatus() != Element::Ghost)
            bnd_trians.push_back(ed_trians[0]);
    }
    BARRIER
}


void IceMesh::AssignIds()
{
    // assign global ids for nodes, edges and trians
    ice_mesh->AssignGlobalID(CELL|EDGE|FACE|NODE);

    // assign id and no_bnd id for nodes
    grid_info[mesh::gridElemType::Node]->id = ice_mesh->CreateTag("id node", INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);
    grid_info[mesh::gridElemType::Node]->id_no_bnd = ice_mesh->CreateTag("id node no bnd", INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);

    int node_id = 0;
    for (int procn = 0; procn < ice_mesh->GetProcessorsNumber(); procn++)
    {
        if (ice_mesh->GetProcessorRank() == procn)
        {
            for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
            {
                if (nodeit->GetStatus() != Element::Ghost)
                    nodeit->Integer(grid_info[mesh::gridElemType::Node]->id) = nodeit->GlobalID();

                if ((nodeit->GetStatus() != Element::Ghost) and
                    (nodeit->Integer(grid_info[mesh::gridElemType::Node]->is_bnd) == 0))
                {
                    nodeit->Integer(grid_info[mesh::gridElemType::Node]->id_no_bnd) = node_id;
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

    // exchange node ids
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Node]->id, INMOST::NODE, 0);
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Node]->id_no_bnd, INMOST::NODE, 0);

    // assign id and no_bnd id for edges
    grid_info[mesh::gridElemType::Edge]->id = ice_mesh->CreateTag("id edge", INMOST::DATA_INTEGER, INMOST::FACE, INMOST::NONE, 1);
    grid_info[mesh::gridElemType::Edge]->id_no_bnd = ice_mesh->CreateTag("id edge no bnd", INMOST::DATA_INTEGER, INMOST::FACE, INMOST::NONE, 1);

    int edge_id = 0;
    for (int procn = 0; procn < ice_mesh->GetProcessorsNumber(); procn++)
    {
        if (ice_mesh->GetProcessorRank() == procn)
        {
            for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
            {
                if (edgeit->GetStatus() != Element::Ghost)
                    edgeit->Integer(grid_info[mesh::gridElemType::Edge]->id) = edgeit->GlobalID();

                if ((edgeit->GetStatus() != Element::Ghost) and
                    (edgeit->Integer(grid_info[mesh::gridElemType::Edge]->is_bnd) == 0))
                {
                    edgeit->Integer(grid_info[mesh::gridElemType::Edge]->id_no_bnd) = edge_id;
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

    // exchange edge ids
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Edge]->id, INMOST::FACE, 0);
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Edge]->id_no_bnd, INMOST::FACE, 0);

    // assign no_bnd ids for trians
    grid_info[mesh::gridElemType::Trian]->id = ice_mesh->CreateTag("id trian", INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 1);
    grid_info[mesh::gridElemType::Trian]->id_no_bnd = ice_mesh->CreateTag("id trian no bnd", INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 1);

    int trian_id = 0;
    for (int procn = 0; procn < ice_mesh->GetProcessorsNumber(); procn++)
    {
        if (ice_mesh->GetProcessorRank() == procn)
        {
            for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
            {
                if (trianit->GetStatus() != Element::Ghost)
                trianit->Integer(grid_info[mesh::gridElemType::Trian]->id) = trianit->GlobalID();

                if ((trianit->GetStatus() != Element::Ghost) and
                    (trianit->Integer(grid_info[mesh::gridElemType::Trian]->is_bnd) == 0))
                {
                    trianit->Integer(grid_info[mesh::gridElemType::Trian]->id_no_bnd) = trian_id;
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

    // exchange trian ids
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Trian]->id, INMOST::CELL, 0);
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Trian]->id_no_bnd, INMOST::CELL, 0);

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
            int pid = nodeit->Integer(grid_info[mesh::gridElemType::Node]->id);

            if(pid < nodeidmin)
                nodeidmin = pid;
            if((pid + 1) > nodeidmax)
                nodeidmax = pid + 1;

            if (nodeit->Integer(grid_info[mesh::gridElemType::Node]->is_bnd) == 0)
            {
                int nobnd_pid = nodeit->Integer(grid_info[mesh::gridElemType::Node]->id_no_bnd);
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
            int pid = edgeit->Integer(grid_info[mesh::gridElemType::Edge]->id);
            if(pid < edgeidmin)
                edgeidmin = pid;
            if((pid + 1) > edgeidmax)
                edgeidmax = pid + 1;
            
            if (edgeit->Integer(grid_info[mesh::gridElemType::Edge]->is_bnd) == 0)
            {
                int nobnd_pid = edgeit->Integer(grid_info[mesh::gridElemType::Edge]->id_no_bnd);
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
            int pid = trianit->Integer(grid_info[mesh::gridElemType::Trian]->id);
            if(pid < trianidmin)
            {
                trianidmin = pid;
            } 
            if((pid + 1) > trianidmax)
            {
                trianidmax = pid + 1;
            } 

            if (trianit->Integer(grid_info[mesh::gridElemType::Trian]->is_bnd) == 0)
            {
                int nobnd_pid = trianit->Integer(grid_info[mesh::gridElemType::Trian]->id_no_bnd);
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
        mesh_log.Log("Boundary edges selected successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // select boundary nodes
    mesh_timer.Launch();
    SelectBndNodes();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Boundary nodes selected successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // select boundary triangles
    mesh_timer.Launch();
    SelectBndTrians();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Boundary triangles selected successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // select boundary edges no ghost
    mesh_timer.Launch();
    SelectBndEdgesNoGhost();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Boundary edges (no ghost) selected successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // make boundary identification tag for nodes, edges and triangles
    grid_info[mesh::gridElemType::Node]->is_bnd = ice_mesh->CreateTag("is node bnd", INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);
    for (size_t i = 0; i < bnd_nodes.size(); ++i)
        bnd_nodes[i].Integer(grid_info[mesh::gridElemType::Node]->is_bnd) = 1;
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Node]->is_bnd, INMOST::NODE, 0);

    grid_info[mesh::gridElemType::Edge]->is_bnd = ice_mesh->CreateTag("is edge bnd", INMOST::DATA_INTEGER, INMOST::FACE, INMOST::NONE, 1);
    for (size_t i = 0; i < bnd_edges.size(); ++i)
        bnd_edges[i].Integer(grid_info[mesh::gridElemType::Edge]->is_bnd) = 1;
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Edge]->is_bnd, INMOST::FACE, 0);

    grid_info[mesh::gridElemType::Trian]->is_bnd = ice_mesh->CreateTag("is trian bnd", INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 1);
    for (size_t i = 0; i < bnd_trians.size(); ++i)
        bnd_trians[i].Integer(grid_info[mesh::gridElemType::Trian]->is_bnd) = 1; 
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Trian]->is_bnd, INMOST::CELL, 0);

    BARRIER

    // assign ids 
    mesh_timer.Launch();
    AssignIds();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Assigning IDS successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // assign id intervals
    mesh_timer.Launch();
    AssignIdIntervals();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();
    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Assigning ID intervals successfully! (" + to_string(duration) + " ms)\n");
    BARRIER
}

void IceMesh::AssignCoords()
{
    // ## model coords ##

    // nodes
    INMOST::Tag node_coords_tag = ice_mesh->CreateTag("model coords node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    grid_info[mesh::gridElemType::Node]->coords[coord::coordType::model] = node_coords_tag;

    for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            nodeit->RealArray(node_coords_tag)[0] = nodeit->Coords()[0];
            nodeit->RealArray(node_coords_tag)[1] = nodeit->Coords()[1];
            nodeit->RealArray(node_coords_tag)[2] = (mesh_info.surface_type == mesh::surfType::sphere)?nodeit->Coords()[2]:0.0;
        }
    }
    ice_mesh->ExchangeData(node_coords_tag, INMOST::NODE, 0);

    // edges
    INMOST::Tag edge_coords_tag = ice_mesh->CreateTag("model coords edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    grid_info[mesh::gridElemType::Edge]->coords[coord::coordType::model] = edge_coords_tag;

    for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        if (edgeit->GetStatus() != Element::Ghost)
        {
            ElementArray<INMOST::Node> adj_nodes = edgeit->getNodes();
            edgeit->RealArray(edge_coords_tag)[0] = 0.5*(adj_nodes[0].Coords()[0] + adj_nodes[1].Coords()[0]);
            edgeit->RealArray(edge_coords_tag)[1] = 0.5*(adj_nodes[0].Coords()[1] + adj_nodes[1].Coords()[1]);
            edgeit->RealArray(edge_coords_tag)[2] = (mesh_info.surface_type == mesh::surfType::sphere)?0.5*(adj_nodes[0].Coords()[2] +
                                                                                                       adj_nodes[1].Coords()[2]):0.0;
        }
    }
    ice_mesh->ExchangeData(edge_coords_tag, INMOST::FACE, 0);

    // trians
    INMOST::Tag trian_coords_tag = ice_mesh->CreateTag("model coords trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    grid_info[mesh::gridElemType::Trian]->coords[coord::coordType::model] = trian_coords_tag;

    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();
            trianit->RealArray(trian_coords_tag)[0] = (adj_nodes[0].Coords()[0] + adj_nodes[1].Coords()[0] + adj_nodes[2].Coords()[0])/3.0;
            trianit->RealArray(trian_coords_tag)[1] = (adj_nodes[0].Coords()[1] + adj_nodes[1].Coords()[1] + adj_nodes[2].Coords()[1])/3.0;
            trianit->RealArray(trian_coords_tag)[2] = (mesh_info.surface_type == mesh::surfType::sphere)?(adj_nodes[0].Coords()[2] +
                                                                                                     adj_nodes[1].Coords()[2] +
                                                                                                     adj_nodes[2].Coords()[2])/3.0:0.0;
        }
    }
    ice_mesh->ExchangeData(trian_coords_tag, INMOST::CELL, 0);

    BARRIER

    // ## cartesian coords ##

    // nodes 
    INMOST::Tag node_cart_tag = ice_mesh->CreateTag("cart coords node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart] = node_cart_tag;

    for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            double model_x = nodeit->RealArray(node_coords_tag)[0];
            double model_y = nodeit->RealArray(node_coords_tag)[1];
            double model_z = nodeit->RealArray(node_coords_tag)[2];

            nodeit->RealArray(node_cart_tag)[0] = (mesh_info.surface_type == mesh::surfType::basin)?cos(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                              :model_x;
            nodeit->RealArray(node_cart_tag)[1] = (mesh_info.surface_type == mesh::surfType::basin)?sin(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                              :model_y;
            nodeit->RealArray(node_cart_tag)[2] = (mesh_info.surface_type == mesh::surfType::basin)?sin(model_y*(M_PI/180.0))
                                                                                              :model_z;
        }
    }
    ice_mesh->ExchangeData(node_cart_tag, INMOST::NODE, 0);

    // edges
    INMOST::Tag edge_cart_tag = ice_mesh->CreateTag("cart coords edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    grid_info[mesh::gridElemType::Edge]->coords[coord::coordType::cart] = edge_cart_tag;

    for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        if (edgeit->GetStatus() != Element::Ghost)
        {
            double model_x = edgeit->RealArray(edge_coords_tag)[0];
            double model_y = edgeit->RealArray(edge_coords_tag)[1];
            double model_z = edgeit->RealArray(edge_coords_tag)[2];

            edgeit->RealArray(edge_cart_tag)[0] = (mesh_info.surface_type == mesh::surfType::basin)?cos(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                              :model_x;
            edgeit->RealArray(edge_cart_tag)[1] = (mesh_info.surface_type == mesh::surfType::basin)?sin(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                              :model_y;
            edgeit->RealArray(edge_cart_tag)[2] = (mesh_info.surface_type == mesh::surfType::basin)?sin(model_y*(M_PI/180.0))
                                                                                              :model_z;
        }
    }
    ice_mesh->ExchangeData(edge_cart_tag, INMOST::FACE, 0);

    // triangles
    INMOST::Tag trian_cart_tag = ice_mesh->CreateTag("cart coords trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    grid_info[mesh::gridElemType::Trian]->coords[coord::coordType::cart] = trian_cart_tag;

    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            double model_x = trianit->RealArray(trian_coords_tag)[0];
            double model_y = trianit->RealArray(trian_coords_tag)[1];
            double model_z = trianit->RealArray(trian_coords_tag)[2];

            trianit->RealArray(trian_cart_tag)[0] = (mesh_info.surface_type == mesh::surfType::basin)?cos(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                              :model_x;
            trianit->RealArray(trian_cart_tag)[1] = (mesh_info.surface_type == mesh::surfType::basin)?sin(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                              :model_y;
            trianit->RealArray(trian_cart_tag)[2] = (mesh_info.surface_type == mesh::surfType::basin)?sin(model_y*(M_PI/180.0))
                                                                                              :model_z;
        }
    }
    ice_mesh->ExchangeData(trian_cart_tag, INMOST::CELL, 0);

    BARRIER

    // ## geo coords ##

    // nodes 
    INMOST::Tag node_geo_tag = ice_mesh->CreateTag("geo coords node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    grid_info[mesh::gridElemType::Node]->coords[coord::coordType::geo] = node_geo_tag;

    for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            if (mesh_info.surface_type == mesh::surfType::sphere)
            {
                double x = nodeit->RealArray(node_cart_tag)[0];
                double y = nodeit->RealArray(node_cart_tag)[1];
                double z = nodeit->RealArray(node_cart_tag)[2];

                double r = 1;
                double lat, lon;


                if ((x >= 0) and (y >= 0))
                {
                    lat = std::asin(z/r);
                    lon = std::atan(y/x);
                }
                else if ((x >= 0) and (y <= 0))
                {
                    lat = std::asin(z/r);
                    lon = 2*M_PI + std::atan(y/x);
                }
                else if ((x <= 0) and (y >= 0))
                {
                    lat = std::asin(z/r);
                    lon = M_PI + std::atan(y/x);
                }
                else
                {
                    lat = std::asin(z/r);
                    lon = M_PI + std::atan(y/x);
                }

                nodeit->RealArray(node_geo_tag)[0] = lon;
                nodeit->RealArray(node_geo_tag)[1] = lat;
                nodeit->RealArray(node_geo_tag)[2] = 1.0;  
            }
            else if (mesh_info.surface_type == mesh::surfType::basin)
            {
                nodeit->RealArray(node_geo_tag)[0] = nodeit->RealArray(node_coords_tag)[0]*(M_PI/180.0);
                nodeit->RealArray(node_geo_tag)[1] = nodeit->RealArray(node_coords_tag)[1]*(M_PI/180.0);
                nodeit->RealArray(node_geo_tag)[2] = 1.0;
            }
            else
            {
                nodeit->RealArray(node_geo_tag)[0] = 0.0;
                nodeit->RealArray(node_geo_tag)[1] = 0.0;
                nodeit->RealArray(node_geo_tag)[2] = 0.0;
            }
        }
    }   
    ice_mesh->ExchangeData(node_geo_tag, INMOST::NODE, 0);

    // edges
    INMOST::Tag edge_geo_tag = ice_mesh->CreateTag("geo coords edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    grid_info[mesh::gridElemType::Edge]->coords[coord::coordType::geo] = edge_geo_tag;

    for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        if (edgeit->GetStatus() != Element::Ghost)
        {
            if (mesh_info.surface_type == mesh::surfType::sphere)
            {
                double x = edgeit->RealArray(edge_cart_tag)[0];
                double y = edgeit->RealArray(edge_cart_tag)[1];
                double z = edgeit->RealArray(edge_cart_tag)[2];

                double r = 1;
                double lat, lon;


                if ((x >= 0) and (y >= 0))
                {
                    lat = std::asin(z/r);
                    lon = std::atan(y/x);
                }
                else if ((x >= 0) and (y <= 0))
                {
                    lat = std::asin(z/r);
                    lon = 2*M_PI + std::atan(y/x);
                }
                else if ((x <= 0) and (y >= 0))
                {
                    lat = std::asin(z/r);
                    lon = M_PI + std::atan(y/x);
                }
                else
                {
                    lat = std::asin(z/r);
                    lon = M_PI + std::atan(y/x);
                }

                edgeit->RealArray(edge_geo_tag)[0] = lon;
                edgeit->RealArray(edge_geo_tag)[1] = lat;
                edgeit->RealArray(edge_geo_tag)[2] = 1.0;  
            }
            else if (mesh_info.surface_type == mesh::surfType::basin)
            {
                edgeit->RealArray(edge_geo_tag)[0] = edgeit->RealArray(edge_coords_tag)[0]*M_PI/180.0;
                edgeit->RealArray(edge_geo_tag)[1] = edgeit->RealArray(edge_coords_tag)[1]*M_PI/180.0;
                edgeit->RealArray(edge_geo_tag)[2] = 1.0;
            }
            else
            {
                edgeit->RealArray(edge_geo_tag)[0] = 0.0;
                edgeit->RealArray(edge_geo_tag)[1] = 0.0;
                edgeit->RealArray(edge_geo_tag)[2] = 0.0;
            }
        }
    }
    ice_mesh->ExchangeData(edge_geo_tag, INMOST::FACE, 0);

    // trians
    INMOST::Tag trian_geo_tag = ice_mesh->CreateTag("geo coords trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    grid_info[mesh::gridElemType::Trian]->coords[coord::coordType::geo] = trian_geo_tag;

    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            if (mesh_info.surface_type == mesh::surfType::sphere)
            {
                double x = trianit->RealArray(trian_cart_tag)[0];
                double y = trianit->RealArray(trian_cart_tag)[1];
                double z = trianit->RealArray(trian_cart_tag)[2];

                double r = 1;
                double lat, lon;


                if ((x >= 0) and (y >= 0))
                {
                    lat = std::asin(z/r);
                    lon = std::atan(y/x);
                }
                else if ((x >= 0) and (y <= 0))
                {
                    lat = std::asin(z/r);
                    lon = 2*M_PI + std::atan(y/x);
                }
                else if ((x <= 0) and (y >= 0))
                {
                    lat = std::asin(z/r);
                    lon = M_PI + std::atan(y/x);
                }
                else
                {
                    lat = std::asin(z/r);
                    lon = M_PI + std::atan(y/x);
                }

                trianit->RealArray(trian_geo_tag)[0] = lon;
                trianit->RealArray(trian_geo_tag)[1] = lat;
                trianit->RealArray(trian_geo_tag)[2] = 1.0;  
            }
            else if (mesh_info.surface_type == mesh::surfType::basin)
            {
                trianit->RealArray(trian_geo_tag)[0] = trianit->RealArray(trian_coords_tag)[0]*M_PI/180.0;
                trianit->RealArray(trian_geo_tag)[1] = trianit->RealArray(trian_coords_tag)[1]*M_PI/180.0;
                trianit->RealArray(trian_geo_tag)[2] = 1.0;
            }
            else
            {
                trianit->RealArray(trian_geo_tag)[0] = 0.0;
                trianit->RealArray(trian_geo_tag)[1] = 0.0;
                trianit->RealArray(trian_geo_tag)[2] = 0.0;
            }
        }
    }
    ice_mesh->ExchangeData(trian_geo_tag, INMOST::CELL, 0);

    BARRIER
}

void IceMesh::SaveVTU(const std::string& meshname) const
{
    SIMUG::Logger mesh_log(cout);
    SIMUG::Timer mesh_timer;
    double duration;

    std::string no_spaces_meshname = meshname;
    SIMUG::tools::rtrim(no_spaces_meshname);

    mesh_timer.Launch();
    no_spaces_meshname += (ice_mesh->GetProcessorsNumber() == 1) ? ".vtu": ".pvtu";
    ice_mesh->Save(mesh_info.output_folder + no_spaces_meshname);
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();

    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Mesh saved to \'" + mesh_info.output_folder + no_spaces_meshname + "\'! (" + to_string(duration) + " ms)\n");
    BARRIER
}

void IceMesh::AssignPrognosticVariables()
{
    // making all prognostic model grid variables
    for (int layer_n = 0; layer_n < mesh_info.num_ice_layers; ++layer_n)
    {
        for (auto& [variable, element]: mesh_info.prog_elems)
            prognostic_data[element][layer_n]->Create(variable, mesh::progDims.at(variable), INMOST::DATA_REAL);
    }

    BARRIER
}

void IceMesh::AssignForcingVariables()
{
    // making all forcing grid variables
    for (auto& [variable, element]: mesh_info.forc_elems)
        forcing_data[element]->Create(variable, mesh::forcDims.at(variable), INMOST::DATA_REAL);

    BARRIER
} 