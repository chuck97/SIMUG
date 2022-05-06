#include "mesh.hpp"

using namespace std;
using namespace INMOST;
using namespace SIMUG::mesh;
using namespace SIMUG::coord;

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
    grid_info[gridElemType::Node] = make_shared<NodeInfo>();
    grid_info[gridElemType::Edge] = make_shared<EdgeInfo>();
    grid_info[gridElemType::Trian] = make_shared<TrianInfo>();

    mesh_timer.Launch();
    ComputeMeshInfo();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();

    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Cumputing mesh information successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    // assign coords for grid elements
    mesh_timer.Launch();
    AssignCoords();
    mesh_timer.Stop();
    duration = mesh_timer.GetMaxTime();
    mesh_timer.Reset();

    if (ice_mesh->GetProcessorRank()==0)
        mesh_log.Log("Assigning coordinates successfully! (" + to_string(duration) + " ms)\n");
    BARRIER

    //AssembleBasisData();
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

    // assign id and no_bnd id for nodes
    grid_info[gridElemType::Node]->id = ice_mesh->CreateTag("id node", INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);
    grid_info[gridElemType::Node]->id_no_bnd = ice_mesh->CreateTag("id node no_bnd", INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);

    int node_id = 0;
    for (int procn = 0; procn < ice_mesh->GetProcessorsNumber(); procn++)
    {
        if (ice_mesh->GetProcessorRank() == procn)
        {
            for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
            {
                if (nodeit->GetStatus() != Element::Ghost)
                    nodeit->Integer(grid_info[gridElemType::Node]->id) = nodeit->GlobalID();

                if ((nodeit->GetStatus() != Element::Ghost) and
                    (nodeit->Integer(grid_info[gridElemType::Node]->is_bnd) == 0))
                {
                    nodeit->Integer(grid_info[gridElemType::Node]->id_no_bnd) = node_id;
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

    // assign id and no_bnd id for edges
    grid_info[gridElemType::Edge]->id = ice_mesh->CreateTag("id edge", INMOST::DATA_INTEGER, INMOST::FACE, INMOST::NONE, 1);
    grid_info[gridElemType::Edge]->id_no_bnd = ice_mesh->CreateTag("id edge no_bnd", INMOST::DATA_INTEGER, INMOST::FACE, INMOST::NONE, 1);

    int edge_id = 0;
    for (int procn = 0; procn < ice_mesh->GetProcessorsNumber(); procn++)
    {
        if (ice_mesh->GetProcessorRank() == procn)
        {
            for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
            {
                if (edgeit->GetStatus() != Element::Ghost)
                    edgeit->Integer(grid_info[gridElemType::Edge]->id) = edgeit->GlobalID();

                if ((edgeit->GetStatus() != Element::Ghost) and
                    (edgeit->Integer(grid_info[gridElemType::Edge]->is_bnd) == 0))
                {
                    edgeit->Integer(grid_info[gridElemType::Edge]->id_no_bnd) = edge_id;
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
    grid_info[gridElemType::Trian]->id = ice_mesh->CreateTag("id trian", INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 1);
    grid_info[gridElemType::Trian]->id_no_bnd = ice_mesh->CreateTag("id trian no_bnd", INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 1);

    int trian_id = 0;
    for (int procn = 0; procn < ice_mesh->GetProcessorsNumber(); procn++)
    {
        if (ice_mesh->GetProcessorRank() == procn)
        {
            for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
            {
                if (trianit->GetStatus() != Element::Ghost)
                trianit->Integer(grid_info[gridElemType::Trian]->id) = trianit->GlobalID();

                if ((trianit->GetStatus() != Element::Ghost) and
                    (trianit->Integer(grid_info[gridElemType::Trian]->is_bnd) == 0))
                {
                    trianit->Integer(grid_info[gridElemType::Trian]->id_no_bnd) = trian_id;
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

#ifdef USE_MPI
    ice_mesh->ExchangeGhost(1, NODE);
#endif

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
            int pid = nodeit->Integer(grid_info[gridElemType::Node]->id);

            if(pid < nodeidmin)
                nodeidmin = pid;
            if((pid + 1) > nodeidmax)
                nodeidmax = pid + 1;

            if (nodeit->Integer(grid_info[gridElemType::Node]->is_bnd) == 0)
            {
                int nobnd_pid = nodeit->Integer(grid_info[gridElemType::Node]->id_no_bnd);
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
            int pid = edgeit->Integer(grid_info[gridElemType::Edge]->id);
            if(pid < edgeidmin)
                edgeidmin = pid;
            if((pid + 1) > edgeidmax)
                edgeidmax = pid + 1;
            
            if (edgeit->Integer(grid_info[gridElemType::Edge]->is_bnd) == 0)
            {
                int nobnd_pid = edgeit->Integer(grid_info[gridElemType::Edge]->id_no_bnd);
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
            int pid = trianit->Integer(grid_info[gridElemType::Trian]->id);
            if(pid < trianidmin)
            {
                trianidmin = pid;
            } 
            if((pid + 1) > trianidmax)
            {
                trianidmax = pid + 1;
            } 

            if (trianit->Integer(grid_info[gridElemType::Trian]->is_bnd) == 0)
            {
                int nobnd_pid = trianit->Integer(grid_info[gridElemType::Trian]->id_no_bnd);
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
    grid_info[gridElemType::Node]->is_bnd = ice_mesh.get()->CreateTag("is node bnd", INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);
    for (size_t i = 0; i < bnd_nodes.size(); ++i)
        bnd_nodes[i].Integer(grid_info[gridElemType::Node]->is_bnd) = 1;

    grid_info[gridElemType::Edge]->is_bnd = ice_mesh.get()->CreateTag("is edge bnd", INMOST::DATA_INTEGER, INMOST::FACE, INMOST::NONE, 1);
    for (size_t i = 0; i < bnd_edges.size(); ++i)
        bnd_edges[i].Integer(grid_info[gridElemType::Edge]->is_bnd) = 1;

    grid_info[gridElemType::Trian]->is_bnd = ice_mesh.get()->CreateTag("is trian bnd", INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 1);
    for (size_t i = 0; i < bnd_trians.size(); ++i)
        bnd_trians[i].Integer(grid_info[gridElemType::Trian]->is_bnd) = 1; 
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

#ifdef USE_MPI
    ice_mesh->ExchangeGhost(1, NODE);
#endif
}

void IceMesh::AssignCoords()
{
    // ## model coords ##

    // nodes
    INMOST::Tag node_coords_tag = ice_mesh->CreateTag("model coords node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    grid_info[gridElemType::Node]->coords[coordType::model] = node_coords_tag;

    for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        nodeit->RealArray(node_coords_tag)[0] = nodeit->Coords()[0];
        nodeit->RealArray(node_coords_tag)[1] = nodeit->Coords()[1];
        nodeit->RealArray(node_coords_tag)[2] = (mesh_info.surface_type == surfType::sphere)?nodeit->Coords()[2]:0.0;
    }

    // edges
    INMOST::Tag edge_coords_tag = ice_mesh->CreateTag("model coords edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    grid_info[gridElemType::Edge]->coords[coordType::model] = edge_coords_tag;

    for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        ElementArray<INMOST::Node> adj_nodes = edgeit->getNodes();
        edgeit->RealArray(edge_coords_tag)[0] = 0.5*(adj_nodes[0].Coords()[0] + adj_nodes[1].Coords()[0]);
        edgeit->RealArray(edge_coords_tag)[1] = 0.5*(adj_nodes[0].Coords()[1] + adj_nodes[1].Coords()[1]);
        edgeit->RealArray(edge_coords_tag)[2] = (mesh_info.surface_type == surfType::sphere)?0.5*(adj_nodes[0].Coords()[2] +
                                                                                                   adj_nodes[1].Coords()[2]):0.0;
    }

    // trians
    INMOST::Tag trian_coords_tag = ice_mesh->CreateTag("model coords trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    grid_info[gridElemType::Trian]->coords[coordType::model] = trian_coords_tag;

    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();
        trianit->RealArray(trian_coords_tag)[0] = (adj_nodes[0].Coords()[0] + adj_nodes[1].Coords()[0] + adj_nodes[2].Coords()[0])/3.0;
        trianit->RealArray(trian_coords_tag)[1] = (adj_nodes[0].Coords()[1] + adj_nodes[1].Coords()[1] + adj_nodes[2].Coords()[1])/3.0;
        trianit->RealArray(trian_coords_tag)[2] = (mesh_info.surface_type == surfType::sphere)?(adj_nodes[0].Coords()[2] +
                                                                                                 adj_nodes[1].Coords()[2] +
                                                                                                 adj_nodes[2].Coords()[2])/3.0:0.0;
    }

    BARRIER

    // ## cartesian coords ##

    // nodes 
    INMOST::Tag node_cart_tag = ice_mesh->CreateTag("cart coords node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    grid_info[gridElemType::Node]->coords[coordType::cart] = node_cart_tag;

    for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        double model_x = nodeit->RealArray(node_coords_tag)[0];
        double model_y = nodeit->RealArray(node_coords_tag)[1];
        double model_z = nodeit->RealArray(node_coords_tag)[2];

        nodeit->RealArray(node_cart_tag)[0] = (mesh_info.surface_type == surfType::basin)?cos(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                          :model_x;
        nodeit->RealArray(node_cart_tag)[1] = (mesh_info.surface_type == surfType::basin)?sin(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                          :model_y;
        nodeit->RealArray(node_cart_tag)[2] = (mesh_info.surface_type == surfType::basin)?sin(model_y*(M_PI/180.0))
                                                                                          :model_z;
    }

    // edges
    INMOST::Tag edge_cart_tag = ice_mesh->CreateTag("cart coords edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    grid_info[gridElemType::Edge]->coords[coordType::cart] = edge_cart_tag;

    for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        double model_x = edgeit->RealArray(edge_coords_tag)[0];
        double model_y = edgeit->RealArray(edge_coords_tag)[1];
        double model_z = edgeit->RealArray(edge_coords_tag)[2];

        edgeit->RealArray(edge_cart_tag)[0] = (mesh_info.surface_type == surfType::basin)?cos(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                          :model_x;
        edgeit->RealArray(edge_cart_tag)[1] = (mesh_info.surface_type == surfType::basin)?sin(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                          :model_y;
        edgeit->RealArray(edge_cart_tag)[2] = (mesh_info.surface_type == surfType::basin)?sin(model_y*(M_PI/180.0))
                                                                                          :model_z;
    }

    // triangles
    INMOST::Tag trian_cart_tag = ice_mesh->CreateTag("cart coords trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    grid_info[gridElemType::Trian]->coords[coordType::cart] = trian_cart_tag;

    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        double model_x = trianit->RealArray(trian_coords_tag)[0];
        double model_y = trianit->RealArray(trian_coords_tag)[1];
        double model_z = trianit->RealArray(trian_coords_tag)[2];

        trianit->RealArray(trian_cart_tag)[0] = (mesh_info.surface_type == surfType::basin)?cos(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                          :model_x;
        trianit->RealArray(trian_cart_tag)[1] = (mesh_info.surface_type == surfType::basin)?sin(model_x*(M_PI/180.0))*cos(model_y*(M_PI/180.0))
                                                                                          :model_y;
        trianit->RealArray(trian_cart_tag)[2] = (mesh_info.surface_type == surfType::basin)?sin(model_y*(M_PI/180.0))
                                                                                          :model_z;
    }

    BARRIER

    // ## geo coords ##

    // nodes 
    INMOST::Tag node_geo_tag = ice_mesh->CreateTag("geo coords node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    grid_info[gridElemType::Node]->coords[coordType::geo] = node_geo_tag;

    for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        if (mesh_info.surface_type == surfType::sphere)
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
        else if (mesh_info.surface_type == surfType::basin)
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

    // edges
    INMOST::Tag edge_geo_tag = ice_mesh->CreateTag("geo coords edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    grid_info[gridElemType::Edge]->coords[coordType::geo] = edge_geo_tag;

    for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        if (mesh_info.surface_type == surfType::sphere)
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
        else if (mesh_info.surface_type == surfType::basin)
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

    // trians
    INMOST::Tag trian_geo_tag = ice_mesh->CreateTag("geo coords trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    grid_info[gridElemType::Trian]->coords[coordType::geo] = trian_geo_tag;

    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        if (mesh_info.surface_type == surfType::sphere)
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
        else if (mesh_info.surface_type == surfType::basin)
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

    BARRIER

#ifdef USE_MPI
    ice_mesh->ExchangeGhost(1, NODE);
#endif
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