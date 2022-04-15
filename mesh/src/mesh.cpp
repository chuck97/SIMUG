#include "mesh.hpp"

using namespace INMOST;
using namespace std;
using namespace SIMUG;

template <typename RT>
mesh::IceMesh::IceMesh(const MeshConfig& mesh_config_,
                       const ModelConfig& model_config_,
                       std::ostream& os_):
    mesh_config(mesh_config_),
    model_config(model_config_),
    os(os_)
{
    Mesh::Initialize(NULL, NULL);

#if defined(USE_PARTITIONER)
    Partitioner::Initialize(NULL, NULL);
#endif

#if defined(USE_SOLVER)
    Solver::Initialize(NULL, NULL, ""); 
#endif

    bool repartition = false;
    ice_mesh = new Mesh();

#if defined(USE_MPI)
    ice_mesh->SetCommunicator(INMOST_MPI_COMM_WORLD);
    os << "###### Parallel realization ######\n";
    os << "###### " << ice_mesh->GetProcessorsNumber() << " processes ######\n";
#else
    os << "###### Serial realization ######\n";
#endif

    // load mesh
    ice_mesh->Load(model_config.GetVar(model::var::mesh))); 
    repartition = true;

    // Make mesh partition    
    Partition();

    // find boundary verticies and setup ids
    SetBoundaryNodes(true);

    // set id interval
    AssignIntervals();

    // assign model variables
    AssignVariables();

    // assign coords
    AssignCoords();

    // compute mesh info
    ComputeMeshInfo();

    // Calculate vectors transition matricies
    AssembleBasisData(true);
};
/*
void IceMesh::Partition()
{
#if defined(USE_PARTITIONER)
	if (ice_mesh->GetProcessorsNumber() > 1) // need repartition
	{ 
	    Partitioner* p = new Partitioner(ice_mesh);
#ifdef USE_PARTITIONER_PARMETIS
        p->SetMethod(Partitioner::Parmetis, Partitioner::Partition);
#elif USE_PARTITIONER_ZOLTAN
        p->SetMethod(Partitioner::Zoltan, Partitioner::Partition);
#else
        p->SetMethod(Partitioner::INNER_KMEANS, Partitioner::Partition);
#endif
        BARRIER
		p->Evaluate();
		delete p;
		BARRIER

		ice_mesh->Redistribute();
		ice_mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
        ice_mesh->ExchangeGhost(1, NODE); 
        BARRIER

        ice_mesh->AssignGlobalID(CELL|EDGE|FACE|NODE);
        BARRIER
        ice_mesh->ExchangeGhost(1, NODE);
        BARRIER
	}
#endif
};

void IceMesh::PrintPMF(const std::string& filename) const
{
    std::string no_spaces_filename = filename;
    rtrim(no_spaces_filename);
    if(no_spaces_filename.substr(no_spaces_filename.size()-4, 4) != ".pmf")
    {
        INMOST_ICE_ERR("Can't write mesh data to this file. Filename should ended by .pmf")
    }
    BARRIER
	ice_mesh->Save(output_info.GetOutputPvtuDirectory() + no_spaces_filename);
	BARRIER
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Mesh saved to " << no_spaces_filename << ";" << std::endl; 
    }
};

void IceMesh::PrintPVTU(const std::string& filename) const
{
    std::string no_spaces_filename = filename;
    rtrim(no_spaces_filename);
    if(no_spaces_filename.substr(no_spaces_filename.size()-5, 5) != ".pvtu")
    {
        INMOST_ICE_ERR("Can't write mesh to this file. Filename should ended by .pvtu")
    }
    BARRIER
	ice_mesh->Save(output_info.GetOutputPvtuDirectory() + no_spaces_filename);
	BARRIER
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Mesh saved to " << no_spaces_filename << ";" << std::endl; 
    }
};

IceMesh::~IceMesh()
{
#if defined(USE_PARTITIONER)
    Partitioner::Finalize();
#endif

    delete ice_mesh;

#if defined(USE_SOLVER)
    Solver::Finalize();
#endif
    
    Mesh::Finalize();
};

INMOST::Mesh* IceMesh::GetMesh()
{
    return ice_mesh;
};

MeshData& IceMesh::GetData()
{
    return data;
};

std::map<NodeCoordsNotation, INMOST::Tag>& IceMesh::GetCoords()
{
    return node_coords;
};

IdInterval& IceMesh::GetIdIntervalGlobal()
{
    return id_interval_global;
};

IdInterval& IceMesh::GetIdIntervalNoBnd()
{
    return id_interval_no_bnd;
};

LocalBasisData& IceMesh::GetLocalBasisData()
{
    return local_basis_data;
};

void IceMesh::AssignVariables()
{
    // assign node variables
    for (auto item: NodeModelVariableNotationList)
    {
        if (count(NodeVectorNotationList.begin(), NodeVectorNotationList.end(), item) == 0)
        {
            INMOST::Tag node_var = ice_mesh->CreateTag(ModelVariableNotationToName[item],
                                                       DATA_REAL,
                                                       NODE,
                                                       NONE,
                                                       1);
            data.NodeData[item] = node_var;
        }
        else
        {
            INMOST::Tag node_var = ice_mesh->CreateTag(ModelVariableNotationToName[item],
                                                       DATA_REAL,
                                                       NODE,
                                                       NONE,
                                                       3);
            data.NodeData[item] = node_var;
        }

        // Display only needed variables
        if (!output_info.GetDisplayedVariables()[item])
        {
            ice_mesh->SetFileOption("Tag:" + ModelVariableNotationToName[item], "nosave");
        }
    }
    BARRIER

    // assign triangle variables
    for (auto item: TriangleModelVariableNotationList)
    {
        INMOST::Tag trian_var = ice_mesh->CreateTag(ModelVariableNotationToName[item],
                                                    DATA_REAL,
                                                    CELL,
                                                    NONE,
                                                    1);
        data.TriangleData[item] = trian_var;
        
        // Display only needed variables
        if (!output_info.GetDisplayedVariables()[item])
        {
            ice_mesh->SetFileOption("Tag:" + ModelVariableNotationToName[item], "nosave");
        }
    }
    BARRIER
};

void IceMesh::AssignCoords()
{
    if (mesh_params.GetCoordsType() == CoordsType::Cartesian3D)
    {
        // assign model coords
        INMOST::Tag model_coords = ice_mesh->CreateTag(ModelCoordsNotationToName[NodeCoordsNotation::model],
                                                        DATA_REAL,
                                                        NODE,
                                                        NONE,
                                                        3);
        node_coords[NodeCoordsNotation::model] = model_coords;

        // calculate model coords
        for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
            nodeit != ice_mesh->EndNode();
            ++nodeit)
        {
            nodeit->RealArray(node_coords[NodeCoordsNotation::model])[0] = nodeit->Coords()[0]*(M_PI/180.0);
            nodeit->RealArray(node_coords[NodeCoordsNotation::model])[1] = nodeit->Coords()[1]*(M_PI/180.0);
            nodeit->RealArray(node_coords[NodeCoordsNotation::model])[2] = 0.0;
        }
        BARRIER

        // assign geo coords
        INMOST::Tag geo_coords = ice_mesh->CreateTag(ModelCoordsNotationToName[NodeCoordsNotation::geo],
                                                     DATA_REAL,
                                                     NODE,
                                                     NONE,
                                                     3);
        node_coords[NodeCoordsNotation::geo] = geo_coords;
    
        // calculate geo coords
        for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
            nodeit != ice_mesh->EndNode();
            ++nodeit)
        {
            double model_lon = nodeit->RealArray(node_coords[NodeCoordsNotation::model])[0]*(180.0/M_PI);
            double model_lat = nodeit->RealArray(node_coords[NodeCoordsNotation::model])[1]*(180.0/M_PI);

            vector<double> geo_c = from_model_2_geo<double>(model_lon, model_lat);
            nodeit->RealArray(node_coords[NodeCoordsNotation::geo])[0] = geo_c[0]*(M_PI/180.0);
            nodeit->RealArray(node_coords[NodeCoordsNotation::geo])[1] = geo_c[1]*(M_PI/180.0);
            nodeit->RealArray(node_coords[NodeCoordsNotation::geo])[2] = 0.0;
        }
        BARRIER

        // assign Cartesian coords
        INMOST::Tag Cartesian_coords = ice_mesh->CreateTag(ModelCoordsNotationToName[NodeCoordsNotation::Cartesian],
                                                           DATA_REAL,
                                                           NODE,
                                                           NONE,
                                                           3);
        node_coords[NodeCoordsNotation::Cartesian] = Cartesian_coords;

        // calculate Cartesian coords
        double Earth_radius = model_params.GetEarthRadius();
        for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
            nodeit != ice_mesh->EndNode();
            ++nodeit)
        {
            double geo_lon = nodeit->RealArray(node_coords[NodeCoordsNotation::geo])[0];
            double geo_lat = nodeit->RealArray(node_coords[NodeCoordsNotation::geo])[1];

            nodeit->RealArray(node_coords[NodeCoordsNotation::Cartesian])[0] = Earth_radius*cos(geo_lon)*cos(geo_lat);
            nodeit->RealArray(node_coords[NodeCoordsNotation::Cartesian])[1] = Earth_radius*sin(geo_lon)*cos(geo_lat);
            nodeit->RealArray(node_coords[NodeCoordsNotation::Cartesian])[2] = Earth_radius*sin(geo_lat);
        }
        BARRIER
    }
    else if (mesh_params.GetCoordsType() == CoordsType::Cartesian2D)
    {
        // assign model coords
        INMOST::Tag model_coords = ice_mesh->CreateTag(ModelCoordsNotationToName[NodeCoordsNotation::model],
                                                        DATA_REAL,
                                                        NODE,
                                                        NONE,
                                                        3);
        node_coords[NodeCoordsNotation::model] = model_coords;

        // calculate model coords
        for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
            nodeit != ice_mesh->EndNode();
            ++nodeit)
        {
            nodeit->RealArray(node_coords[NodeCoordsNotation::model])[0] = nodeit->Coords()[0];
            nodeit->RealArray(node_coords[NodeCoordsNotation::model])[1] = nodeit->Coords()[1];
            nodeit->RealArray(node_coords[NodeCoordsNotation::model])[2] = 0.0;
        }
        BARRIER

        // assign Cartesian coords
        INMOST::Tag Cartesian_coords = ice_mesh->CreateTag(ModelCoordsNotationToName[NodeCoordsNotation::Cartesian],
                                                           DATA_REAL,
                                                           NODE,
                                                           NONE,
                                                           3);
        node_coords[NodeCoordsNotation::Cartesian] = Cartesian_coords;

        // calculate Cartesian coords
        for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
            nodeit != ice_mesh->EndNode();
            ++nodeit)
        {
            nodeit->RealArray(node_coords[NodeCoordsNotation::Cartesian])[0] = nodeit->Coords()[0];
            nodeit->RealArray(node_coords[NodeCoordsNotation::Cartesian])[1] = nodeit->Coords()[1];
            nodeit->RealArray(node_coords[NodeCoordsNotation::Cartesian])[2] = 0.0;
        }
        BARRIER
    }
    else if (mesh_params.GetCoordsType() == CoordsType::RotatedSpherical2D)
    {
        // assign model coords
        INMOST::Tag model_coords = ice_mesh->CreateTag(ModelCoordsNotationToName[NodeCoordsNotation::model],
                                                        DATA_REAL,
                                                        NODE,
                                                        NONE,
                                                        3);
        node_coords[NodeCoordsNotation::model] = model_coords;

        // calculate model coords
        for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
            nodeit != ice_mesh->EndNode();
            ++nodeit)
        {
            nodeit->RealArray(node_coords[NodeCoordsNotation::model])[0] = nodeit->Coords()[0];
            nodeit->RealArray(node_coords[NodeCoordsNotation::model])[1] = nodeit->Coords()[1];
            nodeit->RealArray(node_coords[NodeCoordsNotation::model])[2] = 0.0;
        }
        BARRIER

        // assign geo coords
        INMOST::Tag geo_coords = ice_mesh->CreateTag(ModelCoordsNotationToName[NodeCoordsNotation::geo],
                                                     DATA_REAL,
                                                     NODE,
                                                     NONE,
                                                     3);
        node_coords[NodeCoordsNotation::geo] = geo_coords;
    
        // calculate geo coords
        Euler_rotation_info<double> rotation(ALPHA_DEF, BETA_DEF, GAMMA_DEF); 
        for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
            nodeit != ice_mesh->EndNode();
            ++nodeit)
        {
            double model_lon = nodeit->RealArray(node_coords[NodeCoordsNotation::model])[0] = nodeit->Coords()[0];
            double model_lat = nodeit->RealArray(node_coords[NodeCoordsNotation::model])[1] = nodeit->Coords()[1];

            Spherical_Coords<double> model_forw(model_lon, model_lat);
            Spherical_Coords<double> geo_forw = Rotate_Spherical<double>(model_forw, rotation.Get_FORWARD());
   
            double lon = geo_forw.Get_x();
            double lat = geo_forw.Get_y();

            nodeit->RealArray(node_coords[NodeCoordsNotation::geo])[0] = lon*M_PI/180.0;
            nodeit->RealArray(node_coords[NodeCoordsNotation::geo])[1] = lat*M_PI/180.0;
            nodeit->RealArray(node_coords[NodeCoordsNotation::geo])[2] = 0.0;
        }
        BARRIER

        // assign Cartesian coords
        INMOST::Tag Cartesian_coords = ice_mesh->CreateTag(ModelCoordsNotationToName[NodeCoordsNotation::Cartesian],
                                                           DATA_REAL,
                                                           NODE,
                                                           NONE,
                                                           3);
        node_coords[NodeCoordsNotation::Cartesian] = Cartesian_coords;

        // calculate Cartesian coords (with rotated XYZ system)
        double Earth_radius = model_params.GetEarthRadius();
        for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
            nodeit != ice_mesh->EndNode();
            ++nodeit)
        {
            double model_lon = nodeit->RealArray(node_coords[NodeCoordsNotation::model])[0]*(M_PI/180.0);
            double model_lat = nodeit->RealArray(node_coords[NodeCoordsNotation::model])[1]*(M_PI/180.0);

            nodeit->RealArray(node_coords[NodeCoordsNotation::Cartesian])[0] = Earth_radius*cos(model_lon)*cos(model_lat);
            nodeit->RealArray(node_coords[NodeCoordsNotation::Cartesian])[1] = Earth_radius*sin(model_lon)*cos(model_lat);
            nodeit->RealArray(node_coords[NodeCoordsNotation::Cartesian])[2] = Earth_radius*sin(model_lat);
        }
        BARRIER

        // calculate TOPAZ stereographic coords
        INMOST::Tag topaz_coords = ice_mesh->CreateTag(ModelCoordsNotationToName[NodeCoordsNotation::topaz_stereographic],
                                                        DATA_REAL,
                                                        NODE,
                                                        NONE,
                                                        3);

        node_coords[NodeCoordsNotation::topaz_stereographic] = topaz_coords;

        ice_mesh->SetFileOption("Tag:"+ModelCoordsNotationToName[NodeCoordsNotation::topaz_stereographic],
                                "nosave");

        for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
            nodeit != ice_mesh->EndNode();
            ++nodeit)
        {
            double geo_lon_deg = nodeit->RealArray(node_coords[NodeCoordsNotation::geo])[0]*(180.0/M_PI);
            double geo_lat_deg = nodeit->RealArray(node_coords[NodeCoordsNotation::geo])[1]*(180.0/M_PI);

            vector<double> topaz_c = from_geo_2_topaz(geo_lon_deg, geo_lat_deg);

            nodeit->RealArray(node_coords[NodeCoordsNotation::topaz_stereographic])[0] = topaz_c[0];
            nodeit->RealArray(node_coords[NodeCoordsNotation::topaz_stereographic])[1] = topaz_c[1];
            nodeit->RealArray(node_coords[NodeCoordsNotation::topaz_stereographic])[2] = 0.0;
        }
        BARRIER
    }
    else
    {
        INMOST_ICE_ERR("Unknown type of coordinates");
    }
};

void IceMesh::AssignIdIntervals()
{

    // calculate global interval
    int idmin = std::numeric_limits<int>::max();
    int idmax = std::numeric_limits<int>::min();

    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
        nodeit != ice_mesh->EndNode();
        ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            int pid = nodeit->Integer(data.NodeIdGlobal);
            if(pid < idmin)
            {
                idmin = pid;
            } 
            if((pid + 1) > idmax)
            {
                idmax = pid + 1;
            } 
        }
    }
    BARRIER
    id_interval_global.IdMin = idmin;
    id_interval_global.IdMax = idmax;
    BARRIER

    // calculate no bnd interval
    idmin = std::numeric_limits<int>::max();
    idmax = std::numeric_limits<int>::min();

    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
        nodeit != ice_mesh->EndNode();
        ++nodeit)
    {
        if((nodeit->GetStatus() != Element::Ghost) and
           (nodeit->Integer(data.IsNodeBnd) == 0))
        {
            int pid = nodeit->Integer(data.NodeIdNoBnd);
            if(pid < idmin)
            {
                idmin = pid;
            } 
            if((pid + 1) > idmax)
            {
                idmax = pid + 1;
            } 
        }
    }
    BARRIER
    id_interval_no_bnd.IdMin = idmin;
    id_interval_no_bnd.IdMax = idmax;
    BARRIER
};

void IceMesh::SetBoundaryNodes(bool display)
{
    data.IsNodeBnd = ice_mesh->CreateTag("is node bnd", DATA_INTEGER, NODE, NONE, 1);

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
    
};

void IceMesh::CalculateMeshInfo()
{
    // calculate number of nodes
    int n_nodes_local = 0;
    int n_nodes_global = 0;
    int n_b_nodes_local = 0;
    int n_b_nodes_global = 0;

    for (Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
                nodeit != ice_mesh->EndNode();
                ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            if (nodeit->Integer(data.IsNodeBnd) == 1)
            {
                ++n_b_nodes_local;
            }
            ++n_nodes_local;
        }
    }

#if defined(USE_MPI)
    MPI_Allreduce(&n_b_nodes_local, &n_b_nodes_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
#endif

#if defined(USE_MPI)
    MPI_Allreduce(&n_nodes_local, &n_nodes_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
#endif
    int n_nb_nodes_global = n_nodes_global - n_b_nodes_global;

    BARRIER

    // calculate number of triangles and tr size
    int n_triangles_local = 0;
    int n_triangles_global = 0;

    double min_tr_size_local = std::numeric_limits<double>::max();
    double min_tr_size_global = std::numeric_limits<double>::max();
    double max_tr_size_local = std::numeric_limits<double>::min();
    double max_tr_size_global = std::numeric_limits<double>::min();

    for(Mesh::iteratorCell trianit = ice_mesh->BeginCell();
                trianit != ice_mesh->EndCell();
                ++trianit) 
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            ++n_triangles_local;
             // calculate edge sizes
            ElementArray<Node> trian_nodes = trianit->getNodes();
            vector<double> node_coords_0 = 
            {trian_nodes[0].RealArray(node_coords[NodeCoordsNotation::Cartesian])[0],
            trian_nodes[0].RealArray(node_coords[NodeCoordsNotation::Cartesian])[1]};
            vector<double> node_coords_1 =
            {trian_nodes[1].RealArray(node_coords[NodeCoordsNotation::Cartesian])[0],
             trian_nodes[1].RealArray(node_coords[NodeCoordsNotation::Cartesian])[1]};
            vector<double> node_coords_2 =
            {trian_nodes[2].RealArray(node_coords[NodeCoordsNotation::Cartesian])[0],
             trian_nodes[2].RealArray(node_coords[NodeCoordsNotation::Cartesian])[1]};
            
            vector<vector<double>> node_coords = {node_coords_0, node_coords_1, node_coords_2};
            

            if (mesh_params.GetCoordsType() == CoordsType::Cartesian2D)
            {
                double a0 = sqrt((node_coords[0][0] - node_coords[1][0])*(node_coords[0][0] - node_coords[1][0]) +
                                 (node_coords[0][1] - node_coords[1][1])*(node_coords[0][1] - node_coords[1][1]));
                double a1 = sqrt((node_coords[0][0] - node_coords[2][0])*(node_coords[0][0] - node_coords[2][0]) +
                                 (node_coords[0][1] - node_coords[2][1])*(node_coords[0][1] - node_coords[2][1]));
                double a2 = sqrt((node_coords[1][0] - node_coords[2][0])*(node_coords[1][0] - node_coords[2][0]) +
                                 (node_coords[1][1] - node_coords[2][1])*(node_coords[1][1] - node_coords[2][1]));
                double min_side = min(min(a0, a1), a2);
                double max_side = max(max(a0, a1), a2);

                if (min_side < min_tr_size_local)
                {
                    min_tr_size_local = min_side;
                }

                if (max_side > max_tr_size_local)
                {
                    max_tr_size_local = max_side;
                }
            }
        }
    }

#if defined(USE_MPI)
    MPI_Allreduce(&n_triangles_local, &n_triangles_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
#endif

#if defined(USE_MPI)
    MPI_Allreduce(&min_tr_size_local, &min_tr_size_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

#if defined(USE_MPI)
    MPI_Allreduce(&max_tr_size_local, &max_tr_size_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    BARRIER
    mesh_info = {n_triangles_global,
                 n_nodes_global,
                 n_b_nodes_global,
                 n_nb_nodes_global,
                 min_tr_size_global,
                 max_tr_size_global};
};

MeshInfo& IceMesh::GetMeshInfo()
{
    return mesh_info;
};

void IceMesh::AssembleLocalBasisData(bool display)
{
    // calculate local triangle basis and vectors transition matricies
    local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle = ice_mesh->CreateTag("vec_transition_matricies_n_to_t", DATA_REAL, CELL, NONE, 12);
    local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal = ice_mesh->CreateTag("vec_transition_matricies_t_to_n", DATA_REAL, CELL, NONE, 12);
    local_basis_data.triangle_data.LocalNodeCoords = ice_mesh->CreateTag("local_node_coords", DATA_REAL, CELL, NONE, 6);
    local_basis_data.nodal_data.LocalNodeCoords = ice_mesh->CreateTag("local_node_coords_in_nodal_basis", DATA_REAL, CELL, NONE, 18);
    INMOST::Tag tg = ice_mesh->CreateTag("det", DATA_REAL, CELL, NONE, 3);

    if (!display)
    {
        ice_mesh->SetFileOption("Tag:vec_transition_matricies_n_to_t", "nosave");
        ice_mesh->SetFileOption("Tag:vec_transition_matricies_t_to_n", "nosave");
        ice_mesh->SetFileOption("Tag:local_node_coords", "nosave");
        ice_mesh->SetFileOption("Tag:local_node_coords_in_nodal_basis", "nosave");
    }

    for(Mesh::iteratorCell trianit = ice_mesh->BeginCell();
            trianit != ice_mesh->EndCell();
            ++trianit) 
    {
        // calculate local basis on triangle
        vector<double> basis_i(3), basis_j(3), basis_k(3);
        ElementArray<Node> local_nodes = trianit->getNodes();
        
        vector<double> v0_coords(3), v1_coords(3), v2_coords(3);
        v0_coords[0] = local_nodes[0]->RealArray(node_coords[NodeCoordsNotation::Cartesian])[0];
        v0_coords[1] = local_nodes[0]->RealArray(node_coords[NodeCoordsNotation::Cartesian])[1];
        v0_coords[2] = local_nodes[0]->RealArray(node_coords[NodeCoordsNotation::Cartesian])[2];
        v1_coords[0] = local_nodes[1]->RealArray(node_coords[NodeCoordsNotation::Cartesian])[0];
        v1_coords[1] = local_nodes[1]->RealArray(node_coords[NodeCoordsNotation::Cartesian])[1];
        v1_coords[2] = local_nodes[1]->RealArray(node_coords[NodeCoordsNotation::Cartesian])[2];
        v2_coords[0] = local_nodes[2]->RealArray(node_coords[NodeCoordsNotation::Cartesian])[0];
        v2_coords[1] = local_nodes[2]->RealArray(node_coords[NodeCoordsNotation::Cartesian])[1];
        v2_coords[2] = local_nodes[2]->RealArray(node_coords[NodeCoordsNotation::Cartesian])[2];
        
        std::vector<double> center_coords(3);
        center_coords[0] = (v0_coords[0] + v1_coords[0] + v2_coords[0])/3.0;
        center_coords[1] = (v0_coords[1] + v1_coords[1] + v2_coords[1])/3.0;
        center_coords[2] = (v0_coords[2] + v1_coords[2] + v2_coords[2])/3.0;
        
        // calculating normal on triangle
        std::vector<double> unit_normal_vec = unit_normal(v0_coords, v1_coords, v2_coords);
        basis_k = unit_normal_vec;
        
        // calculating basis i vector
        std::vector<double> c_0_coords = v0_coords - center_coords;
        double c_0_mod = vec_mod(c_0_coords);
        basis_i = {c_0_coords[0]/c_0_mod, c_0_coords[1]/c_0_mod, c_0_coords[2]/c_0_mod}; 
        
        // calculating basis j vector
        basis_j = vec_product(basis_k, basis_i); 
        
        // local node basis vectors
        std::vector<double> basis_i0(3), basis_j0(3), basis_k0(3);
        std::vector<double> basis_i1(3), basis_j1(3), basis_k1(3);
        std::vector<double> basis_i2(3), basis_j2(3), basis_k2(3);

        if (mesh_params.GetCoordsType() == CoordsType::Cartesian3D)
        {
            // calculating local spherical basis
            double lon0 = local_nodes[0]->RealArray(node_coords[NodeCoordsNotation::geo])[0];
            double lat0 = local_nodes[0]->RealArray(node_coords[NodeCoordsNotation::geo])[1];
            double lon1 = local_nodes[1]->RealArray(node_coords[NodeCoordsNotation::geo])[0];
            double lat1 = local_nodes[1]->RealArray(node_coords[NodeCoordsNotation::geo])[1];
            double lon2 = local_nodes[2]->RealArray(node_coords[NodeCoordsNotation::geo])[0];
            double lat2 = local_nodes[2]->RealArray(node_coords[NodeCoordsNotation::geo])[1]; 
            
            // first node
            basis_i0 = {-sin(lon0), cos(lon0), 0.0};
            basis_j0 = {-sin(lat0)*cos(lon0), -sin(lat0)*sin(lon0), cos(lat0)};
            basis_k0 = {cos(lat0)*cos(lon0), cos(lat0)*sin(lon0), sin(lat0)};
            
            // second node
            basis_i1 = {-sin(lon1), cos(lon1), 0.0};
            basis_j1 = {-sin(lat1)*cos(lon1), -sin(lat1)*sin(lon1), cos(lat1)};
            basis_k1 = {cos(lat1)*cos(lon1), cos(lat1)*sin(lon1), sin(lat1)};
            
            // third node
            basis_i2 = {-sin(lon2), cos(lon2), 0.0};
            basis_j2 = {-sin(lat2)*cos(lon2), -sin(lat2)*sin(lon2), cos(lat2)};
            basis_k2 = {cos(lat2)*cos(lon2), cos(lat2)*sin(lon2), sin(lat2)};
        }
        else if (mesh_params.GetCoordsType() == CoordsType::Cartesian2D)
        {
            // first node
            basis_i0 = {1.0, 0.0, 0.0};
            basis_j0 = {0.0, 1.0, 0.0};
            basis_k0 = {0.0, 0.0, 1.0};

            // second node
            basis_i1 = {1.0, 0.0, 0.0};
            basis_j1 = {0.0, 1.0, 0.0};
            basis_k1 = {0.0, 0.0, 1.0};

            // third node
            basis_i2 = {1.0, 0.0, 0.0};
            basis_j2 = {0.0, 1.0, 0.0};
            basis_k2 = {0.0, 0.0, 1.0};
        }
        else if (mesh_params.GetCoordsType() == CoordsType::RotatedSpherical2D)
        {
            // calculating local spherical basis
            double lon0 = local_nodes[0]->RealArray(node_coords[NodeCoordsNotation::model])[0]*(M_PI/180.0);
            double lat0 = local_nodes[0]->RealArray(node_coords[NodeCoordsNotation::model])[1]*(M_PI/180.0);
            double lon1 = local_nodes[1]->RealArray(node_coords[NodeCoordsNotation::model])[0]*(M_PI/180.0);
            double lat1 = local_nodes[1]->RealArray(node_coords[NodeCoordsNotation::model])[1]*(M_PI/180.0);
            double lon2 = local_nodes[2]->RealArray(node_coords[NodeCoordsNotation::model])[0]*(M_PI/180.0);
            double lat2 = local_nodes[2]->RealArray(node_coords[NodeCoordsNotation::model])[1]*(M_PI/180.0); 
            
            // first node
            basis_i0 = {-sin(lon0), cos(lon0), 0.0};
            basis_j0 = {-sin(lat0)*cos(lon0), -sin(lat0)*sin(lon0), cos(lat0)};
            basis_k0 = {cos(lat0)*cos(lon0), cos(lat0)*sin(lon0), sin(lat0)};

            // second node
            basis_i1 = {-sin(lon1), cos(lon1), 0.0};
            basis_j1 = {-sin(lat1)*cos(lon1), -sin(lat1)*sin(lon1), cos(lat1)};
            basis_k1 = {cos(lat1)*cos(lon1), cos(lat1)*sin(lon1), sin(lat1)};
            
            // third node
            basis_i2 = {-sin(lon2), cos(lon2), 0.0};
            basis_j2 = {-sin(lat2)*cos(lon2), -sin(lat2)*sin(lon2), cos(lat2)};
            basis_k2 = {cos(lat2)*cos(lon2), cos(lat2)*sin(lon2), sin(lat2)};
        }
        else
        {
            INMOST_ICE_ERR("Unknown type of coordinates");
        }

        // first transitiion operator

        double a0 = VectorsScalarProduct(basis_i, basis_i0);
        double b0 = VectorsScalarProduct(basis_i, basis_j0);
        double c0 = VectorsScalarProduct(basis_j, basis_i0);
        double d0 = VectorsScalarProduct(basis_j, basis_j0);
        double det0 = a0*d0 - b0*c0;
        
        trianit->RealArray(tg)[0] =  det0;

        // forward
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[0] = a0;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[1] = b0;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[2] = c0;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[3] = d0;

        // reversed
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[0] =  d0/det0; 
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[1] = -b0/det0;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[2] = -c0/det0;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[3] =  a0/det0;

        // second transitiion operator

        double a1 = VectorsScalarProduct(basis_i, basis_i1);
        double b1 = VectorsScalarProduct(basis_i, basis_j1);
        double c1 = VectorsScalarProduct(basis_j, basis_i1);
        double d1 = VectorsScalarProduct(basis_j, basis_j1);
        double det1 = a1*d1 - b1*c1;

        trianit->RealArray(tg)[1] =  det1;

        // forward
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[4] = a1;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[5] = b1;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[6] = c1;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[7] = d1;

        // reversed
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[4] =  d1/det1; 
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[5] = -b1/det1;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[6] = -c1/det1;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[7] =  a1/det1;

        // third transitiion operator

        double a2 = VectorsScalarProduct(basis_i, basis_i2);
        double b2 = VectorsScalarProduct(basis_i, basis_j2);
        double c2 = VectorsScalarProduct(basis_j, basis_i2);
        double d2 = VectorsScalarProduct(basis_j, basis_j2);
        double det2 = a2*d2 - b2*c2;

        trianit->RealArray(tg)[2] =  det2;

        // forward
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[8]  = a2;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[9]  = b2;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[10] = c2;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromNodalToTriangle)[11] = d2;

        // reversed
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[8]  =  d2/det2; 
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[9]  = -b2/det2;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[10] = -c2/det2;
        trianit->RealArray(local_basis_data.triangle_data.VecTransMatriciesFromTriangleToNodal)[11] =  a2/det2;

        // find coords of v0 in local triangle basis
        std::vector<double> v0_vec = v0_coords - center_coords;
        std::vector<double> v0_local_coords = Solve3x3({basis_i, basis_j, basis_k}, v0_vec);
        trianit->RealArray(local_basis_data.triangle_data.LocalNodeCoords)[0] = v0_local_coords[0];
        trianit->RealArray(local_basis_data.triangle_data.LocalNodeCoords)[1] = v0_local_coords[1]; 

        // find coords of v1 in local triangle basis
        std::vector<double> v1_vec = v1_coords - center_coords;
        std::vector<double> v1_local_coords = Solve3x3({basis_i, basis_j, basis_k}, v1_vec);
        trianit->RealArray(local_basis_data.triangle_data.LocalNodeCoords)[2] = v1_local_coords[0];
        trianit->RealArray(local_basis_data.triangle_data.LocalNodeCoords)[3] = v1_local_coords[1];

        // find coords of v2 in local triangle basis
        std::vector<double> v2_vec = v2_coords - center_coords;
        std::vector<double> v2_local_coords = Solve3x3({basis_i, basis_j, basis_k}, v2_vec);
        trianit->RealArray(local_basis_data.triangle_data.LocalNodeCoords)[4] = v2_local_coords[0];
        trianit->RealArray(local_basis_data.triangle_data.LocalNodeCoords)[5] = v2_local_coords[1];

        // find node coords in nodal basis

        // for first nodal basis
        std::vector<double> v00_vec = v0_coords - v0_coords;
        std::vector<double> v01_vec = v1_coords - v0_coords;
        std::vector<double> v02_vec = v2_coords - v0_coords;

        std::vector<double> v0_local_coords_0 = Solve3x3({basis_i0, basis_j0, basis_k0}, v00_vec);
        std::vector<double> v1_local_coords_0 = Solve3x3({basis_i0, basis_j0, basis_k0}, v01_vec);
        std::vector<double> v2_local_coords_0 = Solve3x3({basis_i0, basis_j0, basis_k0}, v02_vec);

        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[0] = v0_local_coords_0[0];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[1] = v0_local_coords_0[1];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[2] = v1_local_coords_0[0];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[3] = v1_local_coords_0[1];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[4] = v2_local_coords_0[0];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[5] = v2_local_coords_0[1];

        // for second nodal basis
        std::vector<double> v10_vec = v0_coords - v1_coords;
        std::vector<double> v11_vec = v1_coords - v1_coords;
        std::vector<double> v12_vec = v2_coords - v1_coords;

        std::vector<double> v0_local_coords_1 = Solve3x3({basis_i1, basis_j1, basis_k1}, v10_vec);
        std::vector<double> v1_local_coords_1 = Solve3x3({basis_i1, basis_j1, basis_k1}, v11_vec);
        std::vector<double> v2_local_coords_1 = Solve3x3({basis_i1, basis_j1, basis_k1}, v12_vec);
        
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[6]  = v0_local_coords_1[0];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[7]  = v0_local_coords_1[1];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[8]  = v1_local_coords_1[0];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[9]  = v1_local_coords_1[1];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[10] = v2_local_coords_1[0];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[11] = v2_local_coords_1[1];

        // for third nodal basis
        std::vector<double> v20_vec = v0_coords - v2_coords;
        std::vector<double> v21_vec = v1_coords - v2_coords;
        std::vector<double> v22_vec = v2_coords - v2_coords;

        std::vector<double> v0_local_coords_2 = Solve3x3({basis_i2, basis_j2, basis_k2}, v20_vec);
        std::vector<double> v1_local_coords_2 = Solve3x3({basis_i2, basis_j2, basis_k2}, v21_vec);
        std::vector<double> v2_local_coords_2 = Solve3x3({basis_i2, basis_j2, basis_k2}, v22_vec);
        
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[12] = v0_local_coords_2[0];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[13] = v0_local_coords_2[1];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[14] = v1_local_coords_2[0];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[15] = v1_local_coords_2[1];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[16] = v2_local_coords_2[0];
        trianit->RealArray(local_basis_data.nodal_data.LocalNodeCoords)[17] = v2_local_coords_2[1];
    }
    BARRIER
};
*/

template class MeshData<float>;
template class MeshData<double>;

template class NodeData<float>;
template class NodeData<double>;

template class EdgeData<float>;
template class EdgeData<double>;

template class ElemData<float>;
template class ElemData<double>;

template class IceMesh<float>;
template class IceMesh<double>;