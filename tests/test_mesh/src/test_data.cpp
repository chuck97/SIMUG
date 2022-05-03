#include "tests.hpp"

using namespace INMOST;
using namespace SIMUG::mesh;
using namespace std;

// Make mesh partition    
void Partition(INMOST::Mesh* ice_mesh)
{
    Partitioner::Initialize(NULL, NULL);

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
};

bool test_data()
{
    
    INMOST::Mesh* ice_mesh = new Mesh();

#ifdef USE_MPI
    ice_mesh->SetCommunicator(INMOST_MPI_COMM_WORLD);
#endif


    if(ice_mesh->isParallelFileFormat(MESH_PATH))
    {
        ice_mesh->Load(MESH_PATH); 

        if (ice_mesh->GetProcessorRank()==0)
        {
            cout << "Parallel realization" << endl;
        }
    }
    else
    {
        if (ice_mesh->GetProcessorRank() == 0)
        {
            ice_mesh->Load(MESH_PATH);
            cout << "Serial realization" << endl;
        }
    }

    if (ice_mesh->GetProcessorRank()==0)
    {
        cout << "Mesh " << MESH_PATH << " loaded successfully " << endl;
        cout << "There are " << ice_mesh->GetProcessorsNumber() << " processes;" << endl;
    }
    BARRIER

#if defined(USE_PARTITIONER)
    Partition(ice_mesh);
#endif

    NodeData node_data(ice_mesh);
    EdgeData edge_data(ice_mesh);
    TrianData triangle_data(ice_mesh);

    node_data.Create(meshVar::mi, meshDim::scalar, INMOST::DATA_REAL);
    node_data.Create("test variable", meshDim::vector, INMOST::DATA_INTEGER);
    node_data.Create("test variable2", 5, INMOST::DATA_INTEGER);
    node_data.Create("test variable3", 10, INMOST::DATA_INTEGER);
    node_data.Create("test variable4", 15, INMOST::DATA_INTEGER);
    node_data.Create("test variable5", 15, INMOST::DATA_INTEGER);
    ice_mesh->DeleteTag(node_data.Get("test variable"), INMOST::NODE);
    node_data.Delete("test variable3");
    node_data.Delete((std::vector<string>){"test variable4", "test variable5"});
    node_data.Delete((std::vector<meshVar>){meshVar::mi});

    edge_data.Create(meshVar::mi, meshDim::scalar, INMOST::DATA_REAL);
    edge_data.Create("test variable", meshDim::vector, INMOST::DATA_INTEGER);
    edge_data.Create("test variable2", 5, INMOST::DATA_INTEGER);
    edge_data.Create("test variable3", 10, INMOST::DATA_INTEGER);
    edge_data.Create("test variable4", 15, INMOST::DATA_INTEGER);
    edge_data.Create("test variable5", 15, INMOST::DATA_INTEGER);
    ice_mesh->DeleteTag(edge_data.Get("test variable"), INMOST::NODE);
    edge_data.Delete("test variable3");
    edge_data.Delete((std::vector<string>){"test variable4", "test variable5"});
    edge_data.Delete((std::vector<meshVar>){meshVar::mi});

    triangle_data.Create(meshVar::mi, meshDim::scalar, INMOST::DATA_REAL);
    triangle_data.Create("test variable", meshDim::vector, INMOST::DATA_INTEGER);
    triangle_data.Create("test variable2", 5, INMOST::DATA_INTEGER);
    triangle_data.Create("test variable3", 10, INMOST::DATA_INTEGER);
    triangle_data.Create("test variable4", 15, INMOST::DATA_INTEGER);
    triangle_data.Create("test variable5", 15, INMOST::DATA_INTEGER);
    ice_mesh->DeleteTag(triangle_data.Get("test variable"), INMOST::NODE);
    triangle_data.Delete("test variable3");
    triangle_data.Delete((std::vector<string>){"test variable4", "test variable5"});
    triangle_data.Delete((std::vector<meshVar>){meshVar::mi});

    // to output faces
    ice_mesh->SetFileOption("VTK_GRID_DIMS", "2");
    ice_mesh->SetFileOption("VTK_OUTPUT_FACES", "1");

#ifdef USE_MPI
    ice_mesh->Save("./res.pvtu");
#else
    ice_mesh->Save("./res.vtu");
#endif

    delete ice_mesh;
    return true;
}