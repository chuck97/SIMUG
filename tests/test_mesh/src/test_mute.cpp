#include "tests.hpp"

using namespace INMOST;
using namespace SIMUG::mesh;
using namespace std;

// Make mesh partition    
static void Partition(INMOST::Mesh* ice_mesh)
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

bool test_mute()
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

#ifdef USE_MPI
    ice_mesh->Save("./res_mute1.pvtu");
#else
    ice_mesh->Save("./res_mute1.vtu");
#endif
    node_data.Mute(meshVar::mi);

#ifdef USE_MPI
    ice_mesh->Save("./res_mute2.pvtu");
#else
    ice_mesh->Save("./res_mute2.vtu");
#endif

    node_data.Unmute(meshVar::mi);

#ifdef USE_MPI
    ice_mesh->Save("./res_mute3.pvtu");
#else
    ice_mesh->Save("./res_mute3.vtu");
#endif

    delete ice_mesh;
    return true;
}