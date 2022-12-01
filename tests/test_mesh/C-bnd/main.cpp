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

bool test_bnd_selection()
{
    IceMesh imesh(MESH_PATH,
                  mesh::surfType::plane,
                  mesh::gridType::Agrid);

    INMOST::Tag tag_bnd_nodes = imesh.GetMesh()->CreateTag("bnd_nodes", DATA_INTEGER, INMOST::NODE, INMOST::NODE, 1);
    INMOST::Tag tag_bnd_adj_trians = imesh.GetMesh()->CreateTag("adj_tr_for_bnd_edges", DATA_INTEGER, INMOST::CELL, INMOST::CELL, 1);
    INMOST::Tag tag_bnd_trians = imesh.GetMesh()->CreateTag("bnd_trians", DATA_INTEGER, INMOST::CELL, INMOST::CELL, 1);


    for (auto bnodeit = imesh.GetBndNodes().begin(); bnodeit != imesh.GetBndNodes().end(); ++bnodeit)
        bnodeit->Integer(tag_bnd_nodes) = 1.0;

    for (auto btrianit = imesh.GetBndTrians().begin(); btrianit != imesh.GetBndTrians().end(); ++btrianit)
        btrianit->Integer(tag_bnd_trians) = 1.0;
    
    for (auto bedgit = imesh.GetBndEdges().begin(); bedgit != imesh.GetBndEdges().end(); ++bedgit)
    { 
        INMOST::ElementArray<INMOST::Cell> adj_tr =  bedgit->getCells();
        adj_tr[0].Integer(tag_bnd_adj_trians) = 1.0;
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

    if (test_bnd_selection())
    {
        if (rank == 0)
            std::cout << "Mesh bnd elems selection test: OK!\n";
    }
    else
        SIMUG_ERR("Mesh bnd elems selection test: FAILED!\n");

    BARRIER

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}