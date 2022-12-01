#include "simug.hpp"
#include "inmost.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <memory>
#include <vector>
#include <string>
#include <iostream>

#define PLANE_PATH "../../../../SIMUG_v0/MESHES/pmf/Box_low_res.pmf"

using namespace INMOST;
using namespace SIMUG;
using SIMUG::operator*;
using namespace std;

std::vector<double> analytical_velocity(std::pair<double, double> coords, double time)
{
    return {1.0, 1.0, 0.0};
};

bool test_transition()
{
    IceMesh mesh_plane(PLANE_PATH, mesh::surfType::plane, mesh::gridType::Agrid);

    mesh_plane.GetDataSingle(mesh::gridElemType::Node)->Create("velocity node", 3, INMOST::DATA_REAL);
    mesh_plane.GetDataSingle(mesh::gridElemType::Edge)->Create("velocity edge", 3, INMOST::DATA_REAL);
    mesh_plane.GetDataSingle(mesh::gridElemType::Trian)->Create("velocity trian", 3, INMOST::DATA_REAL);

    INMOST::Tag vel_node_tag = mesh_plane.GetDataSingle(mesh::gridElemType::Node)->Get("velocity node");
    INMOST::Tag vel_edge_tag = mesh_plane.GetDataSingle(mesh::gridElemType::Edge)->Get("velocity edge");
    INMOST::Tag vel_trian_tag = mesh_plane.GetDataSingle(mesh::gridElemType::Trian)->Get("velocity trian");

    // make simple forcing class
    Forcing test_forcing(&mesh_plane);

    // set analytical velocity
    test_forcing.SetAnalytical("velocity node", mesh::gridElemType::Node, analytical_velocity);
    test_forcing.SetAnalytical("velocity edge", mesh::gridElemType::Edge, analytical_velocity);
    test_forcing.SetAnalytical("velocity trian", mesh::gridElemType::Trian, analytical_velocity);

    // update velocity value according to analytical function
    test_forcing.Update("velocity node", mesh::gridElemType::Node, coord::coordType::model, 0.0);

    // node <-> trian
    for (auto trianit = mesh_plane.GetMesh()->BeginCell(); trianit != mesh_plane.GetMesh()->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            // get adj nodes
            ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();
            std::vector<double> vel_trian_coords = mesh_plane.VecTransition({adj_nodes[0]->RealArray(vel_node_tag)[0], adj_nodes[0]->RealArray(vel_node_tag)[1]},
                                                                                     adj_nodes[0], trianit->getCells()[0]); 
            trianit->RealArray(vel_trian_tag)[0] = vel_trian_coords[0];
            trianit->RealArray(vel_trian_tag)[1] = vel_trian_coords[1];
            trianit->RealArray(vel_trian_tag)[2] = 0.0;

            //std::vector<double> vel_node_coords = mesh_plane.VecTransition({trianit->RealArray(vel_trian_tag)[0], trianit->RealArray(vel_trian_tag)[1]},
                                                                                     //trianit->getCells()[0], adj_nodes[0]);
            //adj_nodes[0].RealArray(vel_node_tag)[0] = vel_node_coords[0];
            //adj_nodes[0].RealArray(vel_node_tag)[1] = vel_node_coords[1];
            //adj_nodes[0].RealArray(vel_node_tag)[2] = 0.0;
        }
    }
    mesh_plane.GetDataSingle(mesh::gridElemType::Trian)->Exchange("velocity trian");
    //mesh_plane.GetDataSingle(mesh::gridElemType::Node)->Exchange("velocity node");
/*
    // trian <-> edge
    test_forcing.Update("velocity trian", mesh::gridElemType::Trian, coord::coordType::model, 0.0);
    for (auto edgeit = mesh_plane.GetMesh()->BeginFace(); edgeit != mesh_plane.GetMesh()->EndFace(); ++edgeit)
    {
        if (edgeit->GetStatus() != Element::Ghost)
        {
            // get adj trians
            ElementArray<INMOST::Cell> adj_trians = edgeit->getCells();

            std::vector<double> vel_edge_coords = mesh_plane.VecTransition({adj_trians[0]->RealArray(vel_trian_tag)[0], adj_trians[0]->RealArray(vel_trian_tag)[1]},
                                                                            adj_trians[0], edgeit->getFaces()[0]); 
            edgeit->RealArray(vel_edge_tag)[0] = vel_edge_coords[0];
            edgeit->RealArray(vel_edge_tag)[1] = vel_edge_coords[1];
            edgeit->RealArray(vel_edge_tag)[2] = 0.0;


            std::vector<double> vel_trian_coords = mesh_plane.VecTransition({edgeit->RealArray(vel_edge_tag)[0], edgeit->RealArray(vel_edge_tag)[1]},
                                                                             edgeit->getFaces()[0], adj_trians[0]); 
            adj_trians[0].RealArray(vel_trian_tag)[0] = vel_trian_coords[0];
            adj_trians[0].RealArray(vel_trian_tag)[1] = vel_trian_coords[1];
            adj_trians[0].RealArray(vel_trian_tag)[2] = 0.0;
        }
    }
    mesh_plane.GetDataSingle(mesh::gridElemType::Trian)->Exchange("velocity trian");
    mesh_plane.GetDataSingle(mesh::gridElemType::Edge)->Exchange("velocity edge");


    // node <-> edge
    test_forcing.Update("velocity node", mesh::gridElemType::Node, coord::coordType::model, 0.0);
    for (auto edgeit = mesh_plane.GetMesh()->BeginFace(); edgeit != mesh_plane.GetMesh()->EndFace(); ++edgeit)
    {
        if (edgeit->GetStatus() != Element::Ghost)
        {
            // get adj nodes
            ElementArray<INMOST::Node> adj_nodes = edgeit->getNodes();

            std::vector<double> vel_edge_coords = mesh_plane.VecTransition({adj_nodes[0]->RealArray(vel_node_tag)[0], adj_nodes[0]->RealArray(vel_node_tag)[1]},
                                                                            adj_nodes[0], edgeit->getFaces()[0]); 
            edgeit->RealArray(vel_edge_tag)[0] = vel_edge_coords[0];
            edgeit->RealArray(vel_edge_tag)[1] = vel_edge_coords[1];
            edgeit->RealArray(vel_edge_tag)[2] = 0.0;


            std::vector<double> vel_node_coords = mesh_plane.VecTransition({edgeit->RealArray(vel_edge_tag)[0], edgeit->RealArray(vel_edge_tag)[1]},
                                                                            edgeit->getFaces()[0], adj_nodes[0]); 
            adj_nodes[0].RealArray(vel_node_tag)[0] = vel_node_coords[0];
            adj_nodes[0].RealArray(vel_node_tag)[1] = vel_node_coords[1];
            adj_nodes[0].RealArray(vel_node_tag)[2] = 0.0;
        }
    }
    mesh_plane.GetDataSingle(mesh::gridElemType::Node)->Exchange("velocity node");
    mesh_plane.GetDataSingle(mesh::gridElemType::Edge)->Exchange("velocity edge");

*/
/*
    mesh_plane.GetDataSingle(mesh::gridElemType::Node)->Create("number of adj trians", 1, INMOST::DATA_INTEGER);
    mesh_plane.GetDataSingle(mesh::gridElemType::Trian)->Create("number of adj nodes", 1, INMOST::DATA_INTEGER);
    mesh_plane.GetDataSingle(mesh::gridElemType::Node)->Create("is ghost node", 1, INMOST::DATA_INTEGER);
    mesh_plane.GetDataSingle(mesh::gridElemType::Trian)->Create("is ghost trian", 1, INMOST::DATA_INTEGER);

    INMOST::Tag n_adj_tr = mesh_plane.GetDataSingle(mesh::gridElemType::Node)->Get("number of adj trians");
    INMOST::Tag n_adj_nd = mesh_plane.GetDataSingle(mesh::gridElemType::Trian)->Get("number of adj nodes");
    INMOST::Tag is_ghost_nd = mesh_plane.GetDataSingle(mesh::gridElemType::Node)->Get("is ghost node");
    INMOST::Tag is_ghost_tr = mesh_plane.GetDataSingle(mesh::gridElemType::Trian)->Get("is ghost trian");

    for (auto trianit = mesh_plane.GetMesh()->BeginCell(); trianit != mesh_plane.GetMesh()->EndCell(); ++trianit)
    {
        ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();
        trianit->Integer(n_adj_nd) = adj_nodes.size();
        if (trianit->GetStatus() == Element::Ghost)
            trianit->Integer(is_ghost_tr) = 1;
    }

    for (auto nodeit = mesh_plane.GetMesh()->BeginNode(); nodeit != mesh_plane.GetMesh()->EndNode(); ++nodeit)
    {
        ElementArray<INMOST::Cell> adj_trians = nodeit->getCells();
        nodeit->Integer(n_adj_tr) = adj_trians.size();
        if (nodeit->GetStatus() == Element::Ghost)
            nodeit->Integer(is_ghost_nd) = 1;
    }
*/
    

    mesh_plane.SaveVTU("plane");
    
    return true;
}

int main()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (test_transition())
    {
        if (rank == 0)
            std::cout << "Transition test: OK!\n";
    }
    else
        SIMUG_ERR("Transition test: FAILED!\n");

    BARRIER

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}