#include "defines.hpp"
#include "inmost.h"
#include "model_var.hpp"
#include "mesh_info.hpp"
#include "data.hpp"
#include "mesh.hpp"
#include "vecmath.hpp"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <memory>
#include <vector>
#include <string>
#include <iostream>

#define MESH_PATH "/home/users/spetrov/SIMUG/SIMUG_v0/MESHES/pmf/square8km.pmf"
#define ARCTIC_PATH "/home/users/spetrov/SIMUG/SIMUG_v0/MESHES/pmf/Arctic.pmf"
#define SPHERE_PATH "/home/users/spetrov/SIMUG/SIMUG_v0/MESHES/pmf/Sphere.pmf"

using namespace INMOST;
using namespace SIMUG::mesh;
using SIMUG::operator*;
using namespace std;

bool test_basis()
{

    // calculate angle between cart and geo normal vectors at every point for plane
    {
    IceMesh mesh_plane(MESH_PATH, surfType::plane, gridType::Agrid);
    mesh_plane.GetProgData(gridElemType::Node, 0)->Create("angle node", 1, DATA_REAL);

    INMOST::Tag geo_basis_tag = mesh_plane.GetGridInfo(gridElemType::Node)->geo_basis;
    INMOST::Tag cart_basis_tag = mesh_plane.GetGridInfo(gridElemType::Node)->cart_basis;

    INMOST::Tag angle_node_tag = mesh_plane.GetProgData(gridElemType::Node, 0)->Get("angle node");

    for (auto nodeit = mesh_plane.GetMesh()->BeginNode(); nodeit != mesh_plane.GetMesh()->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            std::vector<double> basis_k_geo = {nodeit->RealArray(geo_basis_tag)[6],
                                               nodeit->RealArray(geo_basis_tag)[7],
                                               nodeit->RealArray(geo_basis_tag)[8]};

            std::vector<double> basis_k_cart = {nodeit->RealArray(cart_basis_tag)[6],
                                                nodeit->RealArray(cart_basis_tag)[7],
                                                nodeit->RealArray(cart_basis_tag)[8]};

            nodeit->Real(angle_node_tag) = SIMUG::angle_vecs(basis_k_geo, basis_k_cart);
        }

    }
    mesh_plane.GetProgData(gridElemType::Node, 0)->Exchange("angle node");
    mesh_plane.SaveVTU("plane");

    }

    // calculate angle between cart and geo normal vectors at every point for sphere
    {
    IceMesh mesh_sphere(SPHERE_PATH, surfType::sphere, gridType::Agrid);
    mesh_sphere.GetProgData(gridElemType::Node, 0)->Create("angle node", 1, DATA_REAL);

    INMOST::Tag geo_basis_tag = mesh_sphere.GetGridInfo(gridElemType::Node)->geo_basis;
    INMOST::Tag cart_basis_tag = mesh_sphere.GetGridInfo(gridElemType::Node)->cart_basis;

    INMOST::Tag angle_node_tag = mesh_sphere.GetProgData(gridElemType::Node, 0)->Get("angle node");

    for (auto nodeit = mesh_sphere.GetMesh()->BeginNode(); nodeit != mesh_sphere.GetMesh()->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            std::vector<double> basis_k_geo = {nodeit->RealArray(geo_basis_tag)[6],
                                               nodeit->RealArray(geo_basis_tag)[7],
                                               nodeit->RealArray(geo_basis_tag)[8]};

            std::vector<double> basis_k_cart = {nodeit->RealArray(cart_basis_tag)[6],
                                                nodeit->RealArray(cart_basis_tag)[7],
                                                nodeit->RealArray(cart_basis_tag)[8]};

            nodeit->Real(angle_node_tag) = SIMUG::angle_vecs(basis_k_geo, basis_k_cart);
        }

    }
    mesh_sphere.GetProgData(gridElemType::Node, 0)->Exchange("angle node");
    mesh_sphere.SaveVTU("sphere");
    }

    // calculate angle between cart and geo normal vectors at every point for Arctic
    {
    IceMesh mesh_arctic(ARCTIC_PATH, surfType::basin, gridType::Agrid);
    mesh_arctic.GetProgData(gridElemType::Node, 0)->Create("angle node", 1, DATA_REAL);
    mesh_arctic.GetProgData(gridElemType::Node, 0)->Create("geo k norm", 1, DATA_REAL);
    mesh_arctic.GetProgData(gridElemType::Node, 0)->Create("cart k norm", 1, DATA_REAL);
    mesh_arctic.GetProgData(gridElemType::Node, 0)->Create("geo k * cart k", 1, DATA_REAL);

    INMOST::Tag geo_basis_tag = mesh_arctic.GetGridInfo(gridElemType::Node)->geo_basis;
    INMOST::Tag cart_basis_tag = mesh_arctic.GetGridInfo(gridElemType::Node)->cart_basis;

    INMOST::Tag angle_node_tag = mesh_arctic.GetProgData(gridElemType::Node, 0)->Get("angle node");
    INMOST::Tag geo_k_tag = mesh_arctic.GetProgData(gridElemType::Node, 0)->Get("geo k norm");
    INMOST::Tag cart_k_tag = mesh_arctic.GetProgData(gridElemType::Node, 0)->Get("cart k norm");
    INMOST::Tag geo_k_mult_cart_k = mesh_arctic.GetProgData(gridElemType::Node, 0)->Get("geo k * cart k");

    for (auto nodeit = mesh_arctic.GetMesh()->BeginNode(); nodeit != mesh_arctic.GetMesh()->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            std::vector<double> basis_k_geo = {nodeit->RealArray(geo_basis_tag)[6],
                                               nodeit->RealArray(geo_basis_tag)[7],
                                               nodeit->RealArray(geo_basis_tag)[8]};

            std::vector<double> basis_k_cart = {nodeit->RealArray(cart_basis_tag)[6],
                                                nodeit->RealArray(cart_basis_tag)[7],
                                                nodeit->RealArray(cart_basis_tag)[8]};

            
            nodeit->Real(angle_node_tag) = SIMUG::angle_vecs(basis_k_geo, basis_k_cart);
            nodeit->Real(geo_k_tag) = SIMUG::L2_norm_vec(basis_k_geo);
            nodeit->Real(cart_k_tag) = SIMUG::L2_norm_vec(basis_k_cart);
            nodeit->Real(geo_k_mult_cart_k) = basis_k_geo*basis_k_cart;
        }

    }
    mesh_arctic.GetProgData(gridElemType::Node, 0)->Exchange("angle node");
    mesh_arctic.GetProgData(gridElemType::Node, 0)->Exchange("geo k norm");
    mesh_arctic.GetProgData(gridElemType::Node, 0)->Exchange("cart k norm");
    mesh_arctic.GetProgData(gridElemType::Node, 0)->Exchange("geo k * cart k");
    mesh_arctic.SaveVTU("arctic");
    }
    
    return true;
}

int main()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (test_basis())
    {
        if (rank == 0)
            std::cout << "Basis test: OK!\n";
    }
    else
        SIMUG_ERR("Basis test: FAILED!\n");

    BARRIER

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}