#include "mesh.hpp"

using namespace std;
using namespace INMOST;
using namespace SIMUG::mesh;
using namespace SIMUG::coord;

void IceMesh::AssembleGeoElementBasis()
{
    // nodes
    INMOST::Tag node_geo_basis_x_tag = ice_mesh->CreateTag("geo basis x node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    INMOST::Tag node_geo_basis_y_tag = ice_mesh->CreateTag("geo basis y node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    INMOST::Tag node_geo_basis_z_tag = ice_mesh->CreateTag("geo basis z node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);

    grid_info[gridElemType::Node]->geo_basis = {node_geo_basis_x_tag, node_geo_basis_y_tag, node_geo_basis_z_tag};

    for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        // node geo coords
        double lon = nodeit->RealArray(grid_info[gridElemType::Node]->coords[coordType::geo])[0];
        double lat = nodeit->RealArray(grid_info[gridElemType::Node]->coords[coordType::geo])[1];

        // first basis vector
        nodeit->RealArray(node_geo_basis_x_tag)[0] = (mesh_info.surface_type == surfType::plane)?1.0:-sin(lon);
        nodeit->RealArray(node_geo_basis_x_tag)[1] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lon);
        nodeit->RealArray(node_geo_basis_x_tag)[2] = (mesh_info.surface_type == surfType::plane)?0.0:0.0;

        // second basis vector
        nodeit->RealArray(node_geo_basis_y_tag)[0] = (mesh_info.surface_type == surfType::plane)?0.0:-sin(lat)*cos(lon);
        nodeit->RealArray(node_geo_basis_y_tag)[1] = (mesh_info.surface_type == surfType::plane)?1.0:-sin(lat)*sin(lon);
        nodeit->RealArray(node_geo_basis_y_tag)[2] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lat);

        // third basis vector
        nodeit->RealArray(node_geo_basis_z_tag)[0] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lat)*cos(lon);
        nodeit->RealArray(node_geo_basis_z_tag)[1] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lat)*sin(lon);
        nodeit->RealArray(node_geo_basis_z_tag)[2] = (mesh_info.surface_type == surfType::plane)?1.0:sin(lat);
    }

    // edges
    INMOST::Tag edge_geo_basis_x_tag = ice_mesh->CreateTag("geo basis x edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    INMOST::Tag edge_geo_basis_y_tag = ice_mesh->CreateTag("geo basis y edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    INMOST::Tag edge_geo_basis_z_tag = ice_mesh->CreateTag("geo basis z edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);

    grid_info[gridElemType::Edge]->geo_basis = {edge_geo_basis_x_tag, edge_geo_basis_y_tag, edge_geo_basis_z_tag};

    for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        // edge geo coords
        double lon = edgeit->RealArray(grid_info[gridElemType::Edge]->coords[coordType::geo])[0];
        double lat = edgeit->RealArray(grid_info[gridElemType::Edge]->coords[coordType::geo])[1];

        // first basis vector
        edgeit->RealArray(edge_geo_basis_x_tag)[0] = (mesh_info.surface_type == surfType::plane)?1.0:-sin(lon);
        edgeit->RealArray(edge_geo_basis_x_tag)[1] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lon);
        edgeit->RealArray(edge_geo_basis_x_tag)[2] = (mesh_info.surface_type == surfType::plane)?0.0:0.0;

        // second basis vector
        edgeit->RealArray(edge_geo_basis_y_tag)[0] = (mesh_info.surface_type == surfType::plane)?0.0:-sin(lat)*cos(lon);
        edgeit->RealArray(edge_geo_basis_y_tag)[1] = (mesh_info.surface_type == surfType::plane)?1.0:-sin(lat)*sin(lon);
        edgeit->RealArray(edge_geo_basis_y_tag)[2] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lat);

        // third basis vector
        edgeit->RealArray(edge_geo_basis_z_tag)[0] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lat)*cos(lon);
        edgeit->RealArray(edge_geo_basis_z_tag)[1] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lat)*sin(lon);
        edgeit->RealArray(edge_geo_basis_z_tag)[2] = (mesh_info.surface_type == surfType::plane)?1.0:sin(lat);
    }

    // trians
    INMOST::Tag trian_geo_basis_x_tag = ice_mesh->CreateTag("geo basis x trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    INMOST::Tag trian_geo_basis_y_tag = ice_mesh->CreateTag("geo basis y trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    INMOST::Tag trian_geo_basis_z_tag = ice_mesh->CreateTag("geo basis z trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);

    grid_info[gridElemType::Trian]->geo_basis = {trian_geo_basis_x_tag, trian_geo_basis_y_tag, trian_geo_basis_z_tag};

    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        // trian geo coords
        double lon = trianit->RealArray(grid_info[gridElemType::Trian]->coords[coordType::geo])[0];
        double lat = trianit->RealArray(grid_info[gridElemType::Trian]->coords[coordType::geo])[1];

        // first basis vector
        trianit->RealArray(trian_geo_basis_x_tag)[0] = (mesh_info.surface_type == surfType::plane)?1.0:-sin(lon);
        trianit->RealArray(trian_geo_basis_x_tag)[1] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lon);
        trianit->RealArray(trian_geo_basis_x_tag)[2] = (mesh_info.surface_type == surfType::plane)?0.0:0.0;

        // second basis vector
        trianit->RealArray(trian_geo_basis_y_tag)[0] = (mesh_info.surface_type == surfType::plane)?0.0:-sin(lat)*cos(lon);
        trianit->RealArray(trian_geo_basis_y_tag)[1] = (mesh_info.surface_type == surfType::plane)?1.0:-sin(lat)*sin(lon);
        trianit->RealArray(trian_geo_basis_y_tag)[2] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lat);

        // third basis vector
        trianit->RealArray(trian_geo_basis_z_tag)[0] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lat)*cos(lon);
        trianit->RealArray(trian_geo_basis_z_tag)[1] = (mesh_info.surface_type == surfType::plane)?0.0:cos(lat)*sin(lon);
        trianit->RealArray(trian_geo_basis_z_tag)[2] = (mesh_info.surface_type == surfType::plane)?1.0:sin(lat);
    }
    BARRIER
}

void IceMesh::AssembleCartesianElementBasis()
{
    // ## triangles ##
    INMOST::Tag trian_cart_basis_x_tag = ice_mesh->CreateTag("cart basis x trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    INMOST::Tag trian_cart_basis_y_tag = ice_mesh->CreateTag("cart basis y trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    INMOST::Tag trian_cart_basis_z_tag = ice_mesh->CreateTag("cart basis z trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);

    grid_info[gridElemType::Trian]->cart_basis = {trian_cart_basis_x_tag, trian_cart_basis_y_tag, trian_cart_basis_z_tag};

    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        // gain adj nodes to triangle
        ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

        // gain cartesian coords of nodes
        double node0_x = adj_nodes[0].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[0];
        double node0_y = adj_nodes[0].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[1];
        double node0_z = adj_nodes[0].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[2];

        double node1_x = adj_nodes[1].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[0];
        double node1_y = adj_nodes[1].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[1];
        double node1_z = adj_nodes[1].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[2];

        double node2_x = adj_nodes[2].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[0];
        double node2_y = adj_nodes[2].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[1];
        double node2_z = adj_nodes[2].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[2];

        // calculate centroid coordinates
        double centr_x = (node0_x + node1_x + node2_x)/3.0;
        double centr_y = (node0_y + node1_y + node2_y)/3.0;
        double centr_z = (node0_z + node1_z + node2_z)/3.0;

        // calculate unit normal to triangle
        std::vector<double> unit_normal_vec = (mesh_info.surface_type == surfType::basin)? unit_normal<double>({node0_x, node0_y, node0_z},
                                                                                                               {node2_x, node2_y, node2_z},
                                                                                                               {node1_x, node1_y, node1_z}
                                                                                                               ):
                                                                                           unit_normal<double>({node0_x, node0_y, node0_z},
                                                                                                               {node1_x, node1_y, node1_z},
                                                                                                               {node2_x, node2_y, node2_z}
                                                                                                               );

        // assemble cartesian basis
        std::vector<double> basis_i(3), basis_j(3), basis_k(3);

        // basis k
        basis_k = unit_normal_vec;

        // basis i
        std::vector<double> c0_coords = {node0_x - centr_x, node0_y - centr_y, node0_z - centr_z};
        basis_i = c0_coords*(1.0/L2_norm_vec(c0_coords)); 
        
        // calculating basis j vector
        basis_j = basis_k%basis_i;

        // first basis vector
        trianit->RealArray(trian_cart_basis_x_tag)[0] = basis_i[0];
        trianit->RealArray(trian_cart_basis_x_tag)[1] = basis_i[1];
        trianit->RealArray(trian_cart_basis_x_tag)[2] = basis_i[2];

        // second basis vector
        trianit->RealArray(trian_cart_basis_y_tag)[0] = basis_j[0];
        trianit->RealArray(trian_cart_basis_y_tag)[1] = basis_j[1];
        trianit->RealArray(trian_cart_basis_y_tag)[2] = basis_j[2];

        // third basis vector
        trianit->RealArray(trian_cart_basis_z_tag)[0] = basis_k[0];
        trianit->RealArray(trian_cart_basis_z_tag)[1] = basis_k[1];
        trianit->RealArray(trian_cart_basis_z_tag)[2] = basis_k[2];
    }

    // ## edges ##
    INMOST::Tag edge_cart_basis_x_tag = ice_mesh->CreateTag("cart basis x edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    INMOST::Tag edge_cart_basis_y_tag = ice_mesh->CreateTag("cart basis y edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    INMOST::Tag edge_cart_basis_z_tag = ice_mesh->CreateTag("cart basis z edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);

    grid_info[gridElemType::Edge]->cart_basis = {edge_cart_basis_x_tag, edge_cart_basis_y_tag, edge_cart_basis_z_tag};

    for (auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        // assemble cartesian basis
        std::vector<double> basis_i(3), basis_j(3), basis_k(3);

        std::vector<double> basis_k_not_unit = {0.0, 0.0, 0.0};

        // gain adj trians for edge
        ElementArray<INMOST::Cell> adj_trians = edgeit->getCells();
        double sum_S = 0.0;

        for (size_t tr_num = 0; tr_num < adj_trians.size(); ++tr_num)
        {
            // gain nodes of current triangle
            ElementArray<INMOST::Node> curr_tr_nodes = adj_trians[tr_num].getNodes();

            // gain node coords
            std::vector<double> node_coords(9);

            for (size_t i = 0; i < 3; i++)
            {
                node_coords[i*3+0] = curr_tr_nodes[i].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[0];
                node_coords[i*3+1] = curr_tr_nodes[i].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[1];
                node_coords[i*3+2] = curr_tr_nodes[i].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[2];
            }

            // compute area of current trian
            double Scurr = trian_square<double>({node_coords[0], node_coords[1], node_coords[2]},
                                                {node_coords[3], node_coords[4], node_coords[5]},
                                                {node_coords[6], node_coords[7], node_coords[8]});

            // gain unit normal of current trian
            INMOST::Tag tr_cart_basis_tag = grid_info[gridElemType::Trian]->cart_basis[2];
            std::vector<double> basis_k_tr = {adj_trians[tr_num].RealArray(trian_cart_basis_z_tag)[0],
                                              adj_trians[tr_num].RealArray(trian_cart_basis_z_tag)[1],
                                              adj_trians[tr_num].RealArray(trian_cart_basis_z_tag)[2]};

            basis_k_not_unit = basis_k_not_unit + basis_k_tr*Scurr;
            sum_S += Scurr;
        }

        basis_k_not_unit = basis_k_not_unit*(1.0/sum_S);
        basis_k = basis_k_not_unit*(1.0/L2_norm_vec(basis_k_not_unit));

        // compute j-basis vector
        ElementArray<INMOST::Node> adj_nodes = edgeit->getNodes();
        std::vector<double> edg0_coords(3), edg1_coords(3);

        INMOST::Tag node_cart_tag = grid_info[gridElemType::Node]->coords[coordType::cart];

        edg0_coords[0] = adj_nodes[0].RealArray(node_cart_tag)[0]; edg0_coords[1] = adj_nodes[0].RealArray(node_cart_tag)[1]; edg0_coords[2] = adj_nodes[0].RealArray(node_cart_tag)[2];
        edg1_coords[0] = adj_nodes[1].RealArray(node_cart_tag)[0]; edg1_coords[1] = adj_nodes[1].RealArray(node_cart_tag)[1]; edg1_coords[2] = adj_nodes[1].RealArray(node_cart_tag)[2]; 

        std::vector<double> basis_j_not_unit = {edg1_coords[0] - edg0_coords[0],
                                                edg1_coords[1] - edg0_coords[1],
                                                edg1_coords[2] - edg0_coords[2]};

        basis_j = basis_j_not_unit*(1.0/L2_norm_vec(basis_j_not_unit));
        
        //compute i-basis vector
        basis_i = basis_j%basis_k;

        // first basis vector
        edgeit->RealArray(edge_cart_basis_x_tag)[0] = basis_i[0];
        edgeit->RealArray(edge_cart_basis_x_tag)[1] = basis_i[1];
        edgeit->RealArray(edge_cart_basis_x_tag)[2] = basis_i[2];

        // second basis vector
        edgeit->RealArray(edge_cart_basis_y_tag)[0] = basis_j[0];
        edgeit->RealArray(edge_cart_basis_y_tag)[1] = basis_j[1];
        edgeit->RealArray(edge_cart_basis_y_tag)[2] = basis_j[2];

        // third basis vector
        edgeit->RealArray(edge_cart_basis_z_tag)[0] = basis_k[0];
        edgeit->RealArray(edge_cart_basis_z_tag)[1] = basis_k[1];
        edgeit->RealArray(edge_cart_basis_z_tag)[2] = basis_k[2];
    }
    
    // ## nodes ##
    INMOST::Tag node_cart_basis_x_tag = ice_mesh->CreateTag("cart basis x node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    INMOST::Tag node_cart_basis_y_tag = ice_mesh->CreateTag("cart basis y node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    INMOST::Tag node_cart_basis_z_tag = ice_mesh->CreateTag("cart basis z node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);

    grid_info[gridElemType::Node]->cart_basis = {node_cart_basis_x_tag, node_cart_basis_y_tag, node_cart_basis_z_tag};

    for (auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        // assemble cartesian basis
        std::vector<double> basis_i(3), basis_j(3), basis_k(3);

        std::vector<double> basis_k_not_unit = {0.0, 0.0, 0.0};

        // gain adj trians for node
        ElementArray<INMOST::Cell> adj_trians = nodeit->getCells();
        double sum_S = 0.0;

        for (size_t tr_num = 0; tr_num < adj_trians.size(); ++tr_num)
        {
            // gain nodes of current triangle
            ElementArray<INMOST::Node> curr_tr_nodes = adj_trians[tr_num].getNodes();

            // gain node coords
            std::vector<double> node_coords(9);

            for (size_t i = 0; i < 3; i++)
            {
                node_coords[i*3+0] = curr_tr_nodes[i].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[0];
                node_coords[i*3+1] = curr_tr_nodes[i].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[1];
                node_coords[i*3+2] = curr_tr_nodes[i].RealArray(grid_info[gridElemType::Node]->coords[coordType::cart])[2];
            }

            // compute area of current trian
            double Scurr = trian_square<double>({node_coords[0], node_coords[1], node_coords[2]},
                                                {node_coords[3], node_coords[4], node_coords[5]},
                                                {node_coords[6], node_coords[7], node_coords[8]});

            // gain unit normal of current trian
            std::vector<double> basis_k_tr = {adj_trians[tr_num].RealArray(trian_cart_basis_z_tag)[0],
                                              adj_trians[tr_num].RealArray(trian_cart_basis_z_tag)[1],
                                              adj_trians[tr_num].RealArray(trian_cart_basis_z_tag)[2]};

            basis_k_not_unit = basis_k_not_unit + basis_k_tr*Scurr;
            sum_S += Scurr;
        }

        basis_k_not_unit = basis_k_not_unit*(1.0/sum_S);
        basis_k = basis_k_not_unit*(1.0/L2_norm_vec(basis_k_not_unit));

        std::vector<double> basis_j_not_unit = {0.0, -basis_k[2], basis_k[1]};
        basis_j = basis_j_not_unit*(1.0/L2_norm_vec(basis_j_not_unit));

        basis_i = basis_j%basis_k;

        // first basis vector
        nodeit->RealArray(node_cart_basis_x_tag)[0] = basis_i[0];
        nodeit->RealArray(node_cart_basis_x_tag)[1] = basis_i[1];
        nodeit->RealArray(node_cart_basis_x_tag)[2] = basis_i[2];

        // second basis vector
        nodeit->RealArray(node_cart_basis_y_tag)[0] = basis_j[0];
        nodeit->RealArray(node_cart_basis_y_tag)[1] = basis_j[1];
        nodeit->RealArray(node_cart_basis_y_tag)[2] = basis_j[2];

        // third basis vector
        nodeit->RealArray(node_cart_basis_z_tag)[0] = basis_k[0];
        nodeit->RealArray(node_cart_basis_z_tag)[1] = basis_k[1];
        nodeit->RealArray(node_cart_basis_z_tag)[2] = basis_k[2];
    }
    BARRIER
}

void IceMesh::AssembleBasisData()
{
    AssembleGeoElementBasis();
    AssembleCartesianElementBasis();
}