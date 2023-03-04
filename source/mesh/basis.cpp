#include "mesh.hpp"

using namespace std;
using namespace INMOST;

namespace SIMUG
{

void IceMesh::AssembleGeoElementBasis()
{
    // ## nodes ##
    INMOST::Tag node_geo_basis_x_tag = ice_mesh->CreateTag("geo basis x node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    INMOST::Tag node_geo_basis_y_tag = ice_mesh->CreateTag("geo basis y node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    INMOST::Tag node_geo_basis_z_tag = ice_mesh->CreateTag("geo basis z node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);

    for(auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            // node geo coords
            double lon = nodeit->RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::geo])[0];
            double lat = nodeit->RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::geo])[1];

            // first basis vector
            nodeit->RealArray(node_geo_basis_x_tag)[0] = (mesh_info.surface_type == mesh::surfType::plane)?1.0:-sin(lon);
            nodeit->RealArray(node_geo_basis_x_tag)[1] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lon);
            nodeit->RealArray(node_geo_basis_x_tag)[2] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:0.0;

            // second basis vector
            nodeit->RealArray(node_geo_basis_y_tag)[0] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:-sin(lat)*cos(lon);
            nodeit->RealArray(node_geo_basis_y_tag)[1] = (mesh_info.surface_type == mesh::surfType::plane)?1.0:-sin(lat)*sin(lon);
            nodeit->RealArray(node_geo_basis_y_tag)[2] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lat);

            // third basis vector
            nodeit->RealArray(node_geo_basis_z_tag)[0] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lat)*cos(lon);
            nodeit->RealArray(node_geo_basis_z_tag)[1] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lat)*sin(lon);
            nodeit->RealArray(node_geo_basis_z_tag)[2] = (mesh_info.surface_type == mesh::surfType::plane)?1.0:sin(lat);
        }
    }

    //exchange node basis
    ice_mesh->ExchangeData(node_geo_basis_x_tag, INMOST::NODE, 0);
    ice_mesh->ExchangeData(node_geo_basis_y_tag, INMOST::NODE, 0);
    ice_mesh->ExchangeData(node_geo_basis_z_tag, INMOST::NODE, 0);

    grid_info[mesh::gridElemType::Node]->geo_basis = {node_geo_basis_x_tag, node_geo_basis_y_tag, node_geo_basis_z_tag};

    // ## edges ##
    INMOST::Tag edge_geo_basis_x_tag = ice_mesh->CreateTag("geo basis x edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    INMOST::Tag edge_geo_basis_y_tag = ice_mesh->CreateTag("geo basis y edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    INMOST::Tag edge_geo_basis_z_tag = ice_mesh->CreateTag("geo basis z edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);

    for(auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        if (edgeit->GetStatus() != Element::Ghost)
        {
            // edge geo coords
            double lon = edgeit->RealArray(grid_info[mesh::gridElemType::Edge]->coords[coord::coordType::geo])[0];
            double lat = edgeit->RealArray(grid_info[mesh::gridElemType::Edge]->coords[coord::coordType::geo])[1];

            // first basis vector
            edgeit->RealArray(edge_geo_basis_x_tag)[0] = (mesh_info.surface_type == mesh::surfType::plane)?1.0:-sin(lon);
            edgeit->RealArray(edge_geo_basis_x_tag)[1] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lon);
            edgeit->RealArray(edge_geo_basis_x_tag)[2] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:0.0;

            // second basis vector
            edgeit->RealArray(edge_geo_basis_y_tag)[0] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:-sin(lat)*cos(lon);
            edgeit->RealArray(edge_geo_basis_y_tag)[1] = (mesh_info.surface_type == mesh::surfType::plane)?1.0:-sin(lat)*sin(lon);
            edgeit->RealArray(edge_geo_basis_y_tag)[2] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lat);

            // third basis vector
            edgeit->RealArray(edge_geo_basis_z_tag)[0] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lat)*cos(lon);
            edgeit->RealArray(edge_geo_basis_z_tag)[1] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lat)*sin(lon);
            edgeit->RealArray(edge_geo_basis_z_tag)[2] = (mesh_info.surface_type == mesh::surfType::plane)?1.0:sin(lat);
        }
    }

    // exchange edge basis
    ice_mesh->ExchangeData(edge_geo_basis_x_tag, INMOST::FACE, 0);
    ice_mesh->ExchangeData(edge_geo_basis_y_tag, INMOST::FACE, 0);
    ice_mesh->ExchangeData(edge_geo_basis_z_tag, INMOST::FACE, 0);

    grid_info[mesh::gridElemType::Edge]->geo_basis = {edge_geo_basis_x_tag, edge_geo_basis_y_tag, edge_geo_basis_z_tag};

    // ## trians ##
    INMOST::Tag trian_geo_basis_x_tag = ice_mesh->CreateTag("geo basis x trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    INMOST::Tag trian_geo_basis_y_tag = ice_mesh->CreateTag("geo basis y trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    INMOST::Tag trian_geo_basis_z_tag = ice_mesh->CreateTag("geo basis z trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);


    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            // trian geo coords
            double lon = trianit->RealArray(grid_info[mesh::gridElemType::Trian]->coords[coord::coordType::geo])[0];
            double lat = trianit->RealArray(grid_info[mesh::gridElemType::Trian]->coords[coord::coordType::geo])[1];

            // first basis vector
            trianit->RealArray(trian_geo_basis_x_tag)[0] = (mesh_info.surface_type == mesh::surfType::plane)?1.0:-sin(lon);
            trianit->RealArray(trian_geo_basis_x_tag)[1] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lon);
            trianit->RealArray(trian_geo_basis_x_tag)[2] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:0.0;

            // second basis vector
            trianit->RealArray(trian_geo_basis_y_tag)[0] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:-sin(lat)*cos(lon);
            trianit->RealArray(trian_geo_basis_y_tag)[1] = (mesh_info.surface_type == mesh::surfType::plane)?1.0:-sin(lat)*sin(lon);
            trianit->RealArray(trian_geo_basis_y_tag)[2] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lat);

            // third basis vector
            trianit->RealArray(trian_geo_basis_z_tag)[0] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lat)*cos(lon);
            trianit->RealArray(trian_geo_basis_z_tag)[1] = (mesh_info.surface_type == mesh::surfType::plane)?0.0:cos(lat)*sin(lon);
            trianit->RealArray(trian_geo_basis_z_tag)[2] = (mesh_info.surface_type == mesh::surfType::plane)?1.0:sin(lat);
        }
    }

    //exchange trian basis
    ice_mesh->ExchangeData(trian_geo_basis_x_tag, INMOST::CELL, 0);
    ice_mesh->ExchangeData(trian_geo_basis_y_tag, INMOST::CELL, 0);
    ice_mesh->ExchangeData(trian_geo_basis_z_tag, INMOST::CELL, 0);

    grid_info[mesh::gridElemType::Trian]->geo_basis = {trian_geo_basis_x_tag, trian_geo_basis_y_tag, trian_geo_basis_z_tag};

    BARRIER
}

void IceMesh::AssembleCartesianElementBasis()
{
    // ## triangles ##
    INMOST::Tag trian_cart_basis_x_tag = ice_mesh->CreateTag("cart basis x trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    INMOST::Tag trian_cart_basis_y_tag = ice_mesh->CreateTag("cart basis y trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    INMOST::Tag trian_cart_basis_z_tag = ice_mesh->CreateTag("cart basis z trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);

    for(auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            // gain adj nodes to triangle
            ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

            // gain cartesian coords of nodes
            double node0_x = adj_nodes[0].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[0];
            double node0_y = adj_nodes[0].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[1];
            double node0_z = adj_nodes[0].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[2];

            double node1_x = adj_nodes[1].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[0];
            double node1_y = adj_nodes[1].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[1];
            double node1_z = adj_nodes[1].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[2];

            double node2_x = adj_nodes[2].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[0];
            double node2_y = adj_nodes[2].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[1];
            double node2_z = adj_nodes[2].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[2];

            // calculate centroid coordinates
            double centr_x = (node0_x + node1_x + node2_x)/3.0;
            double centr_y = (node0_y + node1_y + node2_y)/3.0;
            double centr_z = (node0_z + node1_z + node2_z)/3.0;

            // calculate unit normal to triangle
            std::vector<double> unit_normal_vec = (mesh_info.surface_type == mesh::surfType::basin)? unit_normal<double>({node0_x, node0_y, node0_z},
                                                                                                                   {node2_x, node2_y, node2_z},
                                                                                                                   {node1_x, node1_y, node1_z}
                                                                                                                   ):
                                                                                               unit_normal<double>({node0_x, node0_y, node0_z},
                                                                                                                   {node1_x, node1_y, node1_z},
                                                                                                                   {node2_x, node2_y, node2_z}
                                                                                                                   );

            // assemble cartesian basis
            std::vector<double> basis_i(3), basis_j(3), basis_k(3);

            basis_k = unit_normal_vec;

            //basis_k = unit_normal_vec;

            // fix unit normal if numeration is bad for plane
            /*
            if (mesh_info.surface_type == mesh::surfType::plane)
            {
                if (basis_k[2] < 0.0)
                {
                    basis_k = (-1.0)*basis_k;
                }
            }
            */

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
    }

    //exchange triangle basis
    ice_mesh->ExchangeData(trian_cart_basis_x_tag, INMOST::CELL, 0);
    ice_mesh->ExchangeData(trian_cart_basis_y_tag, INMOST::CELL, 0);
    ice_mesh->ExchangeData(trian_cart_basis_z_tag, INMOST::CELL, 0);

    grid_info[mesh::gridElemType::Trian]->cart_basis = {trian_cart_basis_x_tag, trian_cart_basis_y_tag, trian_cart_basis_z_tag};

    BARRIER

    // ## edges ##
    INMOST::Tag edge_cart_basis_x_tag = ice_mesh->CreateTag("cart basis x edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    INMOST::Tag edge_cart_basis_y_tag = ice_mesh->CreateTag("cart basis y edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    INMOST::Tag edge_cart_basis_z_tag = ice_mesh->CreateTag("cart basis z edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 3);
    INMOST::Tag is_edge_bnd = grid_info[mesh::gridElemType::Edge]->is_bnd;
    INMOST::Tag node_id_tag = grid_info[mesh::gridElemType::Node]->id;

    for (auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        if (edgeit->GetStatus() != Element::Ghost)
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
                    node_coords[i*3+0] = curr_tr_nodes[i].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[0];
                    node_coords[i*3+1] = curr_tr_nodes[i].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[1];
                    node_coords[i*3+2] = curr_tr_nodes[i].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[2];
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

            // compute j-basis vector
            ElementArray<INMOST::Node> adj_nodes = edgeit->getNodes();
            std::vector<double> edg0_coords(3), edg1_coords(3);

            INMOST::Tag node_cart_tag = grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart];

            edg0_coords[0] = adj_nodes[0].RealArray(node_cart_tag)[0]; edg0_coords[1] = adj_nodes[0].RealArray(node_cart_tag)[1]; edg0_coords[2] = adj_nodes[0].RealArray(node_cart_tag)[2];
            edg1_coords[0] = adj_nodes[1].RealArray(node_cart_tag)[0]; edg1_coords[1] = adj_nodes[1].RealArray(node_cart_tag)[1]; edg1_coords[2] = adj_nodes[1].RealArray(node_cart_tag)[2]; 

            std::vector<double> basis_j_not_unit = {edg1_coords[0] - edg0_coords[0],
                                                    edg1_coords[1] - edg0_coords[1],
                                                    edg1_coords[2] - edg0_coords[2]};

            basis_j = basis_j_not_unit*(1.0/L2_norm_vec(basis_j_not_unit));

            //compute i-basis vector
            basis_i = basis_j%basis_k;
            
            // make sure that basis_i vector looks outside of triangle for bnd edge
            if (edgeit->Integer(is_edge_bnd) == 1)
            {
                // get nodes of adjacent trian
                ElementArray<INMOST::Node> adj_trian_nodes = edgeit->getCells()[0]->getNodes();
                int index;
                
                if ((adj_trian_nodes[0]->Integer(node_id_tag) != edgeit->getNodes()[0]->Integer(node_id_tag)) and 
                    (adj_trian_nodes[0]->Integer(node_id_tag) != edgeit->getNodes()[1]->Integer(node_id_tag)))
                {
                    index = 0;
                }
                else if ((adj_trian_nodes[1]->Integer(node_id_tag) != edgeit->getNodes()[0]->Integer(node_id_tag)) and 
                         (adj_trian_nodes[1]->Integer(node_id_tag) != edgeit->getNodes()[1]->Integer(node_id_tag)))
                {
                    index = 1;
                }
                else
                {
                    index = 2;
                }

                std::vector<double> start_vec = 
                {
                    adj_trian_nodes[index]->RealArray(node_cart_tag)[0],
                    adj_trian_nodes[index]->RealArray(node_cart_tag)[1],
                    adj_trian_nodes[index]->RealArray(node_cart_tag)[2],
                };

                std::vector<double> end_vec = 
                {
                    adj_trian_nodes[(index + 1)%3]->RealArray(node_cart_tag)[0],
                    adj_trian_nodes[(index + 1)%3]->RealArray(node_cart_tag)[1],
                    adj_trian_nodes[(index + 1)%3]->RealArray(node_cart_tag)[2],
                };
                double sc_mult = (end_vec - start_vec)*basis_i;

                if (sc_mult < 0)
                {
                    basis_i = (-1.0)*basis_i;
                    basis_j = (-1.0)*basis_j;
                }
            }

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
    }

    //exchange edge basis
    ice_mesh->ExchangeData(edge_cart_basis_x_tag, INMOST::FACE, 0);
    ice_mesh->ExchangeData(edge_cart_basis_y_tag, INMOST::FACE, 0);
    ice_mesh->ExchangeData(edge_cart_basis_z_tag, INMOST::FACE, 0);

    grid_info[mesh::gridElemType::Edge]->cart_basis = {edge_cart_basis_x_tag, edge_cart_basis_y_tag, edge_cart_basis_z_tag};

    BARRIER
    
    // ## nodes ##
    INMOST::Tag node_cart_basis_x_tag = ice_mesh->CreateTag("cart basis x node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    INMOST::Tag node_cart_basis_y_tag = ice_mesh->CreateTag("cart basis y node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);
    INMOST::Tag node_cart_basis_z_tag = ice_mesh->CreateTag("cart basis z node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 3);

    for (auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
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
                    node_coords[i*3+0] = curr_tr_nodes[i].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[0];
                    node_coords[i*3+1] = curr_tr_nodes[i].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[1];
                    node_coords[i*3+2] = curr_tr_nodes[i].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[2];
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
    }

    // exchange node basis
    ice_mesh->ExchangeData(node_cart_basis_x_tag, INMOST::NODE, 0);
    ice_mesh->ExchangeData(node_cart_basis_y_tag, INMOST::NODE, 0);
    ice_mesh->ExchangeData(node_cart_basis_z_tag, INMOST::NODE, 0);

    grid_info[mesh::gridElemType::Node]->cart_basis = {node_cart_basis_x_tag, node_cart_basis_y_tag, node_cart_basis_z_tag};

    BARRIER
}

void IceMesh::AssembleGeoToElementTransitionMatricies()
{
    // ## nodes ##
    INMOST::Tag geo_to_elem_node_tag = ice_mesh->CreateTag("geo to elem trans matr node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 4);
    INMOST::Tag elem_to_geo_node_tag = ice_mesh->CreateTag("elem to geo trans matr node", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 4);
    grid_info[mesh::gridElemType::Node]->trans_matr_from_geo_to_elem = geo_to_elem_node_tag;
    grid_info[mesh::gridElemType::Node]->trans_matr_from_elem_to_geo = elem_to_geo_node_tag;
    std::vector<INMOST::Tag> geo_basis_node = grid_info[mesh::gridElemType::Node]->geo_basis;
    std::vector<INMOST::Tag> cart_basis_node = grid_info[mesh::gridElemType::Node]->cart_basis;

    for (auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            // gain geo and cart basis vectors
            std::vector<double> basis_geo_x = {nodeit->RealArray(geo_basis_node[0])[0],
                                               nodeit->RealArray(geo_basis_node[0])[1],
                                               nodeit->RealArray(geo_basis_node[0])[2]};

            std::vector<double> basis_geo_y = {nodeit->RealArray(geo_basis_node[1])[0],
                                               nodeit->RealArray(geo_basis_node[1])[1],
                                               nodeit->RealArray(geo_basis_node[1])[2]};

            std::vector<double> basis_cart_x = {nodeit->RealArray(cart_basis_node[0])[0],
                                                nodeit->RealArray(cart_basis_node[0])[1],
                                                nodeit->RealArray(cart_basis_node[0])[2]};

            std::vector<double> basis_cart_y = {nodeit->RealArray(cart_basis_node[1])[0],
                                                nodeit->RealArray(cart_basis_node[1])[1],
                                                nodeit->RealArray(cart_basis_node[1])[2]};

            // assemble 2x2 forward and backward transition matrics
            std::vector<std::vector<double>> forward_matr = {{basis_geo_x*basis_cart_x, basis_geo_x*basis_cart_y},
                                                            {basis_geo_y*basis_cart_x, basis_geo_y*basis_cart_y}};

            std::vector<std::vector<double>> inverse_matr = inv(forward_matr);

            // assign results
            nodeit->RealArray(geo_to_elem_node_tag)[0] = inverse_matr[0][0];
            nodeit->RealArray(geo_to_elem_node_tag)[1] = inverse_matr[0][1];
            nodeit->RealArray(geo_to_elem_node_tag)[2] = inverse_matr[1][0];
            nodeit->RealArray(geo_to_elem_node_tag)[3] = inverse_matr[1][1];

            nodeit->RealArray(elem_to_geo_node_tag)[0] = forward_matr[0][0];
            nodeit->RealArray(elem_to_geo_node_tag)[1] = forward_matr[0][1];
            nodeit->RealArray(elem_to_geo_node_tag)[2] = forward_matr[1][0];
            nodeit->RealArray(elem_to_geo_node_tag)[3] = forward_matr[1][1];
        }
    }

    // exchange node matricies
    ice_mesh->ExchangeData(geo_to_elem_node_tag, INMOST::NODE, 0);
    ice_mesh->ExchangeData(elem_to_geo_node_tag, INMOST::NODE, 0);

    // ## edges ##
    INMOST::Tag geo_to_elem_edge_tag = ice_mesh->CreateTag("geo to elem trans matr edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 4);
    INMOST::Tag elem_to_geo_edge_tag = ice_mesh->CreateTag("elem to geo trans matr edge", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 4);
    grid_info[mesh::gridElemType::Edge]->trans_matr_from_geo_to_elem = geo_to_elem_edge_tag;
    grid_info[mesh::gridElemType::Edge]->trans_matr_from_elem_to_geo = elem_to_geo_edge_tag;
    std::vector<INMOST::Tag> geo_basis_edge = grid_info[mesh::gridElemType::Edge]->geo_basis;
    std::vector<INMOST::Tag> cart_basis_edge = grid_info[mesh::gridElemType::Edge]->cart_basis;

    for (auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        //if (edgeit->GetStatus() != Element::Ghost)
        {
            // gain geo and cart basis vectors
            std::vector<double> basis_geo_x = {edgeit->RealArray(geo_basis_edge[0])[0],
                                               edgeit->RealArray(geo_basis_edge[0])[1],
                                               edgeit->RealArray(geo_basis_edge[0])[2]};

            std::vector<double> basis_geo_y = {edgeit->RealArray(geo_basis_edge[1])[0],
                                               edgeit->RealArray(geo_basis_edge[1])[1],
                                               edgeit->RealArray(geo_basis_edge[1])[2]};

            std::vector<double> basis_cart_x = {edgeit->RealArray(cart_basis_edge[0])[0],
                                                edgeit->RealArray(cart_basis_edge[0])[1],
                                                edgeit->RealArray(cart_basis_edge[0])[2]};

            std::vector<double> basis_cart_y = {edgeit->RealArray(cart_basis_edge[1])[0],
                                                edgeit->RealArray(cart_basis_edge[1])[1],
                                                edgeit->RealArray(cart_basis_edge[1])[2]};

            // assemble 2x2 forward and backward transition matrics
            std::vector<std::vector<double>> forward_matr = {{basis_geo_x*basis_cart_x, basis_geo_x*basis_cart_y},
                                                             {basis_geo_y*basis_cart_x, basis_geo_y*basis_cart_y}};

            std::vector<std::vector<double>> inverse_matr = inv(forward_matr);

            // assign results
            edgeit->RealArray(geo_to_elem_edge_tag)[0] = inverse_matr[0][0];
            edgeit->RealArray(geo_to_elem_edge_tag)[1] = inverse_matr[0][1];
            edgeit->RealArray(geo_to_elem_edge_tag)[2] = inverse_matr[1][0];
            edgeit->RealArray(geo_to_elem_edge_tag)[3] = inverse_matr[1][1];

            edgeit->RealArray(elem_to_geo_edge_tag)[0] = forward_matr[0][0];
            edgeit->RealArray(elem_to_geo_edge_tag)[1] = forward_matr[0][1];
            edgeit->RealArray(elem_to_geo_edge_tag)[2] = forward_matr[1][0];
            edgeit->RealArray(elem_to_geo_edge_tag)[3] = forward_matr[1][1];
        }
    }

    // exchange edge matricies
    ice_mesh->ExchangeData(geo_to_elem_edge_tag, INMOST::FACE, 0);
    ice_mesh->ExchangeData(elem_to_geo_edge_tag, INMOST::FACE, 0);

    // ## triangles ##
    INMOST::Tag geo_to_elem_trian_tag = ice_mesh->CreateTag("geo to elem trans matr trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 4);
    INMOST::Tag elem_to_geo_trian_tag = ice_mesh->CreateTag("elem to geo trans matr trian", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 4);
    grid_info[mesh::gridElemType::Trian]->trans_matr_from_geo_to_elem = geo_to_elem_trian_tag;
    grid_info[mesh::gridElemType::Trian]->trans_matr_from_elem_to_geo = elem_to_geo_trian_tag;
    std::vector<INMOST::Tag> geo_basis_trian = grid_info[mesh::gridElemType::Trian]->geo_basis;
    std::vector<INMOST::Tag> cart_basis_trian = grid_info[mesh::gridElemType::Trian]->cart_basis;

    for (auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            // gain geo and cart basis vectors
            std::vector<double> basis_geo_x = {trianit->RealArray(geo_basis_trian[0])[0],
                                               trianit->RealArray(geo_basis_trian[0])[1],
                                               trianit->RealArray(geo_basis_trian[0])[2]};

            std::vector<double> basis_geo_y = {trianit->RealArray(geo_basis_trian[1])[0],
                                               trianit->RealArray(geo_basis_trian[1])[1],
                                               trianit->RealArray(geo_basis_trian[1])[2]};

            std::vector<double> basis_cart_x = {trianit->RealArray(cart_basis_trian[0])[0],
                                                trianit->RealArray(cart_basis_trian[0])[1],
                                                trianit->RealArray(cart_basis_trian[0])[2]};

            std::vector<double> basis_cart_y = {trianit->RealArray(cart_basis_trian[1])[0],
                                                trianit->RealArray(cart_basis_trian[1])[1],
                                                trianit->RealArray(cart_basis_trian[1])[2]};

            // assemble 2x2 forward and backward transition matrics
            std::vector<std::vector<double>> forward_matr = {{basis_geo_x*basis_cart_x, basis_geo_x*basis_cart_y},
                                                            {basis_geo_y*basis_cart_x, basis_geo_y*basis_cart_y}};

            std::vector<std::vector<double>> inverse_matr = inv(forward_matr);

            // assign results
            trianit->RealArray(geo_to_elem_trian_tag)[0] = inverse_matr[0][0];
            trianit->RealArray(geo_to_elem_trian_tag)[1] = inverse_matr[0][1];
            trianit->RealArray(geo_to_elem_trian_tag)[2] = inverse_matr[1][0];
            trianit->RealArray(geo_to_elem_trian_tag)[3] = inverse_matr[1][1];

            trianit->RealArray(elem_to_geo_trian_tag)[0] = forward_matr[0][0];
            trianit->RealArray(elem_to_geo_trian_tag)[1] = forward_matr[0][1];
            trianit->RealArray(elem_to_geo_trian_tag)[2] = forward_matr[1][0];
            trianit->RealArray(elem_to_geo_trian_tag)[3] = forward_matr[1][1];
        }
    }

    // exchange triangle matricies
    ice_mesh->ExchangeData(geo_to_elem_trian_tag, INMOST::CELL, 0);
    ice_mesh->ExchangeData(elem_to_geo_trian_tag, INMOST::CELL, 0);

    BARRIER
}

void IceMesh::AssembleElementToElementTransitionMatricies()
{
    // get all element cartesian basis tags
    std::vector<INMOST::Tag> cart_basis_node_tag = grid_info[mesh::gridElemType::Node]->cart_basis;
    std::vector<INMOST::Tag> cart_basis_edge_tag = grid_info[mesh::gridElemType::Edge]->cart_basis;
    std::vector<INMOST::Tag> cart_basis_trian_tag = grid_info[mesh::gridElemType::Trian]->cart_basis;

    // get id tag for every element
    INMOST::Tag node_id_tag = grid_info[mesh::gridElemType::Node]->id;
    INMOST::Tag edge_id_tag = grid_info[mesh::gridElemType::Edge]->id;
    INMOST::Tag trian_id_tag = grid_info[mesh::gridElemType::Trian]->id;

    // get all transition matricies tags
    std::vector<INMOST::Tag>& matr_node_to_edge_tag = grid_info[mesh::gridElemType::Node]->GetTransMatrToEdge();
    std::vector<INMOST::Tag>& matr_node_to_trian_tag = grid_info[mesh::gridElemType::Node]->GetTransMatrToTrian();

    std::vector<INMOST::Tag>& matr_edge_to_node_tag = grid_info[mesh::gridElemType::Edge]->GetTransMatrToNode();
    std::vector<INMOST::Tag>& matr_edge_to_trian_tag = grid_info[mesh::gridElemType::Edge]->GetTransMatrToTrian();

    std::vector<INMOST::Tag>& matr_trian_to_node_tag = grid_info[mesh::gridElemType::Trian]->GetTransMatrToNode();
    std::vector<INMOST::Tag>& matr_trian_to_edge_tag = grid_info[mesh::gridElemType::Trian]->GetTransMatrToEdge();

    // create matricies tags for all elements and mute instantly
    for (int i = 0; i < MAX_NUM_ADJ_EDGES; ++i)
    {
        matr_node_to_edge_tag.push_back(ice_mesh->CreateTag("matr node to edge " + to_string(i), INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 4));
        ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"matr node to edge " + to_string(i), "nosave");
    }

    for (int i = 0; i < MAX_NUM_ADJ_TRIANS; ++i)
    {
        matr_node_to_trian_tag.push_back(ice_mesh->CreateTag("matr node to trian " + to_string(i), INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 4));
        ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"matr node to trian " + to_string(i), "nosave");
    }

    for (int i = 0; i < 3; ++i)
    {
        matr_edge_to_node_tag.push_back(ice_mesh->CreateTag("matr edge to node " + to_string(i), INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 4));
        ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"matr edge to node " + to_string(i), "nosave");
    }
    
    for (int i = 0; i < 2; ++i)
    {
        matr_edge_to_trian_tag.push_back(ice_mesh->CreateTag("matr edge to trian " + to_string(i), INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 4));
        ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"matr edge to trian " + to_string(i), "nosave");
    }

    for (int i = 0; i < 3; ++i)
    {
        matr_trian_to_node_tag.push_back(ice_mesh->CreateTag("matr trian to node " + to_string(i), INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 4));
        ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"matr trian to node " + to_string(i), "nosave");
    }
    
    for (int i = 0; i < 3; ++i)
    {
        matr_trian_to_edge_tag.push_back(ice_mesh->CreateTag("matr trian to edge " + to_string(i), INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 4));
        ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"matr trian to edge " + to_string(i), "nosave");
    }

    // ### "triangle <-> node" and "triangle <-> edge" ###
    for (auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        // get cartesian basis of current triangle
        std::vector<double> trian_basis_x = {trianit->RealArray(cart_basis_trian_tag[0])[0],
                                             trianit->RealArray(cart_basis_trian_tag[0])[1],
                                             trianit->RealArray(cart_basis_trian_tag[0])[2]};

        std::vector<double> trian_basis_y = {trianit->RealArray(cart_basis_trian_tag[1])[0],
                                             trianit->RealArray(cart_basis_trian_tag[1])[1],
                                             trianit->RealArray(cart_basis_trian_tag[1])[2]};

        std::vector<double> trian_basis_z = {trianit->RealArray(cart_basis_trian_tag[2])[0],
                                             trianit->RealArray(cart_basis_trian_tag[2])[1],
                                             trianit->RealArray(cart_basis_trian_tag[2])[2]};

        // # assemble matricies "trian <-> node" #

        // get adjacent nodes
        ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

        for (size_t i = 0; i < 3; ++i)
        {
            // get current node cart basis
            std::vector<double> node_basis_x = {adj_nodes[i].RealArray(cart_basis_node_tag[0])[0],
                                                adj_nodes[i].RealArray(cart_basis_node_tag[0])[1],
                                                adj_nodes[i].RealArray(cart_basis_node_tag[0])[2]};

            std::vector<double> node_basis_y = {adj_nodes[i].RealArray(cart_basis_node_tag[1])[0],
                                                adj_nodes[i].RealArray(cart_basis_node_tag[1])[1],
                                                adj_nodes[i].RealArray(cart_basis_node_tag[1])[2]};

            // assemble "trian -> node" transition matrix
            std::vector<std::vector<double>> node_to_trian_matrix = {{trian_basis_x*node_basis_x, trian_basis_x*node_basis_y},
                                                                     {trian_basis_y*node_basis_x, trian_basis_y*node_basis_y}};

            // assemble inverse transition matrix "node -> trian" = "trian -> node"^-1
            std::vector<std::vector<double>> trian_to_node_matrix = inv(node_to_trian_matrix);


            // store matrix for triangles
            trianit->RealArray(matr_trian_to_node_tag[i])[0] = trian_to_node_matrix[0][0];
            trianit->RealArray(matr_trian_to_node_tag[i])[1] = trian_to_node_matrix[0][1];
            trianit->RealArray(matr_trian_to_node_tag[i])[2] = trian_to_node_matrix[1][0];
            trianit->RealArray(matr_trian_to_node_tag[i])[3] = trian_to_node_matrix[1][1];

            // store matrix for nodes
            ElementArray<INMOST::Cell> adj_trians_for_curr_node = adj_nodes[i].getCells();

            for (size_t j = 0; j < adj_trians_for_curr_node.size(); ++j)
            {
                if (adj_trians_for_curr_node[j].Integer(trian_id_tag) == trianit->Integer(trian_id_tag))
                {
                    adj_nodes[i].RealArray(matr_node_to_trian_tag[j])[0] = node_to_trian_matrix[0][0];
                    adj_nodes[i].RealArray(matr_node_to_trian_tag[j])[1] = node_to_trian_matrix[0][1];
                    adj_nodes[i].RealArray(matr_node_to_trian_tag[j])[2] = node_to_trian_matrix[1][0];
                    adj_nodes[i].RealArray(matr_node_to_trian_tag[j])[3] = node_to_trian_matrix[1][1];

                    break;
                }
            }
        }


        // # assemble matricies "trian <-> edge" #

        // get adjacent edges
        ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

        for (size_t i = 0; i < 3; ++i)
        {
            // get full basis of current edge
            std::vector<double> edge_basis_x = {adj_edges[i].RealArray(cart_basis_edge_tag[0])[0],
                                                adj_edges[i].RealArray(cart_basis_edge_tag[0])[1],
                                                adj_edges[i].RealArray(cart_basis_edge_tag[0])[2]};

            std::vector<double> edge_basis_y = {adj_edges[i].RealArray(cart_basis_edge_tag[1])[0],
                                                adj_edges[i].RealArray(cart_basis_edge_tag[1])[1],
                                                adj_edges[i].RealArray(cart_basis_edge_tag[1])[2]};

            std::vector<double> edge_basis_z = {adj_edges[i].RealArray(cart_basis_edge_tag[2])[0],
                                                adj_edges[i].RealArray(cart_basis_edge_tag[2])[1],
                                                adj_edges[i].RealArray(cart_basis_edge_tag[2])[2]};

            // assemble 3-D rotation matrix "trian -> edge"
            std::vector<std::vector<double>> trian_to_edge_3d = {{edge_basis_x*trian_basis_x, edge_basis_x*trian_basis_y, edge_basis_x*trian_basis_z},
                                                                 {edge_basis_y*trian_basis_x, edge_basis_y*trian_basis_y, edge_basis_y*trian_basis_z},
                                                                 {edge_basis_z*trian_basis_x, edge_basis_z*trian_basis_y, edge_basis_z*trian_basis_z}};

            // compute trian basis vecs in edge coordinates
            std::vector<double> trian_basis_x_old = trian_to_edge_3d*std::vector<double>{1.0, 0.0, 0.0};
            std::vector<double> trian_basis_y_old = trian_to_edge_3d*std::vector<double>{0.0, 1.0, 0.0};
            std::vector<double> trian_basis_z_old = trian_to_edge_3d*std::vector<double>{0.0, 0.0, 1.0};


            // compute rotation angle
            double ang_rot = angle_vecs(edge_basis_z, trian_basis_z);

            // assemble rotation matrix around y edge basis vector for positive and negative angle
            std::vector<std::vector<double>> rotation_3d_pos = {{ cos(ang_rot), 0.0, sin(ang_rot)},
                                                                { 0.0         , 1.0, 0.0         },
                                                                {-sin(ang_rot), 0.0, cos(ang_rot)}};

            std::vector<std::vector<double>> rotation_3d_neg = {{ cos(-ang_rot), 0.0, sin(-ang_rot)},
                                                                { 0.0          , 1.0, 0.0          },
                                                                {-sin(-ang_rot), 0.0, cos(-ang_rot)}};

            // perform rotation of z basis vector and chose the lowest angle between new z vec and z-edge basis vec
            std::vector<double> trian_basis_z_new_pos = rotation_3d_pos*trian_basis_z_old;
            std::vector<double> trian_basis_z_new_neg = rotation_3d_neg*trian_basis_z_old;

            std::vector<std::vector<double>> rotation_3d;

            if (angle_vecs(trian_basis_z_new_pos, std::vector<double>{0.0, 0.0, 1.0}) <=
                angle_vecs(trian_basis_z_new_neg, std::vector<double>{0.0, 0.0, 1.0}))
                rotation_3d = rotation_3d_pos;
            else
                rotation_3d = rotation_3d_neg;

            // perform 3d rotation of trian x and y basis vectors
            std::vector<double> trian_basis_x_new = rotation_3d*trian_basis_x_old;
            std::vector<double> trian_basis_y_new = rotation_3d*trian_basis_y_old;

            // finally assemble "trian <-> edge" transition matricies
            std::vector<std::vector<double>> trian_to_edge_matrix = {{std::vector<double>{1.0, 0.0, 0.0}*trian_basis_x_new, std::vector<double>{1.0, 0.0, 0.0}*trian_basis_y_new},
                                                                     {std::vector<double>{0.0, 1.0, 0.0}*trian_basis_x_new, std::vector<double>{0.0, 1.0, 0.0}*trian_basis_y_new}};

            std::vector<std::vector<double>> edge_to_trian_matrix = inv(trian_to_edge_matrix);

            // store matrix for trians
            trianit->RealArray(matr_trian_to_edge_tag[i])[0] = trian_to_edge_matrix[0][0];
            trianit->RealArray(matr_trian_to_edge_tag[i])[1] = trian_to_edge_matrix[0][1];
            trianit->RealArray(matr_trian_to_edge_tag[i])[2] = trian_to_edge_matrix[1][0];
            trianit->RealArray(matr_trian_to_edge_tag[i])[3] = trian_to_edge_matrix[1][1];

            // store matrix for edges
            ElementArray<INMOST::Cell> adj_trians_for_curr_edge = adj_edges[i].getCells();

            for (size_t j = 0; j < adj_trians_for_curr_edge.size(); ++j)
            {
                if (adj_trians_for_curr_edge[j].Integer(trian_id_tag) == trianit->Integer(trian_id_tag))
                {
                    adj_edges[i].RealArray(matr_edge_to_trian_tag[j])[0] = edge_to_trian_matrix[0][0];
                    adj_edges[i].RealArray(matr_edge_to_trian_tag[j])[1] = edge_to_trian_matrix[0][1];
                    adj_edges[i].RealArray(matr_edge_to_trian_tag[j])[2] = edge_to_trian_matrix[1][0];
                    adj_edges[i].RealArray(matr_edge_to_trian_tag[j])[3] = edge_to_trian_matrix[1][1];

                    break;
                }

                if (j == adj_trians_for_curr_edge.size()-1)
                {
                    SIMUG_ERR("cant find adjacent triangle for edge");
                }
            }
        }
    }

    // exchange assembled data
    //for (int i = 0; i < MAX_NUM_ADJ_TRIANS; ++i)
    //    ice_mesh->ExchangeData(matr_node_to_trian_tag[i], INMOST::NODE, 0);
    //    
    //for (int i = 0; i < 2; ++i)
    //    ice_mesh->ExchangeData(matr_edge_to_trian_tag[i], INMOST::FACE, 0);
    //    
    //for (int i = 0; i < 3; ++i)
    //    ice_mesh->ExchangeData(matr_trian_to_node_tag[i], INMOST::CELL, 0);
    //
    //for (int i = 0; i < 3; ++i)
    //    ice_mesh->ExchangeData(matr_trian_to_edge_tag[i], INMOST::CELL, 0);

    BARRIER
    
    // ### "edge <-> node" ###
    for (auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {        
        // get cartesian basis of current edge
        std::vector<double> edge_basis_x = {edgeit->RealArray(cart_basis_edge_tag[0])[0],
                                            edgeit->RealArray(cart_basis_edge_tag[0])[1],
                                            edgeit->RealArray(cart_basis_edge_tag[0])[2]};

        std::vector<double> edge_basis_y = {edgeit->RealArray(cart_basis_edge_tag[1])[0],
                                            edgeit->RealArray(cart_basis_edge_tag[1])[1],
                                            edgeit->RealArray(cart_basis_edge_tag[1])[2]};

        // # assemble matricies "edge <-> node" #

        // get adjacent nodes
        ElementArray<INMOST::Node> adj_nodes = edgeit->getNodes();

        for (size_t i = 0; i < 2; ++i)
        {
            // get current node cart basis
            std::vector<double> node_basis_x = {adj_nodes[i].RealArray(cart_basis_node_tag[0])[0],
                                                adj_nodes[i].RealArray(cart_basis_node_tag[0])[1],
                                                adj_nodes[i].RealArray(cart_basis_node_tag[0])[2]};

            std::vector<double> node_basis_y = {adj_nodes[i].RealArray(cart_basis_node_tag[1])[0],
                                                adj_nodes[i].RealArray(cart_basis_node_tag[1])[1],
                                                adj_nodes[i].RealArray(cart_basis_node_tag[1])[2]};

            // assemble "edge -> node" transition matrix
            std::vector<std::vector<double>> node_to_edge_matrix = {{edge_basis_x*node_basis_x, edge_basis_x*node_basis_y},
                                                                    {edge_basis_y*node_basis_x, edge_basis_y*node_basis_y}};

            // assemble inverse transition matrix "node -> edge" = "edge -> node"^-1
            std::vector<std::vector<double>> edge_to_node_matrix = inv(node_to_edge_matrix);

            // store matrix for edges
            edgeit->RealArray(matr_edge_to_node_tag[i])[0] = edge_to_node_matrix[0][0];
            edgeit->RealArray(matr_edge_to_node_tag[i])[1] = edge_to_node_matrix[0][1];
            edgeit->RealArray(matr_edge_to_node_tag[i])[2] = edge_to_node_matrix[1][0];
            edgeit->RealArray(matr_edge_to_node_tag[i])[3] = edge_to_node_matrix[1][1];

            // store matrix for nodes
            ElementArray<INMOST::Face> adj_edges_for_curr_node = adj_nodes[i].getFaces();
            for (size_t j = 0; j < adj_edges_for_curr_node.size(); ++ j)
            {
                if (adj_edges_for_curr_node[j].Integer(edge_id_tag) == edgeit->Integer(edge_id_tag))
                {
                    adj_nodes[i].RealArray(matr_node_to_edge_tag[j])[0] = node_to_edge_matrix[0][0];
                    adj_nodes[i].RealArray(matr_node_to_edge_tag[j])[1] = node_to_edge_matrix[0][1];
                    adj_nodes[i].RealArray(matr_node_to_edge_tag[j])[2] = node_to_edge_matrix[1][0];
                    adj_nodes[i].RealArray(matr_node_to_edge_tag[j])[3] = node_to_edge_matrix[1][1];

                    break;
                }
            }
        }
    }

    //exchange assembled data
    //for (int i = 0; i < MAX_NUM_ADJ_EDGES; ++i)
    //    ice_mesh->ExchangeData(matr_node_to_edge_tag[i], INMOST::NODE, 0);
    //
    //for (int i = 0; i < 3; ++i)
    //    ice_mesh->ExchangeData(matr_edge_to_node_tag[i], INMOST::FACE, 0);
       
    BARRIER
}

void IceMesh::ComputeNodeCoordsInTrianBasis()
{
    // create tags for local node coords in trian basis
    std::vector<INMOST::Tag>& node_coords_in_trian_basis_tags = grid_info[mesh::gridElemType::Trian]->GetNodeCoordsInTrianBasis();
    for (int i = 0; i < 3; ++i)
    {
        node_coords_in_trian_basis_tags.push_back(ice_mesh->CreateTag("node coords in trian basis " + to_string(i), INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 2));
        ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"node coords in trian basis " + to_string(i), "nosave");
    }

    // get cart coord tag for node
    INMOST::Tag node_cart_coords_tag =  grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart];

    // get trian basis tag for trian
    std::vector<INMOST::Tag> trian_basis_tag = grid_info[mesh::gridElemType::Trian]->cart_basis;

    // calculate node coords in trian basis
    for (auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        //if (trianit->GetStatus() != Element::Ghost)
        //{
            // get triangle basis for current triangle
            std::vector<double> trian_basis_x = 
            {
                trianit->RealArray(trian_basis_tag[0])[0],
                trianit->RealArray(trian_basis_tag[0])[1],
                trianit->RealArray(trian_basis_tag[0])[2]
            };

            std::vector<double> trian_basis_y = 
            {
                trianit->RealArray(trian_basis_tag[1])[0],
                trianit->RealArray(trian_basis_tag[1])[1],
                trianit->RealArray(trian_basis_tag[1])[2]
            };

            // get cartesian coords of current triangle nodes
            ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

            std::vector<double> zero_node_coords = 
            {
                adj_nodes[0].RealArray(node_cart_coords_tag)[0],
                adj_nodes[0].RealArray(node_cart_coords_tag)[1],
                adj_nodes[0].RealArray(node_cart_coords_tag)[2]
            };

            std::vector<double> first_node_coords = 
            {
                adj_nodes[1].RealArray(node_cart_coords_tag)[0],
                adj_nodes[1].RealArray(node_cart_coords_tag)[1],
                adj_nodes[1].RealArray(node_cart_coords_tag)[2]
            };

            std::vector<double> second_node_coords = 
            {
                adj_nodes[2].RealArray(node_cart_coords_tag)[0],
                adj_nodes[2].RealArray(node_cart_coords_tag)[1],
                adj_nodes[2].RealArray(node_cart_coords_tag)[2]
            };

            // calculate centroid coords of current triangle
            std::vector<double> centroid_coords = (zero_node_coords + first_node_coords + second_node_coords)*(1.0/3.0);

            // calculate vectors from centroid to nodes
            std::vector<double> v0 = zero_node_coords - centroid_coords;
            std::vector<double> v1 = first_node_coords - centroid_coords;
            std::vector<double> v2 = second_node_coords - centroid_coords;

            // calculate node coordinates in triangle basis
            std::vector<double> zero_node_coords_in_trian_basis = {v0*trian_basis_x, v0*trian_basis_y};
            std::vector<double> first_node_coords_in_trian_basis = {v1*trian_basis_x, v1*trian_basis_y};
            std::vector<double> second_node_coords_in_trian_basis = {v2*trian_basis_x, v2*trian_basis_y};

            // save local node coords
            trianit->RealArray(node_coords_in_trian_basis_tags[0])[0] = zero_node_coords_in_trian_basis[0];
            trianit->RealArray(node_coords_in_trian_basis_tags[0])[1] = zero_node_coords_in_trian_basis[1];

            trianit->RealArray(node_coords_in_trian_basis_tags[1])[0] = first_node_coords_in_trian_basis[0];
            trianit->RealArray(node_coords_in_trian_basis_tags[1])[1] = first_node_coords_in_trian_basis[1];

            trianit->RealArray(node_coords_in_trian_basis_tags[2])[0] = second_node_coords_in_trian_basis[0];
            trianit->RealArray(node_coords_in_trian_basis_tags[2])[1] = second_node_coords_in_trian_basis[1];
        //}
    }

    // exchange computed data
    for (int i = 0; i < 3; ++i)
        ice_mesh->ExchangeData(node_coords_in_trian_basis_tags[i], INMOST::CELL, 0);
}

void IceMesh::ComputeIfXedgeBasisIsNormalToTrian()
{
    std::vector<INMOST::Tag>& is_x_edge_basis_normal_vec_tags = grid_info[mesh::gridElemType::Trian]->GetIsXedgeBasisIsNormal();

    for (int i = 0; i < 3; ++i)
    {
        is_x_edge_basis_normal_vec_tags.push_back(ice_mesh->CreateTag("is_x_edge_basis_normal_vec " + to_string(i), INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 1));
        ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"is_x_edge_basis_normal_vec " + to_string(i), "nosave");
    }

    // iterate over triangles
    for (auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            // get edges of triangle
            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

            // get nodes of triangle
            ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

            // iterate over adj edges
            for (int ed_num = 0; ed_num < 3; ++ed_num)
            {
                // get nodes of current edge
                ElementArray<INMOST::Node> adj_nodes_for_edge = adj_edges[ed_num].getNodes();

                // get coords of current edge x basis vector
                std::vector<double> ed_basis_x_vec =
                {
                    adj_edges[ed_num].RealArray(grid_info[mesh::gridElemType::Edge]->cart_basis[0])[0],
                    adj_edges[ed_num].RealArray(grid_info[mesh::gridElemType::Edge]->cart_basis[0])[1],
                    adj_edges[ed_num].RealArray(grid_info[mesh::gridElemType::Edge]->cart_basis[0])[2]
                };

                // get Cartesian coordinates of node in current edge (as end of test vector)
                std::vector<double> end_of_test_vec = 
                {
                    adj_nodes_for_edge[0].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[0],
                    adj_nodes_for_edge[0].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[1],
                    adj_nodes_for_edge[0].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[2]
                };

                // get Cartesian coordinates of node outside current edge (as start of test vector)
                std::vector<double> start_of_test_vec;
                
                for (int i = 0; i < 3; ++i)
                {
                    if (adj_nodes_for_edge[0].Integer(grid_info[mesh::gridElemType::Node]->id) != adj_nodes[i].Integer(grid_info[mesh::gridElemType::Node]->id) &&
                        adj_nodes_for_edge[1].Integer(grid_info[mesh::gridElemType::Node]->id) != adj_nodes[i].Integer(grid_info[mesh::gridElemType::Node]->id))
                    {
                        start_of_test_vec = 
                        {
                            adj_nodes[i].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[0],
                            adj_nodes[i].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[1],
                            adj_nodes[i].RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[2]
                        };
                        break;
                    }
                }

                // make test vector that directs outside of triangle 
                std::vector<double> test_vec = end_of_test_vec - start_of_test_vec;
                test_vec = test_vec*(1.0/L2_norm_vec(test_vec));
                
                // figure out if test vector points to the same direction as edge x basis vector
                trianit->Integer(is_x_edge_basis_normal_vec_tags[ed_num]) = (test_vec*ed_basis_x_vec > 0.0) ? 1 : 0;
            }
        }
    }

    // exchange data
    for (int i = 0; i < 3; ++i)
        ice_mesh->ExchangeData(is_x_edge_basis_normal_vec_tags[i], INMOST::CELL, 0);
}

void IceMesh::ComputeElementsCartesianSize()
{
    // the Cartesian size of every triangle is 0.0
    grid_info[mesh::gridElemType::Node]->GetCartesianSize() = ice_mesh->CreateTag("node cart size ", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
    ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"node cart size ", "nosave");

    // compute the Cartesian length of all edges
    grid_info[mesh::gridElemType::Edge]->GetCartesianSize() = ice_mesh->CreateTag("edge cart size ", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 1);
    ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"edge cart size ", "nosave");

    for (auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
    {
        if (edgeit->GetStatus() != Element::Ghost)
        {
            // get adjacent nodes
            ElementArray<INMOST::Node> adj_nodes = edgeit->getNodes();

            // compute edge coords
            std::vector<double> edge_coords = 
            {
                adj_nodes[1]->RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[0] - adj_nodes[0]->RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[0],
                adj_nodes[1]->RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[1] - adj_nodes[0]->RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[1], 
                adj_nodes[1]->RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[2] - adj_nodes[0]->RealArray(grid_info[mesh::gridElemType::Node]->coords[coord::coordType::cart])[2]
            };

            edgeit->Real(grid_info[mesh::gridElemType::Edge]->GetCartesianSize()) = L2_norm_vec(edge_coords);
        }
    }
    // exchange data
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Edge]->GetCartesianSize(), INMOST::FACE, 0);

    // compute Cartesian area of all triangles
    grid_info[mesh::gridElemType::Trian]->GetCartesianSize() = ice_mesh->CreateTag("trian cart size ", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
    //ice_mesh->SetFileOption((std::string)"Tag:" + (std::string)"trian cart size ", "nosave");
    
    for (auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            // get adjacent edges
            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

            // get edge lengths
            std::vector<double> edge_lengths = 
            {
                adj_edges[0]->Real(grid_info[mesh::gridElemType::Edge]->GetCartesianSize()),
                adj_edges[1]->Real(grid_info[mesh::gridElemType::Edge]->GetCartesianSize()),
                adj_edges[2]->Real(grid_info[mesh::gridElemType::Edge]->GetCartesianSize())
            };
            
            // compute semi-perimeter
            double p = 0.5*(edge_lengths[0] + edge_lengths[1] + edge_lengths[2]);
            trianit->Real(grid_info[mesh::gridElemType::Trian]->GetCartesianSize()) = std::sqrt(p*(p - edge_lengths[0])*(p - edge_lengths[1])*(p - edge_lengths[2]));
        }
    }
    // exchange data
    ice_mesh->ExchangeData(grid_info[mesh::gridElemType::Trian]->GetCartesianSize(), INMOST::CELL, 0);
}

void IceMesh::AssembleBasisData()
{
    AssembleGeoElementBasis();
    AssembleCartesianElementBasis();
    AssembleGeoToElementTransitionMatricies();
    AssembleElementToElementTransitionMatricies();
    ComputeNodeCoordsInTrianBasis();
    ComputeIfXedgeBasisIsNormalToTrian();
    ComputeElementsCartesianSize();
}

}