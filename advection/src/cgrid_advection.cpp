#include "advection.hpp"

using namespace INMOST;
using namespace SIMUG;

void CgridAdvectionSolver::Evaluate(velocity_tag vel_tag, scalar_tag scal_tag)
{
    // get global id tag for node
    INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

    // get global id tag for trian
    INMOST::Tag trian_id_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->id;

    // get cart coord tag for node
    INMOST::Tag node_cart_coords_tag =  mesh->GetGridInfo(mesh::gridElemType::Node)->coords[coord::coordType::cart];

    // get trian basis tag for trian
    std::vector<INMOST::Tag> trian_basis_tag =  mesh->GetGridInfo(mesh::gridElemType::Trian)->cart_basis;

    // get transition matrix from edge geo to edge cart basis
    INMOST::Tag edge_geo_to_elem_tag = mesh->GetGridInfo(mesh::gridElemType::Edge)->trans_matr_from_geo_to_elem;

    // get transition matrix from edge to trian basis
    std::vector<INMOST::Tag> edge_to_trian_matr_tags = mesh->GetGridInfo(mesh::gridElemType::Edge)->GetTransMatrToTrian();

    // calculate new scalar on every triangle
    for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
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

            // calculate coefficients of edge velocities basis functions on current triangle
            std::vector<double> coeffs0 = solve_linear_system(std::vector<std::vector<double>>{{zero_node_coords_in_trian_basis[0], zero_node_coords_in_trian_basis[1], 1.0},
                                                                                               {first_node_coords_in_trian_basis[0], first_node_coords_in_trian_basis[1], 1.0},
                                                                                               {second_node_coords_in_trian_basis[0], second_node_coords_in_trian_basis[1], 1.0}},
                                                               std::vector<double>{-1.0, 1.0, 1.0});

            std::vector<double> coeffs1 = solve_linear_system(std::vector<std::vector<double>>{{zero_node_coords_in_trian_basis[0], zero_node_coords_in_trian_basis[1], 1.0},
                                                                                               {first_node_coords_in_trian_basis[0], first_node_coords_in_trian_basis[1], 1.0},
                                                                                               {second_node_coords_in_trian_basis[0], second_node_coords_in_trian_basis[1], 1.0}},
                                                               std::vector<double>{1.0, -1.0, 1.0});
            
            std::vector<double> coeffs2 = solve_linear_system(std::vector<std::vector<double>>{{zero_node_coords_in_trian_basis[0], zero_node_coords_in_trian_basis[1], 1.0},
                                                                                               {first_node_coords_in_trian_basis[0], first_node_coords_in_trian_basis[1], 1.0},
                                                                                               {second_node_coords_in_trian_basis[0], second_node_coords_in_trian_basis[1], 1.0}},
                                                            std::vector<double>{1.0, 1.0, -1.0});

            // calculate gradient of local basis functions
            std::vector<double> grad0 = {coeffs0[0], coeffs0[1]};
            std::vector<double> grad1 = {coeffs1[0], coeffs1[1]};
            std::vector<double> grad2 = {coeffs2[0], coeffs2[1]};

            // get edges of current trian
            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

            // get transition matricies from geo to elem for edges
            std::vector<std::vector<double>> matr_geo_to_edge0 = {
                        {adj_edges[0].RealArray(edge_geo_to_elem_tag)[0], adj_edges[0].RealArray(edge_geo_to_elem_tag)[1]},
                        {adj_edges[0].RealArray(edge_geo_to_elem_tag)[2], adj_edges[0].RealArray(edge_geo_to_elem_tag)[3]}
                                                     };
            
            std::vector<std::vector<double>> matr_geo_to_edge1 = {
                        {adj_edges[1].RealArray(edge_geo_to_elem_tag)[0], adj_edges[1].RealArray(edge_geo_to_elem_tag)[1]},
                        {adj_edges[1].RealArray(edge_geo_to_elem_tag)[2], adj_edges[1].RealArray(edge_geo_to_elem_tag)[3]}
                                                     };
            
            std::vector<std::vector<double>> matr_geo_to_edge2 = {
                        {adj_edges[2].RealArray(edge_geo_to_elem_tag)[0], adj_edges[2].RealArray(edge_geo_to_elem_tag)[1]},
                        {adj_edges[2].RealArray(edge_geo_to_elem_tag)[2], adj_edges[2].RealArray(edge_geo_to_elem_tag)[3]}
                                                     };

            // compute edge velocities in edge basis
            std::vector<double> edge_vel0 = matr_geo_to_edge0*std::vector<double>{adj_edges[0].RealArray(vel_tag)[0], adj_edges[0].RealArray(vel_tag)[1]};
            std::vector<double> edge_vel1 = matr_geo_to_edge1*std::vector<double>{adj_edges[1].RealArray(vel_tag)[0], adj_edges[1].RealArray(vel_tag)[1]};
            std::vector<double> edge_vel2 = matr_geo_to_edge2*std::vector<double>{adj_edges[2].RealArray(vel_tag)[0], adj_edges[2].RealArray(vel_tag)[1]};

            // find current trian num in edge numeration
            std::vector<size_t> curr_tr_num_in_edge_numerations(3);
            curr_tr_num_in_edge_numerations[0] = (adj_edges[0].getCells()[0].Integer(trian_id_tag) == trianit->Integer(trian_id_tag))?0:1;
            curr_tr_num_in_edge_numerations[1] = (adj_edges[1].getCells()[0].Integer(trian_id_tag) == trianit->Integer(trian_id_tag))?0:1;
            curr_tr_num_in_edge_numerations[2] = (adj_edges[2].getCells()[0].Integer(trian_id_tag) == trianit->Integer(trian_id_tag))?0:1;

            // get transition matricies from edge to trian
            std::vector<std::vector<double>> matr_edge_to_trian0 = {
                        {adj_edges[0].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[0]])[0], adj_edges[0].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[0]])[1]},
                        {adj_edges[0].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[0]])[2], adj_edges[0].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[0]])[3]}
                                                     };
            
            std::vector<std::vector<double>> matr_edge_to_trian1 = {
                        {adj_edges[1].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[1]])[0], adj_edges[1].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[1]])[1]},
                        {adj_edges[1].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[1]])[2], adj_edges[1].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[1]])[3]}
                                                     };
            
            std::vector<std::vector<double>> matr_edge_to_trian2 = {
                        {adj_edges[2].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[2]])[0], adj_edges[2].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[2]])[1]},
                        {adj_edges[2].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[2]])[2], adj_edges[2].RealArray(edge_to_trian_matr_tags[curr_tr_num_in_edge_numerations[2]])[3]}
                                                     };

            // compute edge velocities in trian basis
            std::vector<double> edge_vel_in_trian_basis0 = matr_edge_to_trian0*std::vector<double>{adj_edges[0].RealArray(vel_tag)[0], adj_edges[0].RealArray(vel_tag)[1]};
            std::vector<double> edge_vel_in_trian_basis1 = matr_edge_to_trian1*std::vector<double>{adj_edges[1].RealArray(vel_tag)[0], adj_edges[1].RealArray(vel_tag)[1]};
            std::vector<double> edge_vel_in_trian_basis2 = matr_edge_to_trian2*std::vector<double>{adj_edges[2].RealArray(vel_tag)[0], adj_edges[2].RealArray(vel_tag)[1]};

            std::vector<std::vector<double>> all_edge_vel_in_tr_basis = 
            {
                edge_vel_in_trian_basis0,
                edge_vel_in_trian_basis1,
                edge_vel_in_trian_basis2
            };

            // find oppesite edge for every node
            std::vector<size_t> oposite_edge_num_for_node(3);

            for (size_t edgenum = 0; edgenum < 3; ++edgenum)
            {
                ElementArray<INMOST::Node> adj_nodes_for_edge = adj_edges[edgenum].getNodes();
                for (size_t nodenum = 0; nodenum < 3; ++nodenum)
                {
                    if ((adj_nodes[nodenum].Integer(node_id_tag) != adj_nodes_for_edge[0].Integer(node_id_tag) ) and
                        (adj_nodes[nodenum].Integer(node_id_tag) != adj_nodes_for_edge[1].Integer(node_id_tag) ))
                    {
                        oposite_edge_num_for_node[nodenum] = edgenum;
                        break;
                    }
                }
            }

            // compute velocity divergence
            std::vector<double> vel_div = 
            {
                grad0*all_edge_vel_in_tr_basis[oposite_edge_num_for_node[0]],
                grad1*all_edge_vel_in_tr_basis[oposite_edge_num_for_node[1]],
                grad2*all_edge_vel_in_tr_basis[oposite_edge_num_for_node[2]]
            };

            // finally calculate new scalar
            double old_scalar = trianit->Real(scal_tag);
            double temp_scalar = old_scalar - (time_step/2.0)*old_scalar*(vel_div[0] + vel_div[1] + vel_div[2]);
            trianit->Real(scal_tag) = old_scalar - time_step*temp_scalar*(vel_div[0] + vel_div[1] + vel_div[2]);
        }
    }
    
    // exchange scalar
    mesh->GetMesh()->ExchangeData(scal_tag, INMOST::CELL, 0);
    
    BARRIER
};