#include "momentum.hpp"

using namespace INMOST;

namespace SIMUG
{
    CgridMomentumSolver::CgridMomentumSolver(SIMUG::IceMesh* mesh_,
                                             double time_step_,
                                             velocity_tag vel_tag_,
                                             scalar_tag mass_tag_,
                                             scalar_tag conc_tag_,
                                             scalar_tag thick_tag_,
                                             SIMUG::dyn::timeScheme mom_time_scheme_,
                                             SIMUG::dyn::pressParam mom_press_param_,
                                             SIMUG::dyn::bcType mom_bc_type_,
                                             const std::vector<double>& real_params_,
                                             const std::vector<int>& integer_params_):
            MomentumSolver(mesh_,
                           time_step_,
                           vel_tag_,
                           mass_tag_,
                           conc_tag_,
                           thick_tag_,
                           mom_time_scheme_,
                           SIMUG::dyn::spaceScheme::stabCR,
                           mom_press_param_,
                           mom_bc_type_),
            real_params(real_params_),
            integer_params(integer_params_)
    {
        // initialize timer and logger
        SIMUG::Logger mom_log(std::cout);
        SIMUG::Timer mom_timer;
        double duration;

        // assemble edge-based mass matrix
        mass_matrix_entry_tag = mesh->GetMesh()->CreateTag("mass_matrix_entry_tag", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 1);
        //mesh->GetMesh()->SetFileOption("Tag:mass_matrix_entry_tag", "nosave");
        mom_timer.Launch();
        AssembleMassMatrix(mass_matrix_entry_tag);
        mom_timer.Stop();
        duration = mom_timer.GetMaxTime();
        mom_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("Assembling of edge-based mass matrix: OK (" + std::to_string(duration) + " ms)\n");
        }
        BARRIER

        // compute opposite edge for every node of triangle
        opposite_edge_for_node_tags = mesh->GetMesh()->CreateTag("opposite_edge_for_node_tags", INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 3);
        mesh->GetMesh()->SetFileOption("Tag:opposite_edge_for_node_tags", "nosave");
        mom_timer.Launch();
        ComputeOppositeEdges(opposite_edge_for_node_tags);
        mom_timer.Stop();
        duration = mom_timer.GetMaxTime();
        mom_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("Computation of opposite edge number for nodes: OK (" + std::to_string(duration) + " ms)\n");
        }
        BARRIER

        // compute opposite node for every edge of triangle
        opposite_node_for_edge_tags = mesh->GetMesh()->CreateTag("opposite_node_for_edge_mom_tags", INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 3);
        mesh->GetMesh()->SetFileOption("Tag:opposite_node_for_edge_mom_tags", "nosave");
        mom_timer.Launch();
        ComputeOppositeNodes(opposite_node_for_edge_tags);
        mom_timer.Stop();
        duration = mom_timer.GetMaxTime();
        mom_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("Computation of opposite node number for edges: OK (" + std::to_string(duration) + " ms)\n");
        }
        BARRIER

        // compute edge basis in trian coords
        outward_edge_basis_in_trian_coords_tags = mesh->GetMesh()->CreateTag("edge_basis_in_trian_coords_tags", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 12);
        mesh->GetMesh()->SetFileOption("Tag:outward_edge_basis_in_trian_coords_tags", "nosave");
        mom_timer.Launch();
        ComputeEdgeBasisInTrianCoords(outward_edge_basis_in_trian_coords_tags);
        mom_timer.Stop();
        duration = mom_timer.GetMaxTime();
        mom_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("Computation of outward edge basis in trian coords: OK (" + std::to_string(duration) + " ms)\n");
        }
        BARRIER

        // compute trian height to edge
        trian_height_to_edge_tags = mesh->GetMesh()->CreateTag("trian_height_to_edge_tags", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        mesh->GetMesh()->SetFileOption("Tag:trian_height_to_edge_tags", "nosave");
        mom_timer.Launch();
        ComputeTrianHeights(trian_height_to_edge_tags);
        mom_timer.Stop();
        duration = mom_timer.GetMaxTime();
        mom_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("Computation of trian heights: OK (" + std::to_string(duration) + " ms)\n");
        }
        BARRIER

        // compute gradients of basis functions
        grad_basis_func_tags = mesh->GetMesh()->CreateTag("grad_basis_func_tags", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 6);
        mesh->GetMesh()->SetFileOption("Tag:grad_basis_func_tags", "nosave");
        mom_timer.Launch();
        ComputeGradientsBasisFuncs(grad_basis_func_tags);
        mom_timer.Stop();
        duration = mom_timer.GetMaxTime();
        mom_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("Computation of gradients of basis functions: OK (" + std::to_string(duration) + " ms)\n");
        }
        BARRIER
    }

    void CgridMomentumSolver::AssembleMassMatrix(INMOST::Tag mass_matrix_tag)
    {
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            double trian_area = trianit->Real(mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize());
            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();
                
            for (int edge_ind = 0; edge_ind < 3; ++edge_ind)
            {
                adj_edges[edge_ind]->Real(mass_matrix_tag) = adj_edges[edge_ind]->Real(mass_matrix_tag) + trian_area/3.0;
            }
        }
        mesh->GetMesh()->ExchangeData(mass_matrix_tag, INMOST::FACE, 0);
        BARRIER
    }

    void CgridMomentumSolver::ComputeOppositeEdges(INMOST::Tag op_edge_tags)
    {
        // tag for node id
        INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                // get nodes of trian
                ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

                // get edges of trian
                ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

                // iterate over nodes
                for (int node_num = 0; node_num < 3; ++node_num)
                {
                    int opposite_edge_num_for_node = -1;

                    // iterate over nodes and find opposite edge
                    for (int ed_num = 0; ed_num < 3; ++ed_num)
                    {
                        if ((adj_nodes[node_num].Integer(node_id_tag) != adj_edges[ed_num].getNodes()[0].Integer(node_id_tag)) &&
                            (adj_nodes[node_num].Integer(node_id_tag) != adj_edges[ed_num].getNodes()[1].Integer(node_id_tag)))
                        {
                            opposite_edge_num_for_node = ed_num;
                            break;
                        }
                    }

                    if (opposite_edge_num_for_node == -1)
                    {
                        SIMUG_ERR("cant find opposite edge num for node");
                    }

                    // store opposite node num
                    trianit->IntegerArray(op_edge_tags)[node_num] = opposite_edge_num_for_node;
                }
            }
        }
        mesh->GetMesh()->ExchangeData(op_edge_tags, INMOST::CELL, 0);
        BARRIER
    }

    // computation of opposite node number for every edge of triangle
    void CgridMomentumSolver::ComputeOppositeNodes(INMOST::Tag op_node_tags)
    {
        // tag for node id
        INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                // get nodes of trian
                ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

                // get edges of trian
                ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

                // iterate over edges
                for (int edge_num = 0; edge_num < 3; ++edge_num)
                {
                    int opposite_node_num_for_edge = -1;

                    // iterate over nodes and find opposite edge
                    for (int node_num = 0; node_num < 3; ++node_num)
                    {
                        if ((adj_nodes[node_num].Integer(node_id_tag) != adj_edges[edge_num].getNodes()[0].Integer(node_id_tag)) &&
                            (adj_nodes[node_num].Integer(node_id_tag) != adj_edges[edge_num].getNodes()[1].Integer(node_id_tag)))
                        {
                            opposite_node_num_for_edge = node_num;
                            break;
                        }
                    }

                    if (opposite_node_num_for_edge == -1)
                    {
                        SIMUG_ERR("cant find opposite node num for edge");
                    }

                    // store opposite node num
                    trianit->IntegerArray(op_node_tags)[edge_num] = opposite_node_num_for_edge;
                }
            }
        }
        mesh->GetMesh()->ExchangeData(op_node_tags, INMOST::CELL, 0);
        BARRIER
    }

    void CgridMomentumSolver::ComputeEdgeBasisInTrianCoords(INMOST::Tag edge_basis_in_trian_coords_tags)
    {
        std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();
        INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                ElementArray<INMOST::Face> adj_edges = trianit->getFaces();
                ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

                // get node coords in trian basis
                std::vector<std::vector<double>> node_coords;
                
                for(int i = 0; i < 3; ++i)
                {
                    node_coords.push_back
                    (
                        std::vector<double>
                        {
                            trianit->RealArray(node_coords_in_trian_basis_tags[i])[0],
                            trianit->RealArray(node_coords_in_trian_basis_tags[i])[1]
                        }
                    );
                }

                // ### compute outvard normal and tangential vectors ###

                std::vector<std::vector<double>> tangentials;
                std::vector<std::vector<double>> out_normals;

                // compute tangentials
                tangentials.push_back
                (
                    node_coords[1] - node_coords[0]
                );

                tangentials.push_back
                (
                    node_coords[2] - node_coords[1]
                );

                tangentials.push_back
                (
                    node_coords[0] - node_coords[2]
                );

                // normalize tangentials
                for (int i = 0; i < 3; ++i)
                {
                    tangentials[i] = (1.0/L2_norm_vec(tangentials[i]))*tangentials[i];
                }

                // compute outward normals
                for (int i = 0; i < 3; ++i)
                {
                    out_normals.push_back
                    (
                        std::vector<double>
                        {
                            -tangentials[i][1],
                            tangentials[i][0]
                        }
                    );
                }

                // normalize outward normals
                for (int i = 0; i < 3; ++i)
                {
                    out_normals[i] = (1.0/L2_norm_vec(out_normals[i]))*out_normals[i];
                }

                // check if normals outward
                if (out_normals[0]*(node_coords[0] - node_coords[2]) <= 0.0)
                {
                    out_normals[0] = (-1.0)*out_normals[0];
                }

                if (out_normals[1]*(node_coords[2] - node_coords[0]) <= 0.0)
                {
                    out_normals[1] = (-1.0)*out_normals[1];
                }

                if (out_normals[2]*(node_coords[0] - node_coords[1]) <= 0.0)
                {
                    out_normals[2] = (-1.0)*out_normals[2];
                }

                // check if tangentials on the right 
                for (int i = 0; i < 3; ++i)
                {
                    if ((std::vector<double>{tangentials[i][0], tangentials[i][1], 0.0}%
                         std::vector<double>{out_normals[i][0], out_normals[i][1], 0.0})[2] <= 0.0)
                    {
                        tangentials[i] = (-1.0)*tangentials[i];
                    }
                }

                // store data
                for (int ed_num = 0; ed_num < 3; ++ed_num)
                {
                    ElementArray<Node> nodes_of_current_edge = adj_edges[ed_num].getNodes(); 
                    
                    int current_edge_num = 0;

                    if ((nodes_of_current_edge[0]->Integer(node_id_tag) != adj_nodes[2]->Integer(node_id_tag)) &&
                        (nodes_of_current_edge[1]->Integer(node_id_tag) != adj_nodes[2]->Integer(node_id_tag)))
                    {
                        current_edge_num = 0;
                    }
                    else if ((nodes_of_current_edge[0]->Integer(node_id_tag) != adj_nodes[0]->Integer(node_id_tag)) &&
                             (nodes_of_current_edge[1]->Integer(node_id_tag) != adj_nodes[0]->Integer(node_id_tag)))
                    {
                        current_edge_num = 1;
                    }
                    else
                    {
                        current_edge_num = 2;
                    }

                    trianit->RealArray(edge_basis_in_trian_coords_tags)[current_edge_num*4 + 0] = tangentials[current_edge_num][0];
                    trianit->RealArray(edge_basis_in_trian_coords_tags)[current_edge_num*4 + 1] = tangentials[current_edge_num][1];
                    trianit->RealArray(edge_basis_in_trian_coords_tags)[current_edge_num*4 + 2] = out_normals[current_edge_num][0];
                    trianit->RealArray(edge_basis_in_trian_coords_tags)[current_edge_num*4 + 3] = out_normals[current_edge_num][1];
                }
            }
        }
        mesh->GetMesh()->ExchangeData(edge_basis_in_trian_coords_tags, INMOST::CELL, 0);
        BARRIER
    }

    void CgridMomentumSolver::ComputeTrianHeights(INMOST::Tag trian_hight_tags)
    {
        std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();
        INMOST::Tag trian_area_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize();
        INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                ElementArray<INMOST::Face> adj_edges = trianit->getFaces();
                ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

                // get trian area 
                double trian_area = trianit->Real(trian_area_tag);

                // get node coords in trian basis
                std::vector<std::vector<double>> node_coords;
                
                for(int i = 0; i < 3; ++i)
                {
                    node_coords.push_back
                    (
                        std::vector<double>
                        {
                            trianit->RealArray(node_coords_in_trian_basis_tags[i])[0],
                            trianit->RealArray(node_coords_in_trian_basis_tags[i])[1]
                        }
                    );
                }

                // compute edge sizes
                std::vector<double> edge_sizes;

                edge_sizes.push_back(L2_norm_vec(node_coords[1] - node_coords[0]));
                edge_sizes.push_back(L2_norm_vec(node_coords[2] - node_coords[1]));
                edge_sizes.push_back(L2_norm_vec(node_coords[0] - node_coords[2]));

                // compute trian heights
                std::vector<double> trian_heights;
                trian_heights.push_back(2.0*trian_area/edge_sizes[0]);
                trian_heights.push_back(2.0*trian_area/edge_sizes[1]);
                trian_heights.push_back(2.0*trian_area/edge_sizes[2]);

                /// store data
                for (int ed_num = 0; ed_num < 3; ++ed_num)
                {
                    ElementArray<Node> nodes_of_current_edge = adj_edges[ed_num].getNodes(); 
                    int current_edge_num = 0;

                    if ((nodes_of_current_edge[0]->Integer(node_id_tag) != adj_nodes[2]->Integer(node_id_tag)) &&
                        (nodes_of_current_edge[1]->Integer(node_id_tag) != adj_nodes[2]->Integer(node_id_tag)))
                    {
                        current_edge_num = 0;
                    }
                    else if ((nodes_of_current_edge[0]->Integer(node_id_tag) != adj_nodes[0]->Integer(node_id_tag)) &&
                             (nodes_of_current_edge[1]->Integer(node_id_tag) != adj_nodes[0]->Integer(node_id_tag)))
                    {
                        current_edge_num = 1;
                    }
                    else
                    {
                        current_edge_num = 2;
                    }

                    trianit->RealArray(trian_hight_tags)[current_edge_num] = trian_heights[current_edge_num];
                }
            }
        }
        mesh->GetMesh()->ExchangeData(trian_hight_tags, INMOST::CELL, 0);
        BARRIER
    }

    void CgridMomentumSolver::ComputeGradientsBasisFuncs(INMOST::Tag trian_hight_tags)
    {
        std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();

        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

                // get coordinates of nodes in trian basis
                std::vector<std::vector<double>> node_coords_tr_basis =
                {
                    {trianit->RealArray(node_coords_in_trian_basis_tags[0])[0], trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]},
                    {trianit->RealArray(node_coords_in_trian_basis_tags[1])[0], trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]},
                    {trianit->RealArray(node_coords_in_trian_basis_tags[2])[0], trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]}
                };


                for (int ed_num = 0; ed_num < 3; ++ed_num)
                {
                    // get nodal values of current basis function
                    std::vector<double> basis_func_nodal_values = {1.0, 1.0, 1.0};
                    basis_func_nodal_values[trianit->IntegerArray(opposite_node_for_edge_tags)[ed_num]] = -1.0;

                    // compute coefficients of current basis function
                    std::vector<double> coeffs = solve_linear_system
                    (
                        std::vector<std::vector<double>>
                        {
                            {node_coords_tr_basis[0][0], node_coords_tr_basis[0][1], 1.0},
                            {node_coords_tr_basis[1][0], node_coords_tr_basis[1][1], 1.0},
                            {node_coords_tr_basis[2][0], node_coords_tr_basis[2][1], 1.0},
                        },
                        std::vector<double>
                        {basis_func_nodal_values[0], basis_func_nodal_values[1], basis_func_nodal_values[2]}
                    );

                    trianit->RealArray(grad_basis_func_tags)[ed_num*2 + 0] = coeffs[0];
                    trianit->RealArray(grad_basis_func_tags)[ed_num*2 + 1] = coeffs[1];
                }
            }
        }
        mesh->GetMesh()->ExchangeData(grad_basis_func_tags, INMOST::CELL, 0);

        BARRIER
    }

    Cgrid_mEVP_Solver::Cgrid_mEVP_Solver(SIMUG::IceMesh* mesh_,
                                         double time_step_,
                                         velocity_tag vel_tag_,
                                         scalar_tag mass_tag_,
                                         scalar_tag conc_tag_,
                                         scalar_tag thick_tag_,
                                         SIMUG::dyn::pressParam mom_press_param_,
                                         SIMUG::dyn::bcType mom_bc_type_,
                                         const std::vector<double>& real_params_,
                                         const std::vector<int>& integer_params_) :
            CgridMomentumSolver(mesh_,
                                time_step_,
                                vel_tag_,
                                mass_tag_,
                                conc_tag_,
                                thick_tag_,
                                SIMUG::dyn::timeScheme::mEVP,
                                mom_press_param_,
                                mom_bc_type_,
                                real_params_,
                                integer_params_)
    {
        // initialize timer and logger
        SIMUG::Logger mom_log(std::cout);
        SIMUG::Timer mom_timer;

        // create mandatory tags
        P_tag = mesh->GetMesh()->CreateTag("P_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
        vareps_tag = mesh->GetMesh()->CreateTag("vareps_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        delta_tag = mesh->GetMesh()->CreateTag("delta_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
        force_tags = mesh->GetMesh()->CreateTag("force_tags", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 2);
        edge_stab_tags = mesh->GetMesh()->CreateTag("edge_stab_tags", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 2);
        
        // mute tags
        mesh->GetMesh()->SetFileOption("Tag:P_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:vareps_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:delta_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:force_tags", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:edge_stab_tags", "nosave");

        ComputeVelocity();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("=====================================================================\n");
        }
        
        BARRIER
    }


    void Cgrid_mEVP_Solver::ComputeP()
    {
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                ElementArray<Node> nodes = trianit->getNodes();

                double h = trianit->Real(thick_tag);

                double a = trianit->Real(conc_tag);

                double p_str = IceConsts::pstr;
                double C = IceConsts::C;

                trianit->Real(P_tag) = p_str*h*std::exp(-C*(1.0 - a));
            }
        }

        mesh->GetMesh()->ExchangeData(P_tag, CELL, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::ComputeVarepsilonDelta(INMOST::Tag vel_tag)
    {
        double e = IceConsts::e;
        std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();
        std::vector<INMOST::Tag> is_normal_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetIsXedgeBasisIsNormal();

        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                ElementArray<Face> edges = trianit->getFaces();

                std::vector<std::vector<double>> edge_velocities;
                for (int i = 0; i < 3; ++i)
                {
                    edge_velocities.push_back
                    (
                        mesh->VecTransitionToElemBasis
                        (
                            std::vector<double>
                            {
                                edges[i]->RealArray(vel_tag)[0],
                                edges[i]->RealArray(vel_tag)[1]
                            },
                            edges[i]
                        )
                    );
                }

                // move edge velocity components to triangular basis
                std::vector<std::vector<double>> edge_velocities_trian_basis;

                for (int i = 0; i < 3; ++i)
                {
                    edge_velocities_trian_basis.push_back
                    (
                        mesh->VecTransition(edge_velocities[i], edges[i], trianit->getCells()[0])
                    );
                }


                std::vector<double> ue_T = std::vector<double>
                {
                    edge_velocities_trian_basis[0][0],
                    edge_velocities_trian_basis[1][0],
                    edge_velocities_trian_basis[2][0]
                };

                std::vector<double> ve_T = std::vector<double>
                {
                    edge_velocities_trian_basis[0][1],
                    edge_velocities_trian_basis[1][1],
                    edge_velocities_trian_basis[2][1]
                };

                // compute nodal values of velocities
                std::vector<double> un_T =
                {
                    ue_T[0] + ue_T[1] + ue_T[2] - 2*ue_T[trianit->IntegerArray(opposite_edge_for_node_tags)[0]],
                    ue_T[0] + ue_T[1] + ue_T[2] - 2*ue_T[trianit->IntegerArray(opposite_edge_for_node_tags)[1]],
                    ue_T[0] + ue_T[1] + ue_T[2] - 2*ue_T[trianit->IntegerArray(opposite_edge_for_node_tags)[2]]
                };

                std::vector<double> vn_T =
                {
                    ve_T[0] + ve_T[1] + ve_T[2] - 2*ve_T[trianit->IntegerArray(opposite_edge_for_node_tags)[0]],
                    ve_T[0] + ve_T[1] + ve_T[2] - 2*ve_T[trianit->IntegerArray(opposite_edge_for_node_tags)[1]],
                    ve_T[0] + ve_T[1] + ve_T[2] - 2*ve_T[trianit->IntegerArray(opposite_edge_for_node_tags)[2]]
                };

                // get coordinates of nodes in trian basis
                std::vector<std::vector<double>> node_coords_tr_basis =
                {
                    {trianit->RealArray(node_coords_in_trian_basis_tags[0])[0], trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]},
                    {trianit->RealArray(node_coords_in_trian_basis_tags[1])[0], trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]},
                    {trianit->RealArray(node_coords_in_trian_basis_tags[2])[0], trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]}
                };

                // compute gradient of velocity in trian basis
                std::vector<double> sol_u = solve_linear_system
                (
                    std::vector<std::vector<double>>
                    {
                        {node_coords_tr_basis[0][0], node_coords_tr_basis[0][1], 1.0},
                        {node_coords_tr_basis[1][0], node_coords_tr_basis[1][1], 1.0},
                        {node_coords_tr_basis[2][0], node_coords_tr_basis[2][1], 1.0}
                    },
                    std::vector<double>
                    {un_T[0], un_T[1], un_T[2]}
                );


                std::vector<double> sol_v = solve_linear_system
                (
                    std::vector<std::vector<double>>
                    {
                        {node_coords_tr_basis[0][0], node_coords_tr_basis[0][1], 1.0},
                        {node_coords_tr_basis[1][0], node_coords_tr_basis[1][1], 1.0},
                        {node_coords_tr_basis[2][0], node_coords_tr_basis[2][1], 1.0},
                    },
                    std::vector<double>
                    {vn_T[0], vn_T[1], vn_T[2]}
                );

                double du_dx = sol_u[0];
                double du_dy = sol_u[1];

                double dv_dx = sol_v[0];
                double dv_dy = sol_v[1];


                // compute strain rate tensor
                double eps11 = du_dx;
                double eps22 = dv_dy;
                double eps12 = 0.5*(du_dy + dv_dx);

                // compute delta
                double del_min = IceConsts::delmin;

                double delta = sqrt( (eps11*eps11 + eps22*eps22)*(1.0 + 1.0/(e*e)) +
                                      eps12*eps12*(4.0/(e*e)) +
                                      eps11*eps22*2.0*(1.0 - 1.0/(e*e)) + del_min*del_min
                                   );

                // store strain rate and delta
                trianit->RealArray(vareps_tag)[0] = eps11 + eps22;
                trianit->RealArray(vareps_tag)[1] = eps11 - eps22;
                trianit->RealArray(vareps_tag)[2] = eps12;
                trianit->Real(delta_tag) = delta;
            }
        }
        mesh->GetMesh()->ExchangeData(vareps_tag, CELL, 0);
        mesh->GetMesh()->ExchangeData(delta_tag, CELL, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::AssembleForceVector(INMOST::Tag sig_tag)
    {
        INMOST::Tag trian_area_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize();

        // zeroing of force vector
        for (auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            edgeit->RealArray(force_tags)[0] = 0.0;
            edgeit->RealArray(force_tags)[1] = 0.0;
        }

        mesh->GetMesh()->ExchangeData(force_tags, FACE, 0);
        BARRIER

        // recalculate force vector
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                // get area of trian
                double trian_area = trianit->Real(trian_area_tag);

                ElementArray<Face> adj_edges = trianit->getFaces();

                for (int ed_num = 0; ed_num < 3; ++ed_num)
                {
                    // get gradient of basis function in edge basis
                    std::vector<double> grad_basis = 
                    {
                        mesh->VecTransition
                        (
                            std::vector<double>
                            {
                                trianit->RealArray(grad_basis_func_tags)[ed_num*2 + 0],
                                trianit->RealArray(grad_basis_func_tags)[ed_num*2 + 1]
                            },
                            trianit->getCells()[0],
                            adj_edges[ed_num]
                        )
                    };

                    // move sigma components to edge basis
                    double sigma1 = trianit->RealArray(sig_tag)[0];
                    double sigma2 = trianit->RealArray(sig_tag)[1];
                    double sigma12 = trianit->RealArray(sig_tag)[2];

                    double sigma11 = 0.5*(sigma1 + sigma2);
                    double sigma22 = 0.5*(sigma1 - sigma2);

                    std::vector<std::vector<double>> sig_edge = 
                    {
                        mesh->TensTransition
                        (
                            std::vector<std::vector<double>>
                            {
                                {sigma11, sigma12},
                                {sigma12, sigma22}
                            },
                            trianit->getCells()[0],
                            adj_edges[ed_num]
                        )
                    };

                    // add force from trian to edge
                    adj_edges[ed_num].RealArray(force_tags)[0] = adj_edges[ed_num].RealArray(force_tags)[0] +
                                                                 trian_area*(sig_edge[0][0]*grad_basis[0] + sig_edge[0][1]*grad_basis[1])
                                                                 /adj_edges[ed_num]->Real(mass_matrix_entry_tag);

                    adj_edges[ed_num].RealArray(force_tags)[1] = adj_edges[ed_num].RealArray(force_tags)[1] +
                                                                 trian_area*(sig_edge[1][0]*grad_basis[0] + sig_edge[1][1]*grad_basis[1])
                                                                 /adj_edges[ed_num]->Real(mass_matrix_entry_tag);
                }
            }
        }

        mesh->GetMesh()->ExchangeData(force_tags, FACE, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::UpdateSigmaMevp(INMOST::Tag sig_tag)
    {
        double alpha = real_params[0];
        double e = IceConsts::e;

        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                trianit->RealArray(sigma_tag)[0] = (1.0 - 1.0/alpha)*trianit->RealArray(sigma_tag)[0] + 
                                                   (1.0/alpha)*(trianit->Real(P_tag)/trianit->Real(delta_tag))*
                                                   (trianit->RealArray(vareps_tag)[0] - trianit->Real(delta_tag));
                
                trianit->RealArray(sigma_tag)[1] = (1.0 - 1.0/alpha)*trianit->RealArray(sigma_tag)[1] + 
                                                   (1.0/alpha)*(trianit->Real(P_tag)/(trianit->Real(delta_tag)*e*e))*
                                                   trianit->RealArray(vareps_tag)[1];
                
                trianit->RealArray(sigma_tag)[2] = (1.0 - 1.0/alpha)*trianit->RealArray(sigma_tag)[2] + 
                                                   (1.0/alpha)*(trianit->Real(P_tag)/(trianit->Real(delta_tag)*e*e))*
                                                   trianit->RealArray(vareps_tag)[2];
            }
        }

        mesh->GetMesh()->ExchangeData(sigma_tag, CELL, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::ComputeEdgeStabilization(INMOST::Tag vel_tag)
    {
        double alpha_stab = real_params[2];

        INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;
        INMOST::Tag edge_id_tag = mesh->GetGridInfo(mesh::gridElemType::Edge)->id;
        INMOST::Tag trian_area_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize();

        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if ((edgeit->GetStatus() != Element::Ghost) &&
                (edgeit->Integer(mesh->GetGridInfo(mesh::gridElemType::Edge)->is_bnd) == 0))
            {
                // zeroing previous stabilization
                edgeit->RealArray(edge_stab_tags)[0] = 0.0;
                edgeit->RealArray(edge_stab_tags)[1] = 0.0;

                ElementArray<Cell> adj_trians = edgeit->getCells();

                // move components of the first trian edge velocities to trian basis
                ElementArray<Face> first_trian_edges = adj_trians[0].getFaces();
                std::vector<std::vector<double>> first_trian_edge_vel;

                for (int i = 0; i < 3; ++i)
                {
                    first_trian_edge_vel.push_back
                    (
                        mesh->VecTransition
                        (
                            std::vector<double>
                            {
                                first_trian_edges[i]->RealArray(vel_tag)[0],
                                first_trian_edges[i]->RealArray(vel_tag)[1]
                            },
                            first_trian_edges[i],
                            adj_trians[0]
                        )
                    );
                }

                // move components of the second trian edge velocities to trian basis
                ElementArray<Face> second_trian_edges = adj_trians[1].getFaces();
                std::vector<std::vector<double>> second_trian_edge_vel;

                for (int i = 0; i < 3; ++i)
                {
                    second_trian_edge_vel.push_back
                    (
                        mesh->VecTransition
                        (
                            std::vector<double>
                            {
                                second_trian_edges[i]->RealArray(vel_tag)[0],
                                second_trian_edges[i]->RealArray(vel_tag)[1]
                            },
                            second_trian_edges[i],
                            adj_trians[1]
                        )
                    );
                }

                // move components for current edge basis for first and second trian
                std::vector<std::vector<double>> first_trian_vel;
                for (int i = 0; i < 3; ++i)
                {
                    first_trian_vel.push_back
                    (
                        mesh->VecTransition
                        (
                            std::vector<double>
                            {
                                first_trian_edge_vel[i][0],
                                first_trian_edge_vel[i][1]
                            },
                            adj_trians[0],
                            edgeit->getFaces()[0]
                        )
                    );
                }

                std::vector<std::vector<double>> second_trian_vel;
                for (int i = 0; i < 3; ++i)
                {
                    second_trian_vel.push_back
                    (
                        mesh->VecTransition
                        (
                            std::vector<double>
                            {
                                second_trian_edge_vel[i][0],
                                second_trian_edge_vel[i][1]
                            },
                            adj_trians[1],
                            edgeit->getFaces()[0]
                        )
                    );
                }

                // creat element array of propper edges
                ElementArray<Face> edges; 
                std::vector<std::vector<double>> velocities;

                for (int i = 0; i < 3; ++i)
                {
                    if (first_trian_edges[i]->Integer(edge_id_tag) != edgeit->Integer(edge_id_tag))
                    {
                        edges.push_back(first_trian_edges[i]);
                        velocities.push_back(first_trian_vel[i]);
                    }
                }

                for (int i = 0; i < 3; ++i)
                {
                    if (second_trian_edges[i]->Integer(edge_id_tag) != edgeit->Integer(edge_id_tag))
                    {
                        edges.push_back(second_trian_edges[i]);
                        velocities.push_back(second_trian_vel[i]);
                    }
                }

                // main cycle
                for (int i = 0; i < 4; ++i)
                {
                    for (int j = 0; j < 4; ++j)
                    {
                        if (
                                (((edges[i]->getNodes()[0])->Integer(node_id_tag) == (edges[j]->getNodes()[0])->Integer(node_id_tag)) and 
                                 ((edges[i]->getNodes()[1])->Integer(node_id_tag) != (edges[j]->getNodes()[1])->Integer(node_id_tag))) or 
                                (((edges[i]->getNodes()[0])->Integer(node_id_tag) == (edges[j]->getNodes()[1])->Integer(node_id_tag)) and 
                                 ((edges[i]->getNodes()[1])->Integer(node_id_tag) != (edges[j]->getNodes()[0])->Integer(node_id_tag)))
                           )
                        {
                            edgeit->RealArray(edge_stab_tags)[0] -= (1.0/3.0)*velocities[i][0];
                            edgeit->RealArray(edge_stab_tags)[1] -= (1.0/3.0)*velocities[i][1];
                        }
                        else
                        {
                            edgeit->RealArray(edge_stab_tags)[0] += (1.0/3.0)*velocities[i][0];
                            edgeit->RealArray(edge_stab_tags)[1] += (1.0/3.0)*velocities[i][1];
                        }
                    }
                }

                // calculate xi value on edge
                double xi_edge = (adj_trians[0]->Real(trian_area_tag)*(0.5*adj_trians[0]->Real(P_tag)/adj_trians[0]->Real(delta_tag)) +
                                  adj_trians[1]->Real(trian_area_tag)*(0.5*adj_trians[1]->Real(P_tag)/adj_trians[1]->Real(delta_tag)))/
                                  (adj_trians[0]->Real(trian_area_tag) + adj_trians[1]->Real(trian_area_tag));

                edgeit->RealArray(edge_stab_tags)[0] *=  2.0*xi_edge*alpha_stab*(1.0/edgeit->Real(mass_matrix_entry_tag));
                edgeit->RealArray(edge_stab_tags)[1] *=  2.0*xi_edge*alpha_stab*(1.0/edgeit->Real(mass_matrix_entry_tag));
            }
        }

        mesh->GetMesh()->ExchangeData(edge_stab_tags, FACE, 0);
        BARRIER
    }

/*
    void Cgrid_mEVP_Solver::ComputeVarepsilonDelta(INMOST::Tag vel_tag)
    {
        double e = IceConsts::e;
        std::vector<INMOST::Tag> is_normal_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetIsXedgeBasisIsNormal();

        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                ElementArray<Face> adj_edges = trianit->getFaces();

                // get outward edge basis in cuurent trian coords
                std::vector<std::vector<double>> tangentials = 
                {
                    std::vector<double>
                    {
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[0],
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[1]
                    },
                    std::vector<double>
                    {
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[4],
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[5]
                    },
                    std::vector<double>
                    {
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[8],
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[9]
                    }
                };

                std::vector<std::vector<double>> normals = 
                {
                    std::vector<double>
                    {
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[2],
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[3]
                    },
                    std::vector<double>
                    {
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[6],
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[7]
                    },
                    std::vector<double>
                    {
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[10],
                        trianit->RealArray(outward_edge_basis_in_trian_coords_tags)[11]
                    }
                };

                // get trian heights
                std::vector<double> heights = 
                {
                    trianit->RealArray(trian_height_to_edge_tags)[0],
                    trianit->RealArray(trian_height_to_edge_tags)[1],
                    trianit->RealArray(trian_height_to_edge_tags)[2]
                };

                // move velocity components to edge basis
                std::vector<std::vector<double>> edge_velocities;
                for (int i = 0; i < 3; ++i)
                {
                    edge_velocities.push_back
                    (
                        mesh->VecTransitionToElemBasis
                        (
                            std::vector<double>
                            {
                                adj_edges[i]->RealArray(vel_tag)[0],
                                adj_edges[i]->RealArray(vel_tag)[1]
                            },
                            adj_edges[i]
                        )
                    );
                }

                // compute velocity components in outward trian basis
                std::vector<std::vector<double>> velocities;
                for (int i = 0; i < 3; ++i)
                {
                    if (trianit->Integer(is_normal_tags[i]) == 1)
                    {
                        velocities.push_back
                        (
                            std::vector<double>
                            {
                                -edge_velocities[i][1],
                                 edge_velocities[i][0]
                            }
                        );
                    }
                    else
                    {
                        velocities.push_back
                        (
                            std::vector<double>
                            {
                                 edge_velocities[i][1],
                                -edge_velocities[i][0]
                            }
                        );
                    }
                }

                // compute strain rate tensor
                std::vector<std::vector<double>> dot_vareps(2, std::vector<double>(2));

                for (int i = 0; i < 3; ++i)
                {
                    
                    dot_vareps = dot_vareps + (1.0/heights[i])*(velocities[i][0]*vec_transpose_product(tangentials[i], normals[i]) +
                                                                velocities[i][1]*vec_transpose_product(normals[i], normals[i]) +
                                                                velocities[i][0]*vec_transpose_product(normals[i], tangentials[i])
                                                                );
                    
                }


                double eps11 = dot_vareps[0][0];
                double eps22 = dot_vareps[1][1];
                double eps12 = dot_vareps[0][1];

                // compute delta
                double delta = sqrt( (eps11*eps11 + eps22*eps22)*(1.0 + 1.0/(e*e)) +
                                      eps12*eps12*(4.0/(e*e)) +
                                      eps11*eps22*2.0*(1.0 - 1.0/(e*e))
                                   );

                // store strain rate and delta
                trianit->RealArray(vareps_tag)[0] = eps11 + eps22;
                trianit->RealArray(vareps_tag)[1] = eps11 - eps22;
                trianit->RealArray(vareps_tag)[2] = eps12;
                trianit->Real(delta_tag) = delta;
            }
        }

        mesh->GetMesh()->ExchangeData(vareps_tag, CELL, 0);
        mesh->GetMesh()->ExchangeData(delta_tag, CELL, 0);

        BARRIER
    }
*/

    void Cgrid_mEVP_Solver::ComputeVelocity()
    {
        std::cout << "calculations" << std::endl;

        // move velocity vector to edge basis
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if (edgeit->GetStatus() != Element::Ghost)
            {
                mesh->VecTransitionToElemBasis
                (
                    std::vector<double>
                    {
                        edgeit->RealArray(vel_tag)[0],
                        edgeit->RealArray(vel_tag)[1]
                    },
                    edgeit->getFaces()[0]
                );
            }
        }
        mesh->GetMesh()->ExchangeData(vel_tag, FACE, 0);
        BARRIER

        // compute P
        ComputeP();
        BARRIER

        // compute strain rate tensor and delta
        ComputeVarepsilonDelta(vel_tag);
        BARRIER 

        // compute edge stabilization
        ComputeEdgeStabilization(vel_tag);
        BARRIER

        // save mass matrix
        /*
        INMOST::Tag check_mass_matrix = mesh->GetMesh()->CreateTag("check_mass_matrix", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                trianit->RealArray(check_mass_matrix)[0] = (trianit->getFaces())[0]->Real(mass_matrix_entry_tag);
            }
        }
        mesh->GetMesh()->ExchangeData(check_mass_matrix, CELL, 0);
        BARRIER
        */

        // save stabilization
        /*
        INMOST::Tag check_stab_tag = mesh->GetMesh()->CreateTag("check_stab_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 2);
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                trianit->RealArray(check_stab_tag)[0] = (trianit->getFaces())[0]->RealArray(edge_stab_tags)[0];
                trianit->RealArray(check_stab_tag)[1] = (trianit->getFaces())[0]->RealArray(edge_stab_tags)[1];
            }
        }
        mesh->GetMesh()->ExchangeData(check_stab_tag, CELL, 0);
        BARRIER
        */

        // move velocity vector to geo basis
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if (edgeit->GetStatus() != Element::Ghost)
            {
                mesh->VecTransitionToGeoBasis
                (
                    std::vector<double>
                    {
                        edgeit->RealArray(vel_tag)[0],
                        edgeit->RealArray(vel_tag)[1]
                    },
                    edgeit->getFaces()[0]
                );
            }
        }
        mesh->GetMesh()->ExchangeData(vel_tag, FACE, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::PrintProfiling()
    {
        std::cout << "profiling" << std::endl;
    }
}