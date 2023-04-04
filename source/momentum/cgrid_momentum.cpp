#include "momentum.hpp"

using namespace INMOST;

namespace SIMUG
{
    CgridMomentumSolver::CgridMomentumSolver(SIMUG::IceMesh* mesh_,
                                             double time_step_,
                                             velocity_tag vel_tag_,
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

        // compute mass matrix
        mass_matrix_tag = mesh->GetMesh()->CreateTag("mass_matrix_tag", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 6);
        //mesh->GetMesh()->SetFileOption("Tag:mass_matrix_tag", "nosave");
        mom_timer.Launch();
        ComputeMassMatrix(mass_matrix_tag);
        mom_timer.Stop();
        duration = mom_timer.GetMaxTime();
        mom_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("Computation of mass matrix: OK (" + std::to_string(duration) + " ms)\n");
        }
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

    void CgridMomentumSolver::ComputeMassMatrix(INMOST::Tag mm_tag)
    {
        INMOST::Tag trian_area_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize();

        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            ElementArray<Face> adj_edges = trianit->getFaces();
            double area = trianit->Real(trian_area_tag); 

            for (int ed_num = 0; ed_num < 3; ++ed_num)
            {
                if (adj_edges[ed_num]->GetStatus() != Element::Ghost)
                {
                    adj_edges[ed_num]->Real(mass_matrix_tag) += area/3.0;
                }
            }
        }
         mesh->GetMesh()->ExchangeData(mass_matrix_tag, INMOST::FACE, 0);
    }

    Cgrid_mEVP_Solver::Cgrid_mEVP_Solver(SIMUG::IceMesh* mesh_,
                                         double time_step_,
                                         velocity_tag vel_tag_,
                                         scalar_tag conc_tag_,
                                         scalar_tag thick_tag_,
                                         SIMUG::dyn::pressParam mom_press_param_,
                                         SIMUG::dyn::bcType mom_bc_type_,
                                         const std::vector<double>& real_params_,
                                         const std::vector<int>& integer_params_) :
            CgridMomentumSolver(mesh_,
                                time_step_,
                                vel_tag_,
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

        // create mandatory tags
        edge_thick_tag = mesh->GetMesh()->CreateTag("edge_thick_tag", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 1);
        edge_conc_tag = mesh->GetMesh()->CreateTag("edge_conc_tag", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 1);

        ua_tags = mesh->GetDataSingle(mesh::gridElemType::Edge)->Get(mesh::meshVar::ua);
        uw_tags = mesh->GetDataSingle(mesh::gridElemType::Edge)->Get(mesh::meshVar::uw);

        old_vel_tag = mesh->GetMesh()->CreateTag("old_vel_tag", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 2);
        prev_vel_tag = mesh->GetMesh()->CreateTag("prev_vel_tag", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 2);
        new_vel_tag = mesh->GetMesh()->CreateTag("new_vel_tag", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 2);

        prev_sigma_tag = mesh->GetMesh()->CreateTag("prev_sigma_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        new_sigma_tag = mesh->GetMesh()->CreateTag("new_sigma_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        old_sigma_tag = mesh->GetMesh()->CreateTag("old_sigma_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);

        delta_tag = mesh->GetDataSingle(mesh::gridElemType::Trian)->Get(mesh::meshVar::del);
        P_tag = mesh->GetDataSingle(mesh::gridElemType::Trian)->Get(mesh::meshVar::P0);
        vareps_tag = mesh->GetDataSingle(mesh::gridElemType::Trian)->Get(mesh::meshVar::eps);
        xi_edge_tag = mesh->GetMesh()->CreateTag("xi_edge_tag", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 1);
        force_tags = mesh->GetMesh()->CreateTag("force_tags", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 2);
        edge_stab_tags = mesh->GetMesh()->CreateTag("edge_stab_tags", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 2);
        edge_stab_sum_tags = mesh->GetMesh()->CreateTag("edge_stab_sum_tags", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 2);
        disc_level_tags = mesh->GetMesh()->CreateTag("disc_level_tags", INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 2);
        level_tag = mesh->GetDataSingle(mesh::gridElemType::Trian)->Get(mesh::meshVar::hw);

        shear_deformation_tag = mesh->GetMesh()->CreateTag("shear_deformation_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);

        check_vel_tag = mesh->GetMesh()->CreateTag("check_vel_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 6);
        
        // mute tags
        mesh->GetMesh()->SetFileOption("Tag:edge_thick_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:edge_conc_tag", "nosave");

        mesh->GetMesh()->SetFileOption("Tag:old_vel_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:prev_vel_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:new_vel_tag", "nosave");

        mesh->GetMesh()->SetFileOption("Tag:prev_sigma_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:new_sigma_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:old_sigma_tag", "nosave");

        mesh->GetMesh()->SetFileOption("Tag:xi_edge_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:xi_node_tag", "nosave");
        //mesh->GetMesh()->SetFileOption("Tag:force_tags", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:edge_stab_tags", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:edge_stab_sum_tags", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:disc_level_tags", "nosave");

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("=====================================================================\n");
        }
        
        BARRIER
    }

    void Cgrid_mEVP_Solver::UpdateScalars()
    {
        // fix scalars
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                double a = trianit->Real(conc_tag);
                double h = trianit->Real(thick_tag);

                
                if ((a < SIMUG::IceConsts::amin) or (h < SIMUG::IceConsts::hmin))
                {
                    trianit->Real(conc_tag) = SIMUG::IceConsts::amin;
                    trianit->Real(thick_tag) = SIMUG::IceConsts::hmin;
                }

                if (a > 1.0)
                {
                    trianit->Real(conc_tag) = 1.0;
                }
            }
        }

        mesh->GetMesh()->ExchangeData(conc_tag, CELL, 0);
        mesh->GetMesh()->ExchangeData(thick_tag, CELL, 0);
        BARRIER

        // zeroing edge scalars
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if ((edgeit->GetStatus() != Element::Ghost) &&
                (edgeit->Integer(mesh->GetGridInfo(mesh::gridElemType::Edge)->is_bnd) == 0))
            {
                ElementArray<Cell> adj_trians = edgeit->getCells();

                edgeit->Real(edge_thick_tag) = 0.5*(adj_trians[0]->Real(thick_tag) + adj_trians[1]->Real(thick_tag));
                edgeit->Real(edge_conc_tag) = 0.5*(adj_trians[0]->Real(conc_tag) + adj_trians[1]->Real(conc_tag));
            }

            if ((edgeit->GetStatus() != Element::Ghost) &&
                (edgeit->Integer(mesh->GetGridInfo(mesh::gridElemType::Edge)->is_bnd) == 1))
            {
                ElementArray<Cell> adj_trians = edgeit->getCells();

                edgeit->Real(edge_thick_tag) = adj_trians[0]->Real(thick_tag);
                edgeit->Real(edge_conc_tag) = adj_trians[0]->Real(conc_tag);
            }
        }

        mesh->GetMesh()->ExchangeData(edge_thick_tag, FACE, 0);
        mesh->GetMesh()->ExchangeData(edge_conc_tag, FACE, 0);
        BARRIER
    }


    void Cgrid_mEVP_Solver::ComputeP()
    {
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                double h = trianit->Real(thick_tag);

                double a = trianit->Real(conc_tag);
                double delta = trianit->Real(delta_tag);

                double p_str = IceConsts::pstr;
                double C = IceConsts::C;
                double del_min = IceConsts::delmin;

                trianit->Real(P_tag) = p_str*h*std::exp(-C*(1.0 - a));
            }
        }

        mesh->GetMesh()->ExchangeData(P_tag, CELL, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::ComputeVarepsilonDelta()
    {
        double e = IceConsts::e;

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
                        std::vector<double>
                        {
                            edges[i]->RealArray(new_vel_tag)[0],
                            edges[i]->RealArray(new_vel_tag)[1]
                        }
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

                // assemble gradient vectors
                std::vector<double> grad_u(2);
                for (int i = 0; i < 3; ++i)
                {
                    grad_u = grad_u + ue_T[i]*
                    std::vector<double>
                    {
                        trianit->RealArray(grad_basis_func_tags)[i*2 + 0],
                        trianit->RealArray(grad_basis_func_tags)[i*2 + 1]
                    };
                }

                std::vector<double> grad_v(2);
                for (int i = 0; i < 3; ++i)
                {
                    grad_v = grad_v + ve_T[i]*
                    std::vector<double>
                    {
                        trianit->RealArray(grad_basis_func_tags)[i*2 + 0],
                        trianit->RealArray(grad_basis_func_tags)[i*2 + 1]
                    };
                }

                double du_dx = grad_u[0];
                double du_dy = grad_u[1];

                double dv_dx = grad_v[0];
                double dv_dy = grad_v[1];

                // compute strain rate tensor
                double eps11 = du_dx;
                double eps22 = dv_dy;
                double eps12 = 0.5*(du_dy + dv_dx);

                // compute delta
                double delta = std::sqrt( (eps11*eps11 + eps22*eps22)*(1.0 + 1.0/(e*e)) +
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

    void Cgrid_mEVP_Solver::AssembleForceVector()
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
            // get area of trian
            double trian_area = trianit->Real(trian_area_tag);

            ElementArray<Face> adj_edges = trianit->getFaces();

            for (int ed_num = 0; ed_num < 3; ++ed_num)
            {
                if (adj_edges[ed_num]->GetStatus() != Element::Ghost)
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
                    double sigma1 = trianit->RealArray(new_sigma_tag)[0];
                    double sigma2 = trianit->RealArray(new_sigma_tag)[1];
                    double sigma12 = trianit->RealArray(new_sigma_tag)[2];

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
                    adj_edges[ed_num].RealArray(force_tags)[0] += (adj_edges[ed_num]->GetStatus() != Element::Ghost) ? trian_area*(sig_edge[0][0]*grad_basis[0] + sig_edge[0][1]*grad_basis[1]) : 0.0;
                    adj_edges[ed_num].RealArray(force_tags)[1] += (adj_edges[ed_num]->GetStatus() != Element::Ghost) ? trian_area*(sig_edge[1][0]*grad_basis[0] + sig_edge[1][1]*grad_basis[1]) : 0.0;
                }
            }
        }

        mesh->GetMesh()->ExchangeData(force_tags, FACE, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::UpdateSigmaMevp()
    {
        double alpha = real_params[0];
        double e = IceConsts::e;
        double del_min = IceConsts::delmin;

        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                trianit->RealArray(new_sigma_tag)[0] = (1.0 - 1.0/alpha)*trianit->RealArray(prev_sigma_tag)[0] + 
                                                       (1.0/alpha)*(trianit->Real(P_tag)/(trianit->Real(delta_tag) + del_min))*
                                                       (trianit->RealArray(vareps_tag)[0] - trianit->Real(delta_tag));
                
                trianit->RealArray(new_sigma_tag)[1] = (1.0 - 1.0/alpha)*trianit->RealArray(prev_sigma_tag)[1] + 
                                                       (1.0/alpha)*(trianit->Real(P_tag)/((trianit->Real(delta_tag) + del_min)*e*e))*
                                                       trianit->RealArray(vareps_tag)[1];
                
                trianit->RealArray(new_sigma_tag)[2] = (1.0 - 1.0/alpha)*trianit->RealArray(prev_sigma_tag)[2] + 
                                                       (1.0/alpha)*(trianit->Real(P_tag)/((trianit->Real(delta_tag) + del_min)*e*e))*
                                                       trianit->RealArray(vareps_tag)[2];
            }
        }

        mesh->GetMesh()->ExchangeData(new_sigma_tag, CELL, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::ComputeEdgeStabilizationSum()
    {
        INMOST::Tag is_edge_bnd = mesh->GetGridInfo(mesh::gridElemType::Edge)->is_bnd;

        // zerroing stabilization sum
        for (auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
                edgeit->RealArray(edge_stab_sum_tags)[0] = 0.0;
                edgeit->RealArray(edge_stab_sum_tags)[1] = 0.0;
        }
        mesh->GetMesh()->ExchangeData(edge_stab_sum_tags, FACE, 0);
        BARRIER

        // compute new stabilization sum
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            ElementArray<Face> edges = trianit->getFaces();

            // move velocity components to trian basis
            std::vector<std::vector<double>> trian_velocities;
            for (int i = 0; i < 3; ++i)
            {
                trian_velocities.push_back
                (
                    mesh->VecTransition
                    (
                        std::vector<double>
                        {
                            edges[i]->RealArray(new_vel_tag)[0],
                            edges[i]->RealArray(new_vel_tag)[1]
                        },
                        edges[i],
                        trianit->getCells()[0]
                    )
                );
            }

            // compute the difference of velocities on every edge (in edge basis)
            std::vector<double> e0_u_diff = 
            mesh->VecTransition(trian_velocities[1], trianit->getCells()[0], edges[0]) - 
            mesh->VecTransition(trian_velocities[2], trianit->getCells()[0], edges[0]);

            std::vector<double> e1_u_diff = 
            mesh->VecTransition(trian_velocities[2], trianit->getCells()[0], edges[1]) - 
            mesh->VecTransition(trian_velocities[0], trianit->getCells()[0], edges[1]);

            std::vector<double> e2_u_diff = 
            mesh->VecTransition(trian_velocities[0], trianit->getCells()[0], edges[2]) - 
            mesh->VecTransition(trian_velocities[1], trianit->getCells()[0], edges[2]);

            edges[0]->RealArray(edge_stab_sum_tags)[0] += (edges[0]->GetStatus() != Element::Ghost) ? e0_u_diff[0] : 0.0;
            edges[0]->RealArray(edge_stab_sum_tags)[1] += (edges[0]->GetStatus() != Element::Ghost) ? e0_u_diff[1] : 0.0;

            edges[1]->RealArray(edge_stab_sum_tags)[0] += (edges[1]->GetStatus() != Element::Ghost) ? e1_u_diff[0] : 0.0;
            edges[1]->RealArray(edge_stab_sum_tags)[1] += (edges[1]->GetStatus() != Element::Ghost) ? e1_u_diff[1] : 0.0;

            edges[2]->RealArray(edge_stab_sum_tags)[0] += (edges[2]->GetStatus() != Element::Ghost) ? e2_u_diff[0] : 0.0;
            edges[2]->RealArray(edge_stab_sum_tags)[1] += (edges[2]->GetStatus() != Element::Ghost) ? e2_u_diff[1] : 0.0;
        }
        mesh->GetMesh()->ExchangeData(edge_stab_sum_tags, FACE, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::ComputeEdgeStabilization(double alpha)
    {
        double e = IceConsts::e;
        double del_min = IceConsts::delmin;
        INMOST::Tag trian_area_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize();
        INMOST::Tag is_trian_bnd = mesh->GetGridInfo(mesh::gridElemType::Trian)->is_bnd;

        // compute stabilization sum
        ComputeEdgeStabilizationSum();
        BARRIER
        
        // compute edge xi
        for (auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if ((edgeit->GetStatus() != Element::Ghost) &&
                (edgeit->Integer(mesh->GetGridInfo(mesh::gridElemType::Edge)->is_bnd) == 0))
            {
                ElementArray<Cell> adj_trians = edgeit->getCells();
                double P0 = adj_trians[0]->Real(P_tag);
                double P1 = adj_trians[1]->Real(P_tag);
                double delta0 = adj_trians[0]->Real(delta_tag);
                double delta1 = adj_trians[1]->Real(delta_tag);
                double S0 = adj_trians[0]->Real(trian_area_tag);
                double S1 = adj_trians[1]->Real(trian_area_tag);

                //double xi0 = 3.5*0.5*(P0 + P1)*((S0 + S1)/3.0)/(time_step);

                edgeit->Real(xi_edge_tag) = 2.5*0.5*(P0 + P1)*((S0 + S1)/3.0)/time_step;
            }

            if ((edgeit->GetStatus() != Element::Ghost) &&
                (edgeit->Integer(mesh->GetGridInfo(mesh::gridElemType::Edge)->is_bnd) == 1))
            {
                ElementArray<Cell> adj_trians = edgeit->getCells();
                double P = adj_trians[0]->Real(P_tag);
                double delta = adj_trians[0]->Real(delta_tag);
                double S = adj_trians[0]->Real(trian_area_tag);

                edgeit->Real(xi_edge_tag) = 2.5*P*((2.0*S)/3.0)/time_step;
            }
        }
        mesh->GetMesh()->ExchangeData(xi_edge_tag, FACE, 0);
        BARRIER

        // zerroing edge stabilization
        for (auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
                edgeit->RealArray(edge_stab_tags)[0] = 0.0;
                edgeit->RealArray(edge_stab_tags)[1] = 0.0;
        }
        mesh->GetMesh()->ExchangeData(edge_stab_tags, FACE, 0);
        BARRIER

        // compute new edge stabilization
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            ElementArray<Face> edges = trianit->getFaces();

            // move stabilization sum components to trian basis
            std::vector<std::vector<double>> trian_stab_sum;
            for (int i = 0; i < 3; ++i)
            {
                trian_stab_sum.push_back
                (
                    mesh->VecTransition
                    (
                        std::vector<double>
                        {
                            edges[i]->RealArray(edge_stab_sum_tags)[0],
                            edges[i]->RealArray(edge_stab_sum_tags)[1]
                        },
                        edges[i],
                        trianit->getCells()[0]
                    )
                );
            }

            // get edges xi
            std::vector<double> xi;
            for (int i = 0; i < 3; ++i)
            {
                xi.push_back(edges[i]->Real(xi_edge_tag));
            }

            // compute new stabilization
            std::vector<double> e0_stab = 
            - xi[1]*mesh->VecTransition(trian_stab_sum[1], trianit->getCells()[0], edges[0])
            + xi[2]*mesh->VecTransition(trian_stab_sum[2], trianit->getCells()[0], edges[0]);

            std::vector<double> e1_stab = 
            - xi[2]*mesh->VecTransition(trian_stab_sum[2], trianit->getCells()[0], edges[1])
            + xi[0]*mesh->VecTransition(trian_stab_sum[0], trianit->getCells()[0], edges[1]);

            std::vector<double> e2_stab = 
            - xi[0]*mesh->VecTransition(trian_stab_sum[0], trianit->getCells()[0], edges[2])
            + xi[1]*mesh->VecTransition(trian_stab_sum[1], trianit->getCells()[0], edges[2]);

            // store
            edges[0]->RealArray(edge_stab_tags)[0] += (edges[0]->GetStatus() != Element::Ghost) ? (-alpha/3.0)*e0_stab[0]: 0.0;
            edges[0]->RealArray(edge_stab_tags)[1] += (edges[0]->GetStatus() != Element::Ghost) ? (-alpha/3.0)*e0_stab[1]: 0.0;

            edges[1]->RealArray(edge_stab_tags)[0] += (edges[1]->GetStatus() != Element::Ghost) ? (-alpha/3.0)*e1_stab[0]: 0.0;
            edges[1]->RealArray(edge_stab_tags)[1] += (edges[1]->GetStatus() != Element::Ghost) ? (-alpha/3.0)*e1_stab[1]: 0.0;

            edges[2]->RealArray(edge_stab_tags)[0] += (edges[2]->GetStatus() != Element::Ghost) ? (-alpha/3.0)*e2_stab[0]: 0.0;
            edges[2]->RealArray(edge_stab_tags)[1] += (edges[2]->GetStatus() != Element::Ghost) ? (-alpha/3.0)*e2_stab[1]: 0.0;
        }
        mesh->GetMesh()->ExchangeData(edge_stab_tags, FACE, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::ComputeLevelVector()
    {
        double rho_ice = IceConsts::rhoi;
        double g = GenConsts::g;
        INMOST::Tag trian_area_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize();

        // zerroing current discrete level vector
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if ((edgeit->GetStatus() != Element::Ghost))
            {
                edgeit->RealArray(disc_level_tags)[0] = 0.0;
                edgeit->RealArray(disc_level_tags)[1] = 0.0;
            }
        }
        mesh->GetMesh()->ExchangeData(disc_level_tags, FACE, 0);

        // compute new discrete level vector
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            ElementArray<Face> adj_edges = trianit->getFaces();
            double level = trianit->Real(level_tag);
            double thickness = trianit->Real(thick_tag);
            double area = trianit->Real(trian_area_tag);

            for (int ed_num = 0; ed_num < 3; ++ed_num)
            {
                if ((adj_edges[ed_num]->GetStatus() != Element::Ghost))
                {
                    // get gradient in current edge basis
                    std::vector<double> current_grad = mesh->VecTransition
                    (
                        std::vector<double>
                        {
                            trianit->RealArray(grad_basis_func_tags)[2*ed_num + 0],
                            trianit->RealArray(grad_basis_func_tags)[2*ed_num + 1]
                        },
                        trianit->getCells()[0],
                        adj_edges[ed_num]
                    );

                    adj_edges[ed_num]->RealArray(disc_level_tags)[0] += rho_ice*thickness*g*level*current_grad[0]*area;
                    adj_edges[ed_num]->RealArray(disc_level_tags)[1] += rho_ice*thickness*g*level*current_grad[1]*area;
                }
            }
        }
        mesh->GetMesh()->ExchangeData(disc_level_tags, FACE, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::UpdateVelocityMevp()
    {
        double rho_i = IceConsts::rhoi;
        double rho_w = WaterConsts::rhow;
        double rho_a = AirConsts::rhoa;
        double Cw = IceConsts::Cw;
        double Ca = IceConsts::Ca;
        double dt = time_step;
        double beta = real_params[1];
        double f = 0.0;// GenConsts::f;

        INMOST::Tag edge_id_tag = mesh->GetGridInfo(mesh::gridElemType::Edge)->id;
        INMOST::Tag is_node_bnd_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->is_bnd;
        INMOST::Tag is_edge_bnd_tag = mesh->GetGridInfo(mesh::gridElemType::Edge)->is_bnd;
        INMOST::Tag is_trian_bnd_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->is_bnd;

        // assemble rhs, lhs and solve
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if ((edgeit->GetStatus() != Element::Ghost) and (edgeit->Integer(is_edge_bnd_tag) == 0))
            {
                std::vector<double> u_old = 
                {
                    edgeit->RealArray(old_vel_tag)[0],
                    edgeit->RealArray(old_vel_tag)[1]
                };

                std::vector<double> u_prev = 
                {
                    edgeit->RealArray(prev_vel_tag)[0],
                    edgeit->RealArray(prev_vel_tag)[1]
                };

                std::vector<double> u_a = 
                {
                    edgeit->RealArray(ua_tags)[0],
                    edgeit->RealArray(ua_tags)[1]
                };

                std::vector<double> u_w = 
                {
                    edgeit->RealArray(uw_tags)[0],
                    edgeit->RealArray(uw_tags)[1]
                };

                std::vector<double> Force = 
                {
                    edgeit->RealArray(force_tags)[0],
                    edgeit->RealArray(force_tags)[1]
                };

                std::vector<double> Level = 
                {
                    edgeit->RealArray(disc_level_tags)[0],
                    edgeit->RealArray(disc_level_tags)[1]
                };

                // Stabilization
                std::vector<double> Stab = 
                {
                    edgeit->RealArray(edge_stab_tags)[0],
                    edgeit->RealArray(edge_stab_tags)[1]
                };

                double h = edgeit->Real(edge_thick_tag);
                double a = edgeit->Real(edge_conc_tag);

                double mm = edgeit->Real(mass_matrix_tag);

                /*
                double lhs = (beta + 1.0)*rho_i*h/dt + rho_w*Cw*L2_norm_vec(u_prev - u_w);                
                */

                double lhs = (beta + 1.0)*rho_i*h/dt + a*rho_w*Cw*L2_norm_vec(u_prev - u_w);


                std::vector<double> rhs(2);

                /*
                rhs = (rho_i*h/dt)*(beta*u_prev + u_old) + (1.0/mm)*(Level - Force)
                      + rho_a*Ca*L2_norm_vec(u_a)*u_a 
                      + rho_w*Cw*L2_norm_vec(u_prev - u_w)*u_w
                      + rho_i*h*f*std::vector<double>{-u_prev[1], u_prev[0]}
                      + (1.0/mm)*Stab;
                */

                rhs = (rho_i*h/dt)*(beta*u_prev + u_old) + (1.0/mm)*(Level - Force)
                      + a*rho_a*Ca*L2_norm_vec(u_a)*u_a 
                      + a*rho_w*Cw*L2_norm_vec(u_prev - u_w)*u_w
                      + rho_i*h*f*std::vector<double>{-u_prev[1], u_prev[0]}
                      + (1.0/mm)*Stab;


                edgeit->RealArray(new_vel_tag)[0] = rhs[0]/lhs;
                edgeit->RealArray(new_vel_tag)[1] = rhs[1]/lhs;
            }
        }
        mesh->GetMesh()->ExchangeData(new_vel_tag, FACE, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::UpdateVelocity()
    {
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if (edgeit->GetStatus() != Element::Ghost)
            {
                edgeit->RealArray(prev_vel_tag)[0] = edgeit->RealArray(new_vel_tag)[0];
                edgeit->RealArray(prev_vel_tag)[1] = edgeit->RealArray(new_vel_tag)[1];
            }
        }
        mesh->GetMesh()->ExchangeData(prev_vel_tag, FACE, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::UpdateSigma()
    {
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                trianit->RealArray(prev_sigma_tag)[0] = trianit->RealArray(new_sigma_tag)[0];
                trianit->RealArray(prev_sigma_tag)[1] = trianit->RealArray(new_sigma_tag)[1];
                trianit->RealArray(prev_sigma_tag)[2] = trianit->RealArray(new_sigma_tag)[2];
            }
        }
        mesh->GetMesh()->ExchangeData(prev_sigma_tag, CELL, 0);
        BARRIER
    }


    void Cgrid_mEVP_Solver::MoveVectors(bool to_local)
    {
        if (to_local)
        {
            // move velocity vectors to node basis
            for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
            {
                if (edgeit->GetStatus() != Element::Ghost)
                {
                    // ice velocity
                    std::vector<double> edge_ice_vel = 
                    mesh->VecTransitionToElemBasis
                    (
                        std::vector<double>
                        {
                            edgeit->RealArray(vel_tag)[0],
                            edgeit->RealArray(vel_tag)[1]
                        },
                        edgeit->getFaces()[0]
                    );
                    edgeit->RealArray(vel_tag)[0] = edge_ice_vel[0]; edgeit->RealArray(vel_tag)[1] = edge_ice_vel[1]; 

                    // air velocity
                    std::vector<double> edge_air_vel = 
                    mesh->VecTransitionToElemBasis
                    (
                        std::vector<double>
                        {
                            edgeit->RealArray(ua_tags)[0],
                            edgeit->RealArray(ua_tags)[1]
                        },
                        edgeit->getFaces()[0]
                    );

                    edgeit->RealArray(ua_tags)[0] = edge_air_vel[0]; edgeit->RealArray(ua_tags)[1] = edge_air_vel[1]; 

                    // water velocity
                    std::vector<double> edge_water_vel = 
                    mesh->VecTransitionToElemBasis
                    (
                        std::vector<double>
                        {
                            edgeit->RealArray(uw_tags)[0],
                            edgeit->RealArray(uw_tags)[1]
                        },
                        edgeit->getFaces()[0]
                    );
                    edgeit->RealArray(uw_tags)[0] = edge_water_vel[0]; edgeit->RealArray(uw_tags)[1] = edge_water_vel[1]; 
                }
            }
        }
        else
        {
            for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
            {
                if (edgeit->GetStatus() != Element::Ghost)
                {
                    // ice velocity
                    std::vector<double> geo_ice_vel = 
                    mesh->VecTransitionToGeoBasis
                    (
                        std::vector<double>
                        {
                            edgeit->RealArray(vel_tag)[0],
                            edgeit->RealArray(vel_tag)[1]
                        },
                        edgeit->getFaces()[0]
                    );
                    edgeit->RealArray(vel_tag)[0] = geo_ice_vel[0]; edgeit->RealArray(vel_tag)[1] = geo_ice_vel[1]; 

                    // air velocity
                    std::vector<double> geo_air_vel = 
                    mesh->VecTransitionToGeoBasis
                    (
                        std::vector<double>
                        {
                            edgeit->RealArray(ua_tags)[0],
                            edgeit->RealArray(ua_tags)[1]
                        },
                        edgeit->getFaces()[0]
                    );
                    edgeit->RealArray(ua_tags)[0] = geo_air_vel[0]; edgeit->RealArray(ua_tags)[1] = geo_air_vel[1]; 

                    // water velocity
                    std::vector<double> geo_water_vel = 
                    mesh->VecTransitionToGeoBasis
                    (
                        std::vector<double>
                        {
                            edgeit->RealArray(uw_tags)[0],
                            edgeit->RealArray(uw_tags)[1]
                        },
                        edgeit->getFaces()[0]
                    );
                    edgeit->RealArray(uw_tags)[0] = geo_water_vel[0]; edgeit->RealArray(uw_tags)[1] = geo_water_vel[1]; 
                }
            }
        }

        mesh->GetMesh()->ExchangeData(vel_tag, FACE, 0);
        mesh->GetMesh()->ExchangeData(ua_tags, FACE, 0);
        mesh->GetMesh()->ExchangeData(uw_tags, FACE, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::Initialize()
    {
        // assign velocities from the previous step
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if (edgeit->GetStatus() != Element::Ghost)
            {
                edgeit->RealArray(old_vel_tag)[0] = edgeit->RealArray(vel_tag)[0];
                edgeit->RealArray(old_vel_tag)[1] = edgeit->RealArray(vel_tag)[1];

                edgeit->RealArray(prev_vel_tag)[0] = edgeit->RealArray(vel_tag)[0];
                edgeit->RealArray(prev_vel_tag)[1] = edgeit->RealArray(vel_tag)[1];

                edgeit->RealArray(new_vel_tag)[0] = edgeit->RealArray(vel_tag)[0];
                edgeit->RealArray(new_vel_tag)[1] = edgeit->RealArray(vel_tag)[1];
            }
        }
        mesh->GetMesh()->ExchangeData(old_vel_tag, FACE, 0);
        mesh->GetMesh()->ExchangeData(prev_vel_tag, FACE, 0);
        mesh->GetMesh()->ExchangeData(new_vel_tag, FACE, 0);

        // assign sigma from the previous step
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                trianit->RealArray(old_sigma_tag)[0] = trianit->RealArray(new_sigma_tag)[0];
                trianit->RealArray(old_sigma_tag)[1] = trianit->RealArray(new_sigma_tag)[1];
                trianit->RealArray(old_sigma_tag)[2] = trianit->RealArray(new_sigma_tag)[2];

                trianit->RealArray(prev_sigma_tag)[0] = trianit->RealArray(new_sigma_tag)[0];
                trianit->RealArray(prev_sigma_tag)[1] = trianit->RealArray(new_sigma_tag)[1];
                trianit->RealArray(prev_sigma_tag)[2] = trianit->RealArray(new_sigma_tag)[2];
            }
        }
        mesh->GetMesh()->ExchangeData(old_sigma_tag, CELL, 0);
        mesh->GetMesh()->ExchangeData(prev_sigma_tag, CELL, 0);

        BARRIER
    }

    void Cgrid_mEVP_Solver::Finalize()
    {
        // store velocity
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if (edgeit->GetStatus() != Element::Ghost)
            {
                edgeit->RealArray(vel_tag)[0] = edgeit->RealArray(new_vel_tag)[0];
                edgeit->RealArray(vel_tag)[1] = edgeit->RealArray(new_vel_tag)[1];
            }
        }
        mesh->GetMesh()->ExchangeData(vel_tag, FACE, 0);
        BARRIER

        // store sigma
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                trianit->RealArray(old_sigma_tag)[0] = trianit->RealArray(new_sigma_tag)[0];
                trianit->RealArray(old_sigma_tag)[1] = trianit->RealArray(new_sigma_tag)[1];
                trianit->RealArray(old_sigma_tag)[2] = trianit->RealArray(new_sigma_tag)[2];
            }
        }
        mesh->GetMesh()->ExchangeData(old_sigma_tag, CELL, 0);

        BARRIER
    }

    void Cgrid_mEVP_Solver::LogError(int presudostep)
    {
        // compute u diff norm
        double u_diff_norm = 0.0;
        double u_old_norm = 0.0;

        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if (edgeit->GetStatus() != Element::Ghost)
            {
                double current_vec_diff_sqr = L2_norm_vec
                (
                    std::vector<double>
                    {
                        edgeit->RealArray(new_vel_tag)[0] - edgeit->RealArray(prev_vel_tag)[0],
                        edgeit->RealArray(new_vel_tag)[1] - edgeit->RealArray(prev_vel_tag)[1]
                    }
                );
                current_vec_diff_sqr *= current_vec_diff_sqr;
                u_diff_norm += current_vec_diff_sqr;

                double current_old_vec_sqr = L2_norm_vec
                (
                    std::vector<double>
                    {
                        edgeit->RealArray(old_vel_tag)[0],
                        edgeit->RealArray(old_vel_tag)[1]
                    }
                );
                current_old_vec_sqr *= current_old_vec_sqr;
                u_old_norm += current_old_vec_sqr;
            }
        }
        BARRIER

        // compute sum on every process
        double u_diff_norm_all = u_diff_norm;
        double u_old_norm_all = u_old_norm;

#if defined(USE_MPI)
        MPI_Allreduce(&u_diff_norm_all, &u_diff_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        BARRIER

#if defined(USE_MPI)
        MPI_Allreduce(&u_old_norm_all, &u_old_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        BARRIER

        u_diff_norm_all = std::sqrt(u_diff_norm_all);
        u_old_norm_all = std::sqrt(u_old_norm_all);

        u_old_norm_all = (u_old_norm_all < REAL_MIN_ABS_VAL) ? 1e15 : u_old_norm_all;

        // compute sigma diff norm
        double sigma_diff_norm = 0.0;
        double sigma_old_norm = 0.0;

        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                double current_sigma_diff_sqr = L2_norm_vec
                (
                    std::vector<double>
                    {
                        trianit->RealArray(new_sigma_tag)[0] - trianit->RealArray(prev_sigma_tag)[0],
                        trianit->RealArray(new_sigma_tag)[1] - trianit->RealArray(prev_sigma_tag)[1],
                        trianit->RealArray(new_sigma_tag)[2] - trianit->RealArray(prev_sigma_tag)[2],
                    }
                );
                current_sigma_diff_sqr *= current_sigma_diff_sqr;
                sigma_diff_norm += current_sigma_diff_sqr;

                double current_sigma_old_sqr = L2_norm_vec
                (
                    std::vector<double>
                    {
                        trianit->RealArray(old_sigma_tag)[0],
                        trianit->RealArray(old_sigma_tag)[1],
                        trianit->RealArray(old_sigma_tag)[2],
                    }
                );
                current_sigma_old_sqr *= current_sigma_old_sqr;
                sigma_old_norm += current_sigma_old_sqr;
            }
        }
        BARRIER

        // compute sum on every process
        double sigma_diff_norm_all = sigma_diff_norm;
        double sigma_old_norm_all = sigma_old_norm;

#if defined(USE_MPI)
        MPI_Allreduce(&sigma_diff_norm_all, &sigma_diff_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        BARRIER

#if defined(USE_MPI)
        MPI_Allreduce(&sigma_old_norm_all, &sigma_old_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        BARRIER

        sigma_diff_norm_all = std::sqrt(sigma_diff_norm_all);
        sigma_old_norm_all = std::sqrt(sigma_old_norm_all);

        sigma_old_norm_all = (sigma_old_norm_all < REAL_MIN_ABS_VAL) ? 1e20 : sigma_old_norm_all;

        double max_vel = GetMaxVelocity();
        BARRIER

        // log
        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            std::cout << "Pseudoit " << presudostep << ", vel error = " << u_diff_norm_all/u_old_norm_all 
                                                    << ", sigma error = " << sigma_diff_norm_all/sigma_old_norm_all 
                                                    << std::endl;
            
            std::cout << "max vel = " <<  max_vel << std::endl;
        }

        BARRIER
    }

    void Cgrid_mEVP_Solver::ComputeShearDeformation()
    {
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                double e11 = 0.5*(trianit->RealArray(vareps_tag)[0] + trianit->RealArray(vareps_tag)[1]); 
                double e22 = 0.5*(trianit->RealArray(vareps_tag)[0] - trianit->RealArray(vareps_tag)[1]);
                double e12 = trianit->RealArray(vareps_tag)[2];

                trianit->Real(shear_deformation_tag) = std::sqrt((e11 - e22)*(e11 - e22) + 4.0*e12*e12);
            }
        }
        mesh->GetMesh()->ExchangeData(shear_deformation_tag, CELL, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::ComputeVelocity()
    {
        double alpha_stab = real_params[2];
        int err_frequency = integer_params[0]/5;

        SIMUG::Timer mom_timer;
        SIMUG::Logger mom_log(std::cout);
        double duration;

        if (mesh->GetMesh()->GetProcessorRank() == 0)
            mom_log.Log("\nmEVP solver info \n");

        // compute edge-based scalars
        UpdateScalars();

        // move ice, air and water velocity to edge basis
        MoveVectors(true);

        // initialize mEVP solver
        Initialize();

        // compute ice strength
        ComputeP();

        // compute discrete level vector
        ComputeLevelVector();

        // log initial error
        //LogError(0);

        for (int pseudostep = 1; pseudostep < (integer_params[0] + 1); ++pseudostep)
        {
            // update vareps and delta
            mom_timer.Launch();
            ComputeVarepsilonDelta();
            mom_timer.Stop();
            duration = mom_timer.GetMaxTime();
            mom_timer.Reset();
            strain_rate_computation_time += duration;

            // update sigma
            UpdateSigmaMevp();

            // assemble force vector
            mom_timer.Launch();
            AssembleForceVector();
            mom_timer.Stop();
            duration = mom_timer.GetMaxTime();
            mom_timer.Reset();
            force_assembling_time += duration;

            // assemble edge stabilization
            mom_timer.Launch();
            ComputeEdgeStabilization(alpha_stab);
            mom_timer.Stop();
            duration = mom_timer.GetMaxTime();
            mom_timer.Reset();
            stabilization_assembling_time += duration;
            
            // update velocity
            mom_timer.Launch();
            UpdateVelocityMevp();
            mom_timer.Stop();
            duration = mom_timer.GetMaxTime();
            mom_timer.Reset();
            velocity_computation_time += duration;

            // Log error
            if (pseudostep%err_frequency == 0)
            {
                LogError(pseudostep);
            }

            // update velocity value
            UpdateVelocity();

            // update sigma value
            UpdateSigma();
        
            BARRIER
        }

        // compute shear deformation
        ComputeShearDeformation();

        // finalize mEVP solver
        Finalize();

        // move ice, air and water velocity back to geo basis
        MoveVectors(false);


        // storre edge velocities in trians
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                trianit->RealArray(check_vel_tag)[0] = (trianit->getFaces())[0]->RealArray(vel_tag)[0];
                trianit->RealArray(check_vel_tag)[1] = (trianit->getFaces())[0]->RealArray(vel_tag)[1];
                trianit->RealArray(check_vel_tag)[2] = (trianit->getFaces())[1]->RealArray(vel_tag)[0];
                trianit->RealArray(check_vel_tag)[3] = (trianit->getFaces())[1]->RealArray(vel_tag)[1];
                trianit->RealArray(check_vel_tag)[4] = (trianit->getFaces())[2]->RealArray(vel_tag)[0];
                trianit->RealArray(check_vel_tag)[5] = (trianit->getFaces())[2]->RealArray(vel_tag)[1];
            }
        }
        mesh->GetMesh()->ExchangeData(check_vel_tag, CELL, 0);
        BARRIER

        // Print profiling info
        PrintProfiling();
        BARRIER
    }

    void Cgrid_mEVP_Solver::PrintProfiling()
    {
        SIMUG::Logger mom_log(std::cout);

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("## Momentum solver profiling ##\n");
            mom_log.Log("Strain rate computation time: " + std::to_string(strain_rate_computation_time) + " ms\n");
            mom_log.Log("Force assembling time: " + std::to_string(force_assembling_time) + " ms\n");
            mom_log.Log("Velocity computation time: " + std::to_string(velocity_computation_time) + " ms\n");
            mom_log.Log("Stabilization assembling time: " + std::to_string(stabilization_assembling_time) + " ms\n");
        }

        strain_rate_computation_time = 0.0;
        force_assembling_time = 0.0;
        velocity_computation_time = 0.0;
        stabilization_assembling_time = 0.0;
        BARRIER
    }

    double Cgrid_mEVP_Solver::GetMaxVelocity()
    {

        double max_local = -1.0;
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            if (edgeit->GetStatus() != Element::Ghost)
            {
                double u = edgeit->RealArray(new_vel_tag)[0];
                double v = edgeit->RealArray(new_vel_tag)[1];
                double cur = std::sqrt(u*u + v*v);
                if (cur > max_local)
                {
                    max_local = cur;
                }
            }
        }
        BARRIER
        double max_vel_for_all = max_local;
        MPI_Allreduce(&max_local, &max_vel_for_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        BARRIER
        return max_vel_for_all;
    }
}