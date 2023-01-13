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
        mesh->GetMesh()->SetFileOption("Tag:mass_matrix_entry_tag", "nosave");
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

        // compute opposite node for every edge of triangle
        opposite_edge_for_node_tags = mesh->GetMesh()->CreateTag("opposite_edge_for_node_tags", INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 3);
        //mesh->GetMesh()->SetFileOption("Tag:opposite_node_for_edge_tags", "nosave");
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
                    int opposite_edge_num_for_node = 0;

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

                    // store opposite node num
                    trianit->IntegerArray(op_edge_tags)[node_num] = opposite_edge_num_for_node;
                }
            }
        }
        mesh->GetMesh()->ExchangeData(op_edge_tags, INMOST::CELL, 0);
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

        // create tag for P, strain rate and delta
        P_tag = mesh->GetMesh()->CreateTag("P_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
        vareps_tag = mesh->GetMesh()->CreateTag("vareps_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        delta_tag = mesh->GetMesh()->CreateTag("delta_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
        
        // mute tags
        //mesh->GetMesh()->SetFileOption("Tag:P_tag", "nosave");
        //mesh->GetMesh()->SetFileOption("Tag:vareps_tag", "nosave");
        //mesh->GetMesh()->SetFileOption("Tag:delta_tag", "nosave");

        ComputeVarepsilonDelta(vel_tag);

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

        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                ElementArray<Face> edges = trianit->getFaces();

                std::vector<double> ue0_e0, ue1_e1, ue2_e2;

                // compute components on edge velocities in edge basis
                if (trianit->Integer(mesh->GetGridInfo(mesh::gridElemType::Trian)->GetIsXedgeBasisIsNormal()[0]) == 0)
                {
                    ue0_e0 = std::vector<double> {edges[0]->RealArray(vel_tag)[1], -edges[0]->RealArray(vel_tag)[0]};
                }
                else
                {
                    ue0_e0 = std::vector<double> {-edges[0]->RealArray(vel_tag)[1], edges[0]->RealArray(vel_tag)[0]};
                }

                if (trianit->Integer(mesh->GetGridInfo(mesh::gridElemType::Trian)->GetIsXedgeBasisIsNormal()[1]) == 0)
                {
                    ue1_e1 = std::vector<double> {edges[1]->RealArray(vel_tag)[1], -edges[1]->RealArray(vel_tag)[0]};
                }
                else
                {
                    ue1_e1 = std::vector<double> {-edges[1]->RealArray(vel_tag)[1], edges[1]->RealArray(vel_tag)[0]};
                }

                if (trianit->Integer(mesh->GetGridInfo(mesh::gridElemType::Trian)->GetIsXedgeBasisIsNormal()[2]) == 0)
                {
                    ue2_e2 = std::vector<double> {edges[2]->RealArray(vel_tag)[1], -edges[2]->RealArray(vel_tag)[0]};
                }
                else
                {
                    ue2_e2 = std::vector<double> {-edges[2]->RealArray(vel_tag)[1], edges[2]->RealArray(vel_tag)[0]};
                }
                

                // move edge velocity components to triangular basis
                std::vector<double> ue0_T = mesh->VecTransition(ue0_e0, edges[0], trianit->getCells()[0]);
                std::vector<double> ue1_T = mesh->VecTransition(ue1_e1, edges[1], trianit->getCells()[0]);
                std::vector<double> ue2_T = mesh->VecTransition(ue2_e2, edges[2], trianit->getCells()[0]);

                std::vector<double> ue_T = {ue0_T[0], ue1_T[0], ue2_T[0]};
                std::vector<double> ve_T = {ue0_T[1], ue1_T[1], ue2_T[1]};

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
        BARRIER
        mesh->GetMesh()->ExchangeData(vareps_tag, CELL, 0);
        BARRIER
        mesh->GetMesh()->ExchangeData(delta_tag, CELL, 0);
        BARRIER
    }

    void Cgrid_mEVP_Solver::ComputeVelocity()
    {
        std::cout << "calculations" << std::endl;
    }

    void Cgrid_mEVP_Solver::PrintProfiling()
    {
        std::cout << "profiling" << std::endl;
    }
}