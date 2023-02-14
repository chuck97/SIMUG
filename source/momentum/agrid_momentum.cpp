#include "momentum.hpp"

using namespace INMOST;

namespace SIMUG
{
    AgridMomentumSolver::AgridMomentumSolver(SIMUG::IceMesh* mesh_,
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
                           SIMUG::dyn::spaceScheme::CFE,
                           mom_press_param_,
                           mom_bc_type_),
                           real_params(real_params_),
                           integer_params(integer_params_)
    {
        // initialize timer and logger
        SIMUG::Logger mom_log(std::cout);
        SIMUG::Timer mom_timer;
        double duration;

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
            mom_log.Log("Computation of basis functions gradients: OK (" + std::to_string(duration) + " ms)\n");
        }
        BARRIER

        // compute lumped mass matrix components
        lumped_mass_matrix_tag = mesh->GetMesh()->CreateTag("lumped_mass_matrix_tag", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
        mesh->GetMesh()->SetFileOption("Tag:lumped_mass_matrix_tag", "nosave");
        mom_timer.Launch();
        ComputeMassMatrix(lumped_mass_matrix_tag);
        mom_timer.Stop();
        duration = mom_timer.GetMaxTime();
        mom_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("Computation of lumped mass matrix components: OK (" + std::to_string(duration) + " ms)\n");
        }
        BARRIER
    }

    void AgridMomentumSolver::ComputeGradientsBasisFuncs(INMOST::Tag grad_bas_tags)
    {
        std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();

        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

                // get coordinates of nodes in trian basis
                std::vector<std::vector<double>> node_coords_tr_basis =
                {
                    {trianit->RealArray(node_coords_in_trian_basis_tags[0])[0], trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]},
                    {trianit->RealArray(node_coords_in_trian_basis_tags[1])[0], trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]},
                    {trianit->RealArray(node_coords_in_trian_basis_tags[2])[0], trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]}
                };


                for (int node_num = 0; node_num < 3; ++node_num)
                {
                    // get nodal values of current basis function
                    std::vector<double> basis_func_nodal_values = {0.0, 0.0, 0.0};
                    basis_func_nodal_values[node_num] = 1.0;
                    

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

                    trianit->RealArray(grad_basis_func_tags)[node_num*2 + 0] = coeffs[0];
                    trianit->RealArray(grad_basis_func_tags)[node_num*2 + 1] = coeffs[1];
                }
            }
        }
        mesh->GetMesh()->ExchangeData(grad_basis_func_tags, INMOST::CELL, 0);
        BARRIER
    }

    void AgridMomentumSolver::ComputeMassMatrix(INMOST::Tag mm_tag)
    {
        INMOST::Tag trian_area_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize();

        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

            double area = trianit->Real(trian_area_tag);

            for (int node_num = 0; node_num < 3; ++node_num)
            {
                if (adj_nodes[node_num]->GetStatus() != Element::Ghost)
                {
                    adj_nodes[node_num]->Real(lumped_mass_matrix_tag) += (area/3.0);
                }
            }
        }
        mesh->GetMesh()->ExchangeData(lumped_mass_matrix_tag, INMOST::NODE, 0);
        BARRIER
    }

    Agrid_mEVP_Solver::Agrid_mEVP_Solver(SIMUG::IceMesh* mesh_,
                                         double time_step_,
                                         velocity_tag vel_tag_,
                                         scalar_tag conc_tag_,
                                         scalar_tag thick_tag_,
                                         SIMUG::dyn::pressParam mom_press_param_,
                                         SIMUG::dyn::bcType mom_bc_type_,
                                         const std::vector<double>& real_params_,
                                         const std::vector<int>& integer_params_) :
            AgridMomentumSolver(mesh_,
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
        ua_tags = mesh->GetDataSingle(mesh::gridElemType::Node)->Get(mesh::meshVar::ua);
        uw_tags = mesh->GetDataSingle(mesh::gridElemType::Node)->Get(mesh::meshVar::uw);

        old_vel_tag = mesh->GetMesh()->CreateTag("old_vel_tag", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 2);
        prev_vel_tag = mesh->GetMesh()->CreateTag("prev_vel_tag", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 2);
        new_vel_tag = mesh->GetMesh()->CreateTag("new_vel_tag", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 2);

        prev_sigma_tag = mesh->GetMesh()->CreateTag("prev_sigma_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        new_sigma_tag = mesh->GetMesh()->CreateTag("new_sigma_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        old_sigma_tag = mesh->GetMesh()->CreateTag("old_sigma_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);

        delta_tag = mesh->GetDataSingle(mesh::gridElemType::Trian)->Get(mesh::meshVar::del);
        P_tag = mesh->GetDataSingle(mesh::gridElemType::Trian)->Get(mesh::meshVar::P0);
        vareps_tag = mesh->GetDataSingle(mesh::gridElemType::Trian)->Get(mesh::meshVar::eps);
        force_tags = mesh->GetMesh()->CreateTag("force_tags", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 2);
        disc_level_tags = mesh->GetMesh()->CreateTag("disc_level_tags", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 2);
        level_tag = mesh->GetDataSingle(mesh::gridElemType::Node)->Get(mesh::meshVar::hw);

        shear_deformation_tag = mesh->GetMesh()->CreateTag("shear_deformation_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);

        // mute tags
        mesh->GetMesh()->SetFileOption("Tag:old_vel_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:prev_vel_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:new_vel_tag", "nosave");

        mesh->GetMesh()->SetFileOption("Tag:prev_sigma_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:new_sigma_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:old_sigma_tag", "nosave");

        mesh->GetMesh()->SetFileOption("Tag:force_tags", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:disc_level_tags", "nosave");

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("=====================================================================\n");
        }
        BARRIER
    }

    void Agrid_mEVP_Solver::UpdateScalars()
    {
        // fix scalars
        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
        {
            if (nodeit->GetStatus() != Element::Ghost)
            {
                double a = nodeit->Real(conc_tag);
                double h = nodeit->Real(thick_tag);

                if ((a < SIMUG::IceConsts::amin) or (h < SIMUG::IceConsts::hmin))
                {
                    nodeit->Real(conc_tag) = SIMUG::IceConsts::amin;
                    nodeit->Real(thick_tag) = SIMUG::IceConsts::hmin;
                }

                if (a > 1.0)
                {
                    nodeit->Real(conc_tag) = 1.0;
                }
            }
        }

        mesh->GetMesh()->ExchangeData(conc_tag, NODE, 0);
        mesh->GetMesh()->ExchangeData(thick_tag, NODE, 0);
        BARRIER
    }

    void Agrid_mEVP_Solver::ComputeP()
    {
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                ElementArray<Node> nodes = trianit->getNodes();

                double h = (1.0/3.0)*(nodes[0]->Real(thick_tag) + nodes[1]->Real(thick_tag) + nodes[2]->Real(thick_tag));

                double a = (1.0/3.0)*(nodes[0]->Real(conc_tag) + nodes[1]->Real(conc_tag) + nodes[2]->Real(conc_tag));;

                double p_str = IceConsts::pstr;
                double C = IceConsts::C;

                trianit->Real(P_tag) = p_str*h*std::exp(-C*(1.0 - a));
            }
        }

        mesh->GetMesh()->ExchangeData(P_tag, CELL, 0);
        BARRIER
    }

    void Agrid_mEVP_Solver::ComputeVarepsilonDelta()
    {
        double e = IceConsts::e;

        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                ElementArray<Node> nodes = trianit->getNodes();

                // det node velocity components in node basis
                std::vector<std::vector<double>> node_velocities;
                for (int i = 0; i < 3; ++i)
                {
                    node_velocities.push_back
                    (
                        std::vector<double>
                        {
                            nodes[i]->RealArray(new_vel_tag)[0],
                            nodes[i]->RealArray(new_vel_tag)[1]
                        }
                    );
                }

                // move node velocity components to triangular basis
                std::vector<std::vector<double>> node_velocities_trian_basis;

                for (int i = 0; i < 3; ++i)
                {
                    node_velocities_trian_basis.push_back
                    (
                        mesh->VecTransition(node_velocities[i], nodes[i], trianit->getCells()[0])
                    );
                }

                // store u and v components
                std::vector<double> ue_T = std::vector<double>
                {
                    node_velocities_trian_basis[0][0],
                    node_velocities_trian_basis[1][0],
                    node_velocities_trian_basis[2][0]
                };

                std::vector<double> ve_T = std::vector<double>
                {
                    node_velocities_trian_basis[0][1],
                    node_velocities_trian_basis[1][1],
                    node_velocities_trian_basis[2][1]
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
                double del_min = IceConsts::delmin;

                double delta = std::sqrt( (eps11*eps11 + eps22*eps22)*(1.0 + 1.0/(e*e)) +
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

     void Agrid_mEVP_Solver::AssembleForceVector()
    {
        INMOST::Tag trian_area_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize();

        // zeroing of force vector
        for (auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
        {
            nodeit->RealArray(force_tags)[0] = 0.0;
            nodeit->RealArray(force_tags)[1] = 0.0;
        }

        mesh->GetMesh()->ExchangeData(force_tags, NODE, 0);
        BARRIER

        // recalculate force vector
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            // get area of trian
            double trian_area = trianit->Real(trian_area_tag);

            ElementArray<Node> adj_nodes = trianit->getNodes();

            for (int node_num = 0; node_num < 3; ++node_num)
            {
                if (adj_nodes[node_num]->GetStatus() != Element::Ghost)
                {
                    // get gradient of basis function in node basis
                    std::vector<double> grad_basis = 
                    {
                        mesh->VecTransition
                        (
                            std::vector<double>
                            {
                                trianit->RealArray(grad_basis_func_tags)[node_num*2 + 0],
                                trianit->RealArray(grad_basis_func_tags)[node_num*2 + 1]
                            },
                            trianit->getCells()[0],
                            adj_nodes[node_num]
                        )
                    };

                    // move sigma components to node basis
                    double sigma1 = trianit->RealArray(new_sigma_tag)[0];
                    double sigma2 = trianit->RealArray(new_sigma_tag)[1];
                    double sigma12 = trianit->RealArray(new_sigma_tag)[2];

                    double sigma11 = 0.5*(sigma1 + sigma2);
                    double sigma22 = 0.5*(sigma1 - sigma2);

                    std::vector<std::vector<double>> sig_node = 
                    {
                        mesh->TensTransition
                        (
                            std::vector<std::vector<double>>
                            {
                                {sigma11, sigma12},
                                {sigma12, sigma22}
                            },
                            trianit->getCells()[0],
                            adj_nodes[node_num]
                        )
                    };

                    // add force from trian to node
                    adj_nodes[node_num].RealArray(force_tags)[0] += trian_area*(sig_node[0][0]*grad_basis[0] + sig_node[0][1]*grad_basis[1]);
                    adj_nodes[node_num].RealArray(force_tags)[1] += trian_area*(sig_node[1][0]*grad_basis[0] + sig_node[1][1]*grad_basis[1]);
                }
            }
        }
        mesh->GetMesh()->ExchangeData(force_tags, NODE, 0);
        BARRIER
    }

    void Agrid_mEVP_Solver::UpdateSigmaMevp()
    {
        double alpha = real_params[0];
        double e = IceConsts::e;

        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                trianit->RealArray(new_sigma_tag)[0] = (1.0 - 1.0/alpha)*trianit->RealArray(prev_sigma_tag)[0] + 
                                                       (1.0/alpha)*(trianit->Real(P_tag)/trianit->Real(delta_tag))*
                                                       (trianit->RealArray(vareps_tag)[0] - trianit->Real(delta_tag));
                
                trianit->RealArray(new_sigma_tag)[1] = (1.0 - 1.0/alpha)*trianit->RealArray(prev_sigma_tag)[1] + 
                                                       (1.0/alpha)*(trianit->Real(P_tag)/(trianit->Real(delta_tag)*e*e))*
                                                       trianit->RealArray(vareps_tag)[1];
                
                trianit->RealArray(new_sigma_tag)[2] = (1.0 - 1.0/alpha)*trianit->RealArray(prev_sigma_tag)[2] + 
                                                       (1.0/alpha)*(trianit->Real(P_tag)/(trianit->Real(delta_tag)*e*e))*
                                                       trianit->RealArray(vareps_tag)[2];
            }
        }

        mesh->GetMesh()->ExchangeData(new_sigma_tag, CELL, 0);
        BARRIER
    }

    void Agrid_mEVP_Solver::ComputeLevelVector()
    {
        double rho_ice = IceConsts::rhoi;
        double g = GenConsts::g;
        INMOST::Tag trian_area_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize();

        // zerroing current discrete level vector
        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
        {
            if ((nodeit->GetStatus() != Element::Ghost) &&
                (nodeit->Integer(mesh->GetGridInfo(mesh::gridElemType::Node)->is_bnd) == 0))
            {
                nodeit->RealArray(disc_level_tags)[0] = 0.0;
                nodeit->RealArray(disc_level_tags)[1] = 0.0;
            }
        }
        mesh->GetMesh()->ExchangeData(disc_level_tags, NODE, 0);

        // compute new discrete level vector
        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            ElementArray<Node> adj_nodes = trianit->getNodes();
            double area = trianit->Real(trian_area_tag);

            for (int node_num = 0; node_num < 3; ++node_num)
            {
                if ((adj_nodes[node_num]->Integer(mesh->GetGridInfo(mesh::gridElemType::Node)->is_bnd) == 0) &&
                    (adj_nodes[node_num]->GetStatus() != Element::Ghost))
                {
                    // get nodal scalars
                    double level = adj_nodes[node_num]->Real(level_tag);
                    double thickness = adj_nodes[node_num]->Real(thick_tag);

                    // get gradient of basis function in current node basis
                    std::vector<double> current_grad = mesh->VecTransition
                    (
                        std::vector<double>
                        {
                            trianit->RealArray(grad_basis_func_tags)[2*node_num + 0],
                            trianit->RealArray(grad_basis_func_tags)[2*node_num + 1]
                        },
                        trianit->getCells()[0],
                        adj_nodes[node_num]
                    );

                    adj_nodes[node_num]->RealArray(disc_level_tags)[0] += rho_ice*thickness*g*level*current_grad[0]*area/3.0;
                    adj_nodes[node_num]->RealArray(disc_level_tags)[1] += rho_ice*thickness*g*level*current_grad[1]*area/3.0;
                }
            }
        }
        mesh->GetMesh()->ExchangeData(disc_level_tags, NODE, 0);
        BARRIER
    }

    void Agrid_mEVP_Solver::UpdateVelocityMevp()
    {
        double rho_i = IceConsts::rhoi;
        double rho_w = WaterConsts::rhow;
        double rho_a = AirConsts::rhoa;
        double Cw = IceConsts::Cw;
        double Ca = IceConsts::Ca;
        double dt = time_step;
        double beta = real_params[1];
        double f = GenConsts::f;

        // assemble rhs, lhs and solve
        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
        {
            if ((nodeit->GetStatus() != Element::Ghost) &&
                (nodeit->Integer(mesh->GetGridInfo(mesh::gridElemType::Node)->is_bnd) == 0))
            {
                std::vector<double> u_old = 
                {
                    nodeit->RealArray(old_vel_tag)[0],
                    nodeit->RealArray(old_vel_tag)[1]
                };

                std::vector<double> u_prev = 
                {
                    nodeit->RealArray(prev_vel_tag)[0],
                    nodeit->RealArray(prev_vel_tag)[1]
                };

                std::vector<double> u_a = 
                {
                    nodeit->RealArray(ua_tags)[0],
                    nodeit->RealArray(ua_tags)[1]
                };

                std::vector<double> u_w = 
                {
                    nodeit->RealArray(uw_tags)[0],
                    nodeit->RealArray(uw_tags)[1]
                };

                std::vector<double> Force = 
                {
                    nodeit->RealArray(force_tags)[0],
                    nodeit->RealArray(force_tags)[1]
                };

                std::vector<double> Level = 
                {
                    nodeit->RealArray(disc_level_tags)[0],
                    nodeit->RealArray(disc_level_tags)[1]
                };

                double h = nodeit->Real(thick_tag);
                double a = nodeit->Real(conc_tag);

                double mm = nodeit->Real(lumped_mass_matrix_tag);

                std::vector<double> lhs = 
                {
                    (beta + 1.0)*rho_i*h/dt + a*rho_w*Cw*L2_norm_vec(u_prev - u_w),
                    rho_i*h*f
                };

                std::vector<double> rhs(2);

                rhs = (rho_i*h/dt)*(beta*u_prev + u_old) + (1.0/mm)*(Level - Force)
                      + a*rho_a*Ca*L2_norm_vec(u_a)*u_a 
                      + a*rho_w*Cw*L2_norm_vec(u_prev - u_w)*u_w;
                
                nodeit->RealArray(new_vel_tag)[0] = (lhs[0]*rhs[0] + lhs[1]*rhs[1])/(lhs[0]*lhs[0] + lhs[1]*lhs[1]);
                nodeit->RealArray(new_vel_tag)[1] = (lhs[0]*rhs[1] - lhs[1]*rhs[0])/(lhs[0]*lhs[0] + lhs[1]*lhs[1]);                
            }
        }
        mesh->GetMesh()->ExchangeData(new_vel_tag, NODE, 0);
        BARRIER
    }

    void Agrid_mEVP_Solver::UpdateVelocity()
    {
        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
        {
            if (nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->RealArray(prev_vel_tag)[0] = nodeit->RealArray(new_vel_tag)[0];
                nodeit->RealArray(prev_vel_tag)[1] = nodeit->RealArray(new_vel_tag)[1];
            }
        }
        mesh->GetMesh()->ExchangeData(prev_vel_tag, NODE, 0);
        BARRIER
    }

    void Agrid_mEVP_Solver::UpdateSigma()
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

    void Agrid_mEVP_Solver::MoveVectors(bool to_local)
    {
        if (to_local)
        {
            // move velocity vectors to node basis
            for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
            {
                if (nodeit->GetStatus() != Element::Ghost)
                {
                    // ice velocity
                    std::vector<double> node_ice_vel = 
                    mesh->VecTransitionToElemBasis
                    (
                        std::vector<double>
                        {
                            nodeit->RealArray(vel_tag)[0],
                            nodeit->RealArray(vel_tag)[1]
                        },
                        nodeit->getNodes()[0]
                    );
                    nodeit->RealArray(vel_tag)[0] = node_ice_vel[0]; nodeit->RealArray(vel_tag)[1] = node_ice_vel[1]; 

                    // air velocity
                    std::vector<double> node_air_vel = 
                    mesh->VecTransitionToElemBasis
                    (
                        std::vector<double>
                        {
                            nodeit->RealArray(ua_tags)[0],
                            nodeit->RealArray(ua_tags)[1]
                        },
                        nodeit->getNodes()[0]
                    );

                    nodeit->RealArray(ua_tags)[0] = node_air_vel[0]; nodeit->RealArray(ua_tags)[1] = node_air_vel[1]; 

                    // water velocity
                    std::vector<double> node_water_vel = 
                    mesh->VecTransitionToElemBasis
                    (
                        std::vector<double>
                        {
                            nodeit->RealArray(uw_tags)[0],
                            nodeit->RealArray(uw_tags)[1]
                        },
                        nodeit->getNodes()[0]
                    );
                    nodeit->RealArray(uw_tags)[0] = node_water_vel[0]; nodeit->RealArray(uw_tags)[1] = node_water_vel[1]; 
                }
            }
        }
        else
        {
            for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
            {
                if (nodeit->GetStatus() != Element::Ghost)
                {
                    // ice velocity
                    std::vector<double> geo_ice_vel = 
                    mesh->VecTransitionToGeoBasis
                    (
                        std::vector<double>
                        {
                            nodeit->RealArray(vel_tag)[0],
                            nodeit->RealArray(vel_tag)[1]
                        },
                        nodeit->getNodes()[0]
                    );
                    nodeit->RealArray(vel_tag)[0] = geo_ice_vel[0]; nodeit->RealArray(vel_tag)[1] = geo_ice_vel[1]; 

                    // air velocity
                    std::vector<double> geo_air_vel = 
                    mesh->VecTransitionToGeoBasis
                    (
                        std::vector<double>
                        {
                            nodeit->RealArray(ua_tags)[0],
                            nodeit->RealArray(ua_tags)[1]
                        },
                        nodeit->getNodes()[0]
                    );
                    nodeit->RealArray(ua_tags)[0] = geo_air_vel[0]; nodeit->RealArray(ua_tags)[1] = geo_air_vel[1]; 

                    // water velocity
                    std::vector<double> geo_water_vel = 
                    mesh->VecTransitionToGeoBasis
                    (
                        std::vector<double>
                        {
                            nodeit->RealArray(uw_tags)[0],
                            nodeit->RealArray(uw_tags)[1]
                        },
                        nodeit->getNodes()[0]
                    );
                    nodeit->RealArray(uw_tags)[0] = geo_water_vel[0]; nodeit->RealArray(uw_tags)[1] = geo_water_vel[1]; 
                }
            }
        }

        mesh->GetMesh()->ExchangeData(vel_tag, NODE, 0);
        mesh->GetMesh()->ExchangeData(ua_tags, NODE, 0);
        mesh->GetMesh()->ExchangeData(uw_tags, NODE, 0);
        BARRIER
    }

    void Agrid_mEVP_Solver::Initialize()
    {
        // assign velocities from the previous step
        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
        {
            if (nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->RealArray(old_vel_tag)[0] = nodeit->RealArray(vel_tag)[0];
                nodeit->RealArray(old_vel_tag)[1] = nodeit->RealArray(vel_tag)[1];

                nodeit->RealArray(prev_vel_tag)[0] = nodeit->RealArray(vel_tag)[0];
                nodeit->RealArray(prev_vel_tag)[1] = nodeit->RealArray(vel_tag)[1];

                nodeit->RealArray(new_vel_tag)[0] = nodeit->RealArray(vel_tag)[0];
                nodeit->RealArray(new_vel_tag)[1] = nodeit->RealArray(vel_tag)[1];
            }
        }
        mesh->GetMesh()->ExchangeData(old_vel_tag, NODE, 0);
        mesh->GetMesh()->ExchangeData(prev_vel_tag, NODE, 0);
        mesh->GetMesh()->ExchangeData(new_vel_tag, NODE, 0);

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

    void Agrid_mEVP_Solver::Finalize()
    {
        // store velocity
        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
        {
            if (nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->RealArray(vel_tag)[0] = nodeit->RealArray(new_vel_tag)[0];
                nodeit->RealArray(vel_tag)[1] = nodeit->RealArray(new_vel_tag)[1];
            }
        }
        mesh->GetMesh()->ExchangeData(vel_tag, NODE, 0);
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

    void Agrid_mEVP_Solver::LogError(int presudostep)
    {
        // compute u diff norm
        double u_diff_norm = 0.0;
        double u_old_norm = 0.0;

        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
        {
            if (nodeit->GetStatus() != Element::Ghost)
            {
                double current_vec_diff_sqr = L2_norm_vec
                (
                    std::vector<double>
                    {
                        nodeit->RealArray(new_vel_tag)[0] - nodeit->RealArray(prev_vel_tag)[0],
                        nodeit->RealArray(new_vel_tag)[1] - nodeit->RealArray(prev_vel_tag)[1]
                    }
                );
                current_vec_diff_sqr *= current_vec_diff_sqr;
                u_diff_norm += current_vec_diff_sqr;

                double current_old_vec_sqr = L2_norm_vec
                (
                    std::vector<double>
                    {
                        nodeit->RealArray(old_vel_tag)[0],
                        nodeit->RealArray(old_vel_tag)[1]
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

        // log
        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            std::cout << "Pseudoit " << presudostep << ", vel error = " << u_diff_norm_all/u_old_norm_all 
                                                    << ", sigma error = " << sigma_diff_norm_all/sigma_old_norm_all 
                                                    << std::endl;
        }

        BARRIER
    }

    void Agrid_mEVP_Solver::ComputeShearDeformation()
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

    void Agrid_mEVP_Solver::ComputeVelocity()
    {
        int err_frequency = integer_params[0]/5;

        SIMUG::Timer mom_timer;
        SIMUG::Logger mom_log(std::cout);
        double duration;

        if (mesh->GetMesh()->GetProcessorRank() == 0)
            mom_log.Log("\nmEVP solver info \n");

        // compute node-based scalars
        UpdateScalars();

        // move ice, air and water velocity to node basis
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

        // log profiling
        PrintProfiling();
        BARRIER
    }

    void Agrid_mEVP_Solver::PrintProfiling()
    {
        SIMUG::Logger mom_log(std::cout);

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            mom_log.Log("## Momentum solver profiling ##\n");
            mom_log.Log("Strain rate computation time: " + std::to_string(strain_rate_computation_time) + " ms\n");
            mom_log.Log("Force assembling time: " + std::to_string(force_assembling_time) + " ms\n");
            mom_log.Log("Velocity computation time: " + std::to_string(velocity_computation_time) + " ms\n");
        }

        strain_rate_computation_time = 0.0;
        force_assembling_time = 0.0;
        velocity_computation_time = 0.0;
        BARRIER
    }
}