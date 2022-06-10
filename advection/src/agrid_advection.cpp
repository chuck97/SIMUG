#include "advection.hpp"

using namespace INMOST;
using namespace SIMUG;

AgridAdvectionSolver::AgridAdvectionSolver(SIMUG::IceMesh* mesh_,
                                           double time_step_,
                                           velocity_tag vel_tag_,
                                           INMOST::Solver* slae_solver_,
                                           adv::timeScheme adv_time_scheme_,
                                           adv::spaceScheme adv_space_scheme_,
                                           adv::advFilter adv_filter_,
                                           const std::vector<double>& params_):
            AdvectionSolver(mesh_, vel_tag_, time_step_),
            slae_solver(slae_solver_),
            adv_time_scheme(adv_time_scheme_),
            adv_space_scheme(adv_space_scheme_),
            adv_filter(adv_filter_),
            params(params_)
{
    // initialize timer and logger
    SIMUG::Logger adv_log(std::cout);
    SIMUG::Timer adv_timer;
    double duration;

    // log constructor
    if (mesh->GetMesh()->GetProcessorRank()==0)
    {
        adv_log.Log("================== Advection solver initialization ==================\n");
        adv_log.Log("Advection time scheme: " + adv::advTimeSchemeName.at(adv_time_scheme) + "\n");
        adv_log.Log("Advection space scheme: " + adv::advSpaceSchemeName.at(adv_space_scheme) + "\n");
        adv_log.Log("Advection filter: " + adv::advFilterName.at(adv_filter) + "\n");
    }
    BARRIER

    // compute maximal Courant number
    adv_timer.Launch();
    double max_courant = GetMaxCourant();
    adv_timer.Stop();
    duration = adv_timer.GetMaxTime();
    adv_timer.Reset();

    if (mesh->GetMesh()->GetProcessorRank()==0)
    {
        adv_log.Log("Maximal Courant number: " + std::to_string(max_courant) + " (" + std::to_string(duration) + " ms)\n");
        if (max_courant > 1.0)
            adv_log.Log("Courant number is greater than 1.0 - advection scheme may be unstable!\n");
    }
    
    // assemble LHS and setup linear solver
    adv_timer.Launch();
    AssembleLHS();
    adv_timer.Stop();
    duration = adv_timer.GetMaxTime();
    adv_timer.Reset();

    if (mesh->GetMesh()->GetProcessorRank()==0)
        adv_log.Log("LHS for advection assembled successfully! (" + std::to_string(duration) + " ms)\n");

    if (mesh->GetMesh()->GetProcessorRank()==0)
        adv_log.Log("=====================================================================\n");
    BARRIER
}

double AgridAdvectionSolver::GetMaxCourant()
{
    // get cart coord tag for node
    INMOST::Tag node_cart_coords_tag =  mesh->GetGridInfo(mesh::gridElemType::Node)->coords[coord::coordType::cart];

    double max_courant = 0.0;

    for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
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

            double a0 = L2_norm_vec(first_node_coords - zero_node_coords);
            double a1 = L2_norm_vec(second_node_coords - zero_node_coords);
            double a2 = L2_norm_vec(second_node_coords - first_node_coords);

            double amin = std::min(std::min(a0, a1), a2);

            double u0 = L2_norm_vec(std::vector<double>{adj_nodes[0].RealArray(vel_tag)[0], adj_nodes[0].RealArray(vel_tag)[1]});
            double u1 = L2_norm_vec(std::vector<double>{adj_nodes[1].RealArray(vel_tag)[0], adj_nodes[1].RealArray(vel_tag)[1]});
            double u2 = L2_norm_vec(std::vector<double>{adj_nodes[2].RealArray(vel_tag)[0], adj_nodes[2].RealArray(vel_tag)[1]});

            double umax = std::max(std::max(u0, u1), u2);
            double courant = umax*time_step/amin;

            if (courant > max_courant)
                max_courant = courant;
        }
    }
    BARRIER
        
    // get max courant for all processors
    double max_courant_for_all;
    MPI_Allreduce(&max_courant, &max_courant_for_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return max_courant_for_all;
}

double AgridAdvectionSolver::GetMaxUdivDx()
{
    // get cart coord tag for node
    INMOST::Tag node_cart_coords_tag =  mesh->GetGridInfo(mesh::gridElemType::Node)->coords[coord::coordType::cart];

    double max_u_div_dx = 0.0;

    for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
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

            double a0 = L2_norm_vec(first_node_coords - zero_node_coords);
            double a1 = L2_norm_vec(second_node_coords - zero_node_coords);
            double a2 = L2_norm_vec(second_node_coords - first_node_coords);

            double amin = std::min(std::min(a0, a1), a2);

            double u0 = L2_norm_vec(std::vector<double>{adj_nodes[0].RealArray(vel_tag)[0], adj_nodes[0].RealArray(vel_tag)[1]});
            double u1 = L2_norm_vec(std::vector<double>{adj_nodes[1].RealArray(vel_tag)[0], adj_nodes[1].RealArray(vel_tag)[1]});
            double u2 = L2_norm_vec(std::vector<double>{adj_nodes[2].RealArray(vel_tag)[0], adj_nodes[2].RealArray(vel_tag)[1]});

            double umax = std::max(std::max(u0, u1), u2);
            double u_div_dx = umax/amin;

            if (u_div_dx > max_u_div_dx)
                max_u_div_dx = u_div_dx;
        }
    }
    BARRIER
        
    // get max u_div_dx for all processors if MPI is used
    double max_u_div_dx_for_all = max_u_div_dx;

#ifdef USE_MPI
    MPI_Allreduce(&max_u_div_dx, &max_u_div_dx_for_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    return max_u_div_dx_for_all;
}

void AgridAdvectionSolver::AssembleLHS()
{
    if (adv_space_scheme == adv::spaceScheme::CFE)
    {
        // get min and max node id
        int idmin = mesh->GetMeshInfo().id_interval_nodes.id_min;
        int idmax = mesh->GetMeshInfo().id_interval_nodes.id_max;

        // setup LHS partition
        LHS.SetInterval(idmin, idmax);

        // setup LHS_low for Zalesak filter
        if (adv_filter == adv::advFilter::Zalesak)
            LHS_low.SetInterval(idmin, idmax);

        // global node id tag
        INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

        // cart coords of nodes tag
        std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();

        // create tag for lumped mass matrix entry
        mesh->GetDataSingle(mesh::gridElemType::Node)->Create("M_L_entry", 1, INMOST::DATA_REAL);
        INMOST::Tag M_L_entry = mesh->GetDataSingle(mesh::gridElemType::Node)->Get("M_L_entry");

        for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit) 
        {
            // get nodes of current trian
            ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

            // get node coords of trian (in trian basis)
            std::vector<std::vector<double>> local_node_coords = 
            {
                {trianit->RealArray(node_coords_in_trian_basis_tags[0])[0], trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]},
                {trianit->RealArray(node_coords_in_trian_basis_tags[1])[0], trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]},
                {trianit->RealArray(node_coords_in_trian_basis_tags[2])[0], trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]}
            };

            // find global ids of local nodes
            std::vector<int> adj_node_global_ids = 
            { 
                adj_nodes[0].Integer(node_id_tag),
                adj_nodes[1].Integer(node_id_tag),
                adj_nodes[2].Integer(node_id_tag)
            }; 

            // assemble local 3x3 mass matrix
            std::vector<std::vector<double>> localLHS = LocalMassMatrixAssembling(local_node_coords);

            // assign local M_C_minus_M_L with M_C
            if (adv_filter == adv::advFilter::Zalesak)
            {
                if (trianit->GetStatus()!= Element::Ghost)
                    M_C_minus_M_L.push_back(localLHS);
            }

            // assemble global LHS mass matrix
            for (int i = 0; i < 3; ++i)
            {
                if(adj_nodes[i].GetStatus() != Element::Ghost)
                {
                    for (int j = 0; j < 3; ++j)
                        LHS[adj_node_global_ids[i]][adj_node_global_ids[j]] += localLHS[i][j];
                }
            }

            // assemble global lumped LHS mass matrix for scheme with Zalesak filter
            if (adv_filter == adv::advFilter::Zalesak)
            {
                for (int i = 0; i < 3; ++i)
                {
                    if(adj_nodes[i].GetStatus() != Element::Ghost)
                    {
                        for (int j = 0; j < 3; ++j)
                        {
                            LHS_low[adj_node_global_ids[i]][adj_node_global_ids[i]] += localLHS[i][j];
                            adj_nodes[i]->Real(M_L_entry) += localLHS[i][j];
                        }
                    }
                }

                // substract local M_L from local M_C_minus_M_L
                if (trianit->GetStatus() != Element::Ghost)
                {
                    for (int i = 0; i < 3; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                            M_C_minus_M_L.back()[i][i] -= localLHS[i][j];
                    }
                }
            } 
        }   
        BARRIER

        // exchange global lumped mass matrix entries
        if (adv_filter == adv::advFilter::Zalesak)
        {
            mesh->GetDataSingle(mesh::gridElemType::Node)->Exchange("M_L_entry");
            BARRIER
        }

        // multiply MC_minus_ML by (-ML^-1) for scheme with Zalesak filter
        if (adv_filter == adv::advFilter::Zalesak)
        {
            int t_num = 0;
            for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit) 
            {
                if (trianit->GetStatus() != Element::Ghost)
                {
                    INMOST::ElementArray<Node> adj_nodes = trianit->getNodes(); 
                    for (int i = 0; i < 3; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                            M_C_minus_M_L[t_num][i][j] /= -adj_nodes[i]->Real(M_L_entry);
                    }
                    ++t_num;
                }
            }
        }
        BARRIER

        // Delete temporal tag for lumped mass matrix
        mesh->GetDataSingle(mesh::gridElemType::Node)->Delete("M_L_entry");
    }
    else
    {
        SIMUG_ERR("only CFE space discretization availible for Agrid now!");
    }

    BARRIER

    // Set computed LHS matrix for linear solver
    slae_solver->SetMatrix(LHS);

    BARRIER
}

void AgridAdvectionSolver::AssembleSingleStepRHS(velocity_tag vel_tag, scalar_tag scal_tag)
{

    // node id tag
    INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

    // tag for node coords in local trian basis
    std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();

    for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit) 
    {
        // get nodes of current trian
        ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

        // get node coords of trian (in trian basis)
        std::vector<std::vector<double>> local_node_coords = 
        {
            {trianit->RealArray(node_coords_in_trian_basis_tags[0])[0], trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]},
            {trianit->RealArray(node_coords_in_trian_basis_tags[1])[0], trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]},
            {trianit->RealArray(node_coords_in_trian_basis_tags[2])[0], trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]}
        };

        // get velocity values of trian
        std::vector<std::vector<double>> local_vel_geo = 
        {
            {adj_nodes[0].RealArray(vel_tag)[0], adj_nodes[0].RealArray(vel_tag)[1]},
            {adj_nodes[1].RealArray(vel_tag)[0], adj_nodes[1].RealArray(vel_tag)[1]},
            {adj_nodes[2].RealArray(vel_tag)[0], adj_nodes[2].RealArray(vel_tag)[1]}
        };

        // move velocity components to nodal basis
        std::vector<std::vector<double>> local_vel_nodal = 
        {
            mesh->VecTransitionToElemBasis({local_vel_geo[0][0], local_vel_geo[0][1]}, adj_nodes[0]),
            mesh->VecTransitionToElemBasis({local_vel_geo[1][0], local_vel_geo[1][1]}, adj_nodes[1]),
            mesh->VecTransitionToElemBasis({local_vel_geo[2][0], local_vel_geo[2][1]}, adj_nodes[2])
        };

        // move velocity components to trian basis
        std::vector<std::vector<double>> local_vel_trian = 
        {
            mesh->VecTransition({local_vel_nodal[0][0], local_vel_nodal[0][1]}, adj_nodes[0], trianit->getCells()[0]),
            mesh->VecTransition({local_vel_nodal[1][0], local_vel_nodal[1][1]}, adj_nodes[1], trianit->getCells()[0]),
            mesh->VecTransition({local_vel_nodal[2][0], local_vel_nodal[2][1]}, adj_nodes[2], trianit->getCells()[0])
        };

        // get scalar values of trian
        std::vector<double> local_scal_trian = 
        {
            adj_nodes[0].Real(scal_tag),
            adj_nodes[1].Real(scal_tag),
            adj_nodes[2].Real(scal_tag)
        };

        // assemble local RHS vector of size 3
        std::vector<double> localRHS;

        if (adv_time_scheme == adv::timeScheme::TG2)
        {
            localRHS = LocalTG2RhsAssembling(local_node_coords,
                                             local_vel_trian,
                                             local_scal_trian,
                                             time_step);
        }
        else
        {
            SIMUG_ERR("Only avalible single step solvers: TG2");
        }

        // assemble global RHS vector
        for (int i = 0; i < 3; ++i)
        {
            if(adj_nodes[i].GetStatus() != Element::Ghost)
            {
                RHS[adj_nodes[i].Integer(node_id_tag)] += localRHS[i];
            }
        }
    }   
    BARRIER
}   

void AgridAdvectionSolver::AssembleDoubleStepRHS(velocity_tag vel_tag, scalar_tag scal_tag, scalar_tag scal_half_tag,  StepNumber step_num)
{
    // node id tag
    INMOST::Tag node_id = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

    // tag for node coords in local trian basis
    std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();

    for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit) 
    {
         // get nodes of current trian
        ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

        // get node coords of trian (in trian basis)
        std::vector<std::vector<double>> local_node_coords = 
        {
            {trianit->RealArray(node_coords_in_trian_basis_tags[0])[0], trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]},
            {trianit->RealArray(node_coords_in_trian_basis_tags[1])[0], trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]},
            {trianit->RealArray(node_coords_in_trian_basis_tags[2])[0], trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]}
        };


        // get velocity values of trian
        std::vector<std::vector<double>> local_vel_geo = 
        {
            {adj_nodes[0].RealArray(vel_tag)[0], adj_nodes[0].RealArray(vel_tag)[1]},
            {adj_nodes[1].RealArray(vel_tag)[0], adj_nodes[1].RealArray(vel_tag)[1]},
            {adj_nodes[2].RealArray(vel_tag)[0], adj_nodes[2].RealArray(vel_tag)[1]}
        };

        // move velocity components to nodal basis
        std::vector<std::vector<double>> local_vel_nodal = 
        {
            mesh->VecTransitionToElemBasis({local_vel_geo[0][0], local_vel_geo[0][1]}, adj_nodes[0]),
            mesh->VecTransitionToElemBasis({local_vel_geo[1][0], local_vel_geo[1][1]}, adj_nodes[1]),
            mesh->VecTransitionToElemBasis({local_vel_geo[2][0], local_vel_geo[2][1]}, adj_nodes[2])
        };

        // move velocity components to trian basis
        std::vector<std::vector<double>> local_vel_trian = 
        {
            mesh->VecTransition({local_vel_nodal[0][0], local_vel_nodal[0][1]}, adj_nodes[0], trianit->getCells()[0]),
            mesh->VecTransition({local_vel_nodal[1][0], local_vel_nodal[1][1]}, adj_nodes[1], trianit->getCells()[0]),
            mesh->VecTransition({local_vel_nodal[2][0], local_vel_nodal[2][1]}, adj_nodes[2], trianit->getCells()[0])
        };

        // get scalar values of trian
        std::vector<double> local_scal_trian = 
        {
            adj_nodes[0].Real(scal_tag),
            adj_nodes[1].Real(scal_tag),
            adj_nodes[2].Real(scal_tag)
        };

        // get scalar value of trian from previous step
        std::vector<double> local_scal_half_trian;
        if (step_num == StepNumber::second)
        {
            local_scal_half_trian = 
            {
                adj_nodes[0].Real(scal_half_tag),
                adj_nodes[1].Real(scal_half_tag),
                adj_nodes[2].Real(scal_half_tag)
            };
        }
        
        // assemble local load vector of size 3
        std::vector<double> localRHS;
        if (adv_time_scheme == adv::timeScheme::TTG2)
        {
            
            localRHS = (step_num == StepNumber::first) ? LocalTTG2RhsAssembling(local_node_coords,
                                                                                local_vel_trian,
                                                                                local_scal_trian,
                                                                                local_scal_trian,
                                                                                time_step,
                                                                                1)
                                                       : LocalTTG2RhsAssembling(local_node_coords,
                                                                                local_vel_trian,
                                                                                local_scal_trian,
                                                                                local_scal_half_trian,
                                                                                time_step,
                                                                                1);
        }
        else if (adv_time_scheme == adv::timeScheme::TTG3)
        {
            localRHS = (step_num == StepNumber::first) ? LocalTTG3RhsAssembling(local_node_coords,
                                                                                local_vel_trian,
                                                                                local_scal_trian,
                                                                                local_scal_trian,
                                                                                time_step,
                                                                                1)
                                                       : LocalTTG3RhsAssembling(local_node_coords,
                                                                                local_vel_trian,
                                                                                local_scal_trian,
                                                                                local_scal_half_trian,
                                                                                time_step,
                                                                                1);
        }
        else if (adv_time_scheme == adv::timeScheme::TTG4)
        {
            localRHS = (step_num == StepNumber::first) ? LocalTTG4RhsAssembling(local_node_coords,
                                                                                local_vel_trian,
                                                                                local_scal_trian,
                                                                                local_scal_trian,
                                                                                time_step,
                                                                                1)
                                                       : LocalTTG4RhsAssembling(local_node_coords,
                                                                                local_vel_trian,
                                                                                local_scal_trian,
                                                                                local_scal_half_trian,
                                                                                time_step,
                                                                                1);
        }
        else
        {
            SIMUG_ERR("Only avalible double step solvers: TTG2, TTG3, TTG4");
        }
            
        // assemble global RHS vector
        for (int i = 0; i < 3; ++i)
        {
            if(adj_nodes[i].GetStatus() != Element::Ghost)
            {
                RHS[adj_nodes[i].Integer(node_id)] += localRHS[i];
            }
        }
    }
    BARRIER
}

void AgridAdvectionSolver::Evaluate(velocity_tag vel_tag, scalar_tag scal_tag)
{
    // timer and logger initialization
    double duration;
    SIMUG::Timer adv_timer;

    // create temporarily tags for scalar
    INMOST::Tag scal_low_tag, scal_half_tag, scal_high_tag;

    mesh->GetDataSingle(mesh::gridElemType::Node)->Create("scal_high_tag", 1, INMOST::DATA_REAL);
    scal_high_tag = mesh->GetDataSingle(mesh::gridElemType::Node)->Get("scal_high_tag");

    if (adv_filter == adv::advFilter::Zalesak)
    {
        mesh->GetDataSingle(mesh::gridElemType::Node)->Create("scal_low_tag", 1, INMOST::DATA_REAL);
        scal_low_tag = mesh->GetDataSingle(mesh::gridElemType::Node)->Get("scal_low_tag");
    }

    if (!adv::is_single_step.at(adv_time_scheme))
    {
        mesh->GetDataSingle(mesh::gridElemType::Node)->Create("scal_half_tag", 1, INMOST::DATA_REAL);
        scal_half_tag = mesh->GetDataSingle(mesh::gridElemType::Node)->Get("scal_half_tag");
    }

    // node id tag
    INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

    // create sparse vector for solution
    adv_timer.Launch();

    int idmin = mesh->GetMeshInfo().id_interval_nodes.id_min;
    int idmax = mesh->GetMeshInfo().id_interval_nodes.id_max;
    INMOST::Sparse::Vector SOL;
    SOL.SetInterval(idmin, idmax);

    // make partition for RHS vector
    RHS.SetInterval(idmin, idmax);

    // assemble RHS vector
    (!adv::is_single_step.at(adv_time_scheme)) ? AssembleDoubleStepRHS(vel_tag, scal_tag, scal_tag, StepNumber::first)
                                               : AssembleSingleStepRHS(vel_tag, scal_tag);
    BARRIER
    adv_timer.Stop();
    duration = adv_timer.GetMaxTime();
    adv_timer.Reset();
    RHS_assembling_time += duration;



    // Solve lianear system
    adv_timer.Launch();
    slae_solver->Solve(RHS, SOL);
    adv_timer.Stop();
    duration = adv_timer.GetMaxTime();
    matrix_invertion_time += duration;
    adv_timer.Reset();
    BARRIER

    // update scal high and (scal half optional) value
    for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
            (!adv::is_single_step.at(adv_time_scheme)) ? nodeit->Real(scal_half_tag) = SOL[nodeit->Integer(node_id_tag)] 
                                                       : nodeit->Real(scal_high_tag) = SOL[nodeit->Integer(node_id_tag)];
    }
    BARRIER

    // excchange scal high (scal half)
    (!adv::is_single_step.at(adv_time_scheme)) ? mesh->GetDataSingle(mesh::gridElemType::Node)->Exchange("scal_half_tag")
                                               : mesh->GetDataSingle(mesh::gridElemType::Node)->Exchange("scal_high_tag");
    BARRIER
    
    // make second step if TTG2, TTG3 or TTG4 is used
    if (!adv::is_single_step.at(adv_time_scheme))
    {
        // reset RHS and SOL vector
        adv_timer.Launch();
        SOL.Clear();
        RHS.Clear();
        SOL.SetInterval(idmin, idmax);
        RHS.SetInterval(idmin, idmax);
        BARRIER

        // assemble new RHS vector
        AssembleDoubleStepRHS(vel_tag, scal_tag, scal_half_tag, StepNumber::second);
        BARRIER

        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        RHS_assembling_time += duration;
        adv_timer.Reset();
        
        adv_timer.Launch();
        // solve new slae
        slae_solver->Solve(RHS, SOL);
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        matrix_invertion_time += duration;
        adv_timer.Reset();

        // finally update m_high value
        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
                nodeit->Real(scal_high_tag) = SOL[nodeit->Integer(node_id_tag)];
        }
        BARRIER

        // exchange scal high tag
        mesh->GetDataSingle(mesh::gridElemType::Node)->Exchange("scal_high_tag");
    }

    // Reset SOL and RHS vector
    SOL.Clear();
    RHS.Clear();
    SOL.SetInterval(idmin, idmax);
    RHS.SetInterval(idmin, idmax);
    BARRIER

    /*
    
    // Calculate low order solution if FCT is used
    if (params.is_fct)
    {
        // make old MASS vector and tmp vectors
        Sparse::Vector MASS, tmp_vec_1, tmp_vec_2, tmp_vec_3;
        MASS.SetInterval(idmin, idmax);
        tmp_vec_1.SetInterval(idmin, idmax);
        tmp_vec_2.SetInterval(idmin, idmax);
        tmp_vec_3.SetInterval(idmin, idmax);

        // fill mass vector
        for( Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                MASS[nodeit->Integer(id)] = nodeit->Real(m_tag);
            }
        }
        BARRIER
        mat_mult_vec(LHS, MASS, tmp_vec_1);
        BARRIER
        mat_mult_vec(LHS_low, MASS, tmp_vec_2);
        BARRIER
        vec_mult_num(tmp_vec_1, (params.fct_cd - 1.0), idmin, idmax);
        BARRIER
        vec_mult_num(tmp_vec_2, (1.0 - params.fct_cd), idmin, idmax);
        BARRIER
        vec_plus_vec(tmp_vec_1, tmp_vec_2, tmp_vec_3, idmin, idmax);
        BARRIER
        vec_plus_vec(tmp_vec_3, RHS, RHS_low, idmin, idmax);
        BARRIER


        // solve low order system 
        Sol.SetMatrix(LHS_low);
        Sol.Solve(RHS_low, SOL);

        // update m_low value
        for( Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->Real(m_low_tag) = SOL[nodeit->Integer(id)];
            }
        }
        BARRIER
        n.GetMesh()->ExchangeData(m_low_tag, NODE, 0);

        // clear RHS_low vector
        RHS_low.Clear();
        RHS_low.SetInterval(idmin, idmax); 
    }
    

    // clear RHS and SOL vector
    SOL.Clear();
    RHS.Clear();
    */

    /*
    // apply FCT (if used) and update m finally
    if (params.is_fct)
    {
        FCT_Procedure(n, M_C_minus_M_L, m_tag, m_low_tag, m_high_tag, params.fct_cd);
    }
    else
    {
       for( Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->Real(m_tag) = nodeit->Real(m_high_tag);
            }
        }
        n.GetMesh()->ExchangeData(m_tag, NODE, 0); 
    }
    BARRIER

    
    // clear RHS and SOL vector
    SOL.Clear();
    RHS.Clear();
    SOL.SetInterval(idmin, idmax);
    RHS.SetInterval(idmin, idmax);

*/
    // update scalar value
    for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
            nodeit->Real(scal_tag) = nodeit->Real(scal_high_tag);
    }
    BARRIER

    // exchange scal
    mesh->GetMesh()->ExchangeData(scal_tag, NODE, 0); 

    // clear all temporal tags
    mesh->GetDataSingle(mesh::gridElemType::Node)->Delete("scal_high_tag");
    if (adv_filter == adv::advFilter::Zalesak)
        mesh->GetDataSingle(mesh::gridElemType::Node)->Delete("scal_low_tag");

    if (!adv::is_single_step.at(adv_time_scheme))
         mesh->GetDataSingle(mesh::gridElemType::Node)->Delete("scal_half_tag");
    BARRIER
}

void AgridAdvectionSolver::PrintProfiling()
{
    SIMUG::Logger adv_log(std::cout);

    if (mesh->GetMesh()->GetProcessorRank() == 0)
    {
        adv_log.Log("## Profiling info ##\n");
        adv_log.Log("Total RHS assembling time: " + std::to_string(RHS_assembling_time) + " ms\n");
        adv_log.Log("Total matrix invertion time: " + std::to_string(matrix_invertion_time) + " ms\n");
        adv_log.Log("Total limiter time: " + std::to_string(limiter_time) + " ms\n");
    }
    RHS_assembling_time = 0.0;
    limiter_time = 0.0;
    matrix_invertion_time = 0.0;
    BARRIER
}