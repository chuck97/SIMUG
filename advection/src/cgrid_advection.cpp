#include "advection.hpp"

using namespace INMOST;
using namespace SIMUG;

CgridAdvectionSolver::CgridAdvectionSolver(SIMUG::IceMesh* mesh_,
                                           double time_step_,
                                           velocity_tag vel_tag_,
                                           SIMUG::adv::timeScheme adv_time_scheme_,
                                           SIMUG::adv::spaceScheme adv_space_scheme_,
                                           SIMUG::adv::advFilter adv_filter_,
                                           const std::vector<double>& params_):
            AdvectionSolver(mesh_, time_step_, vel_tag_, adv_time_scheme_, adv_space_scheme_, adv_filter_),
            params(params_)
{
    // initialize timer and logger
    SIMUG::Logger adv_log(std::cout);
    SIMUG::Timer adv_timer;
    double duration;

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
    BARRIER

    // create tag for rhs on triangles
    triangle_rhs_tag = mesh->GetMesh()->CreateTag("triangle_rhs_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
    temp_tag = mesh->GetMesh()->CreateTag("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);

    if (mesh->GetMesh()->GetProcessorRank()==0)
        adv_log.Log("=====================================================================\n");
    BARRIER
}

double CgridAdvectionSolver::GetMaxCourant()
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

            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

            double u0 = L2_norm_vec(std::vector<double>{adj_edges[0].RealArray(vel_tag)[0], adj_edges[0].RealArray(vel_tag)[1]});
            double u1 = L2_norm_vec(std::vector<double>{adj_edges[1].RealArray(vel_tag)[0], adj_edges[1].RealArray(vel_tag)[1]});
            double u2 = L2_norm_vec(std::vector<double>{adj_edges[2].RealArray(vel_tag)[0], adj_edges[2].RealArray(vel_tag)[1]});

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

double CgridAdvectionSolver::GetMaxUdivDx()
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

            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

            double u0 = L2_norm_vec(std::vector<double>{adj_edges[0].RealArray(vel_tag)[0], adj_edges[0].RealArray(vel_tag)[1]});
            double u1 = L2_norm_vec(std::vector<double>{adj_edges[1].RealArray(vel_tag)[0], adj_edges[1].RealArray(vel_tag)[1]});
            double u2 = L2_norm_vec(std::vector<double>{adj_edges[2].RealArray(vel_tag)[0], adj_edges[2].RealArray(vel_tag)[1]});

            double umax = std::max(std::max(u0, u1), u2);
            double u_div_dx = umax/amin;

            if (u_div_dx > max_u_div_dx)
                max_u_div_dx = u_div_dx;
        }
    }
    BARRIER
        
    // get max courant for all processors if MPI is used
    double max_u_div_dx_for_all = max_u_div_dx;

#ifdef USE_MPI
    MPI_Allreduce(&max_u_div_dx, &max_u_div_dx_for_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    return max_u_div_dx_for_all;
}

void CgridAdvectionSolver::ComputeRHS(INMOST::Tag scalar_tag)
{
    if (adv_space_scheme == adv::spaceScheme::FVupwind)
    {
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                // current triangle RHS variable
                double current_trian_rhs = 0.0;

                // iterate over current triangle edges
                ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

                for (size_t ed_num = 0; ed_num < 3; ++ed_num)
                {
                    INMOST::Cell adj_trian;
                    // figure out if adjacent triangle for current edge exists
                    if (adj_edges[ed_num].getCells().size() != 2)
                    {
                        continue;
                    }
                    else
                    {
                        // get adjacent triangle
                        adj_trian = (trianit->Integer(mesh->GetGridInfo(mesh::gridElemType::Trian)->id) == 
                                     adj_edges[ed_num].getCells()[0]->Integer(mesh->GetGridInfo(mesh::gridElemType::Trian)->id)) 
                                     ? adj_edges[ed_num].getCells()[1]
                                     : adj_edges[ed_num].getCells()[0];
                        
                        //std::cout << trianit->Integer(mesh->GetGridInfo(mesh::gridElemType::Trian)->id) << " " << adj_trian->Integer(mesh->GetGridInfo(mesh::gridElemType::Trian)->id) << std::endl;
                    }

                    // move velocity components from geo to edge basis
                    std::vector<double> edge_geo_velocity = 
                    {
                        adj_edges[ed_num].RealArray(vel_tag)[0],
                        adj_edges[ed_num].RealArray(vel_tag)[1]
                    };

                    std::vector<double> edge_edge_velocity = mesh->VecTransitionToElemBasis(edge_geo_velocity, adj_edges[ed_num]);

                    // figure out if edge basis x vector collinear to triangle normal and recalculate normal component of velocity
                    double normal_edge_velocity_component = edge_edge_velocity[0];
                    
                    if (trianit->Integer(mesh->GetGridInfo(mesh::gridElemType::Trian)->GetIsXedgeBasisIsNormal()[ed_num]) == 0)
                        normal_edge_velocity_component = normal_edge_velocity_component*(-1.0);

                    // get the length of the edge
                    double edge_len = adj_edges[ed_num]->Real(mesh->GetGridInfo(mesh::gridElemType::Edge)->GetCartesianSize());

                    // compute simple FV upwind rhs for triangle
                    current_trian_rhs += (normal_edge_velocity_component > 0.0) ? trianit->Real(scalar_tag)*normal_edge_velocity_component*edge_len
                                                                                : adj_trian->Real(scalar_tag)*normal_edge_velocity_component*edge_len;
                    
                }
                // store the value of computed rhs
                trianit->Real(triangle_rhs_tag) = current_trian_rhs;
            }
        }
    }
    else
    {
        SIMUG_ERR("only FV upwind space discretization is implemented for triangle C grid yet!");
    }
}

void CgridAdvectionSolver::Evaluate(velocity_tag vel_tag, scalar_tag scal_tag)
{
    // initialize timer and logger
    SIMUG::Logger adv_log(std::cout);
    SIMUG::Timer adv_timer;
    double duration;

    // tag for node coordinates (in trian basis)
    std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();

    // update scalar variable according to time scheme
    if (adv_time_scheme == adv::timeScheme::Euler)
    {
        // calculate rhs for every triangle
        adv_timer.Launch();
        ComputeRHS(scal_tag);
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        flux_computation_time += duration;
        adv_timer.Reset();

        // update scalar according to Euler scheme
        adv_timer.Launch();
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                // get area of triangle
                double trian_area = trianit->Real(mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize());
                
                // compute new scalar value 
                //std::cout << trianit->Real(scal_tag) << " " << trianit->Real(scal_tag) - (time_step/trian_area)*trianit->Real(triangle_rhs_tag) << std::endl;
                trianit->Real(scal_tag) = trianit->Real(scal_tag) - (time_step/trian_area)*trianit->Real(triangle_rhs_tag);
                //trianit->Real(temp_tag) = -(time_step/trian_area)*trianit->Real(triangle_rhs_tag);
            }
        }
        // exchange scalar value
        mesh->GetMesh()->ExchangeData(scal_tag, INMOST::CELL, 0);
        mesh->GetMesh()->ExchangeData(temp_tag, INMOST::CELL, 0);

        // update step computation time
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        step_computation_time += duration;
        adv_timer.Reset();
    }
    else if (adv_time_scheme == adv::timeScheme::TRK2)
    {
        // create tag for temp scalar
        mesh->GetDataSingle(mesh::gridElemType::Trian)->Create("temp scalar", 1, INMOST::DATA_REAL);
        INMOST::Tag temp_scal_tag = mesh->GetDataSingle(mesh::gridElemType::Trian)->Get("temp scalar");
        
        // calculate rhs for every triangle (first step)
        adv_timer.Launch();
        ComputeRHS(scal_tag);
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        flux_computation_time += duration;
        adv_timer.Reset();

        // update scalar according to TRK2 scheme (first step)
        adv_timer.Launch();
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                // compute area of triangle
                std::vector<std::vector<double>> local_node_coords = 
                {
                    {trianit->RealArray(node_coords_in_trian_basis_tags[0])[0], trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]},
                    {trianit->RealArray(node_coords_in_trian_basis_tags[1])[0], trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]},
                    {trianit->RealArray(node_coords_in_trian_basis_tags[2])[0], trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]}
                };

                double trian_area = trian_square(local_node_coords[0], local_node_coords[1], local_node_coords[2]);
                
                // compute temp scalar value 
                trianit->Real(temp_scal_tag) = trianit->Real(scal_tag) - (time_step/(2.0*trian_area))*trianit->Real(triangle_rhs_tag);
            }
        }
        // exchange scalar value
        mesh->GetMesh()->ExchangeData(temp_scal_tag, INMOST::CELL, 0);

        // update step computation time
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        step_computation_time += duration;
        adv_timer.Reset();

        // calculate rhs for every triangle (second step)
        adv_timer.Launch();
        ComputeRHS(temp_scal_tag);
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        flux_computation_time += duration;
        adv_timer.Reset();

        // update scalar according to TRK2 scheme (second step)
        adv_timer.Launch();
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                // compute area of triangle
                std::vector<std::vector<double>> local_node_coords = 
                {
                    {trianit->RealArray(node_coords_in_trian_basis_tags[0])[0], trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]},
                    {trianit->RealArray(node_coords_in_trian_basis_tags[1])[0], trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]},
                    {trianit->RealArray(node_coords_in_trian_basis_tags[2])[0], trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]}
                };

                double trian_area = trian_square(local_node_coords[0], local_node_coords[1], local_node_coords[2]);
                
                // compute temp scalar value 
                trianit->Real(scal_tag) = trianit->Real(temp_scal_tag) - (time_step/trian_area)*trianit->Real(triangle_rhs_tag);
            }
        }
        // exchange scalar value
        mesh->GetMesh()->ExchangeData(scal_tag, INMOST::CELL, 0);

        // update step computation time
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        step_computation_time += duration;
        adv_timer.Reset();

        // delete tag for temp scalar
        mesh->GetDataSingle(mesh::gridElemType::Trian)->Delete("temp scalar");
    }
    else
    {
        SIMUG_ERR("currently available only Euler and 2-step Runge-Kutta of 2nd order time scheme!");
    }
    BARRIER
}

void CgridAdvectionSolver::PrintProfiling()
{
    SIMUG::Logger adv_log(std::cout);

    if (mesh->GetMesh()->GetProcessorRank() == 0)
    {
        adv_log.Log("## Profiling info ##\n");
        adv_log.Log("Total flux assembling time: " + std::to_string(flux_computation_time) + " ms\n");
        adv_log.Log("Total step computation time: " + std::to_string(step_computation_time) + " ms\n");
        adv_log.Log("Total limiter time: " + std::to_string(limiter_time) + " ms\n");
    }
    limiter_time = 0.0;
    flux_computation_time = 0.0;
    step_computation_time = 0.0;
    BARRIER
}