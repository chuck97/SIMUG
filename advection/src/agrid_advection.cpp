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
    // get min and max node id
    int idmin = mesh->GetMeshInfo().id_interval_nodes.id_min;
    int idmax = mesh->GetMeshInfo().id_interval_nodes.id_max;

    // setup LHS and RHS partition
    LHS.SetInterval(idmin, idmax);
    RHS.SetInterval(idmin, idmax);

    // setup LHS_low and RHS_low for Zalesak filter
    if (adv_filter == adv::advFilter::Zalesak)
    {
        LHS_low.SetInterval(idmin, idmax);
        RHS_low.SetInterval(idmin, idmax);
    }

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
    }

    if (max_courant > 1.0)
    {
        if (mesh->GetMesh()->GetProcessorRank()==0)
        {
            adv_log.Log("Courant number is greater than 1.0 - advection scheme may be unstable!\n");
        }
    }
    
    /*
    // assemble LHS
    adv_timer.Launch();
    AssembleLHS();
    adv_timer.Stop();
    duration = adv_timer.GetMaxTime();
    adv_timer.Reset();

    if (mesh->GetMesh()->GetProcessorRank()==0)
    {
        adv_log.Log("LHS for advection assembled successfully! (" + to_string(duration) + " ms)\n");
    }
    */
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

            double amax = std::max(std::max(a0, a1), a2);

            double u0 = L2_norm_vec(std::vector<double>{adj_nodes[0].RealArray(vel_tag)[0], adj_nodes[0].RealArray(vel_tag)[1]});
            double u1 = L2_norm_vec(std::vector<double>{adj_nodes[1].RealArray(vel_tag)[0], adj_nodes[1].RealArray(vel_tag)[1]});
            double u2 = L2_norm_vec(std::vector<double>{adj_nodes[2].RealArray(vel_tag)[0], adj_nodes[2].RealArray(vel_tag)[1]});

            double umax = std::max(std::max(u0, u1), u2);
            double courant = umax*time_step/amax;

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

            double amax = std::max(std::max(a0, a1), a2);

            double u0 = L2_norm_vec(std::vector<double>{adj_nodes[0].RealArray(vel_tag)[0], adj_nodes[0].RealArray(vel_tag)[1]});
            double u1 = L2_norm_vec(std::vector<double>{adj_nodes[1].RealArray(vel_tag)[0], adj_nodes[1].RealArray(vel_tag)[1]});
            double u2 = L2_norm_vec(std::vector<double>{adj_nodes[2].RealArray(vel_tag)[0], adj_nodes[2].RealArray(vel_tag)[1]});

            double umax = std::max(std::max(u0, u1), u2);
            double u_div_dx = umax/amax;

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

/*
void AgridAdvectionSolver::AssembleLHS()
{
    // global node id tag
    INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;
    
    // cart coords of node tag
    INMOST::Tag node_cart_coords_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->coords[coord::coordType::cart];

    // reate tag for lumped mass matrix entire
    mesh->GetForcData(mesh::gridElemType::Node)->Create("M_L_entire", 1, INMOST::INMOST_DATA_REAL);
    mesh->GetForcData(mesh::gridElemType::Node)->Mute("M_L_entire");
    INMOST::Tag M_L_entire = mesh->GetForcData(mesh::gridElemType::Node)->Get("M_L_entire");

    for(auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit) 
    {
        // get node coords of trian
        INMOST::ElementArray<Node> adj_nodes = trianit->getNodes();

        double v0_x = trianit->RealArray(node_cart_coords_tag)[0];
        double v0_y = trianit->RealArray(node_cart_coords_tag)[1];
        double v1_x = trianit->RealArray(node_cart_coords_tag)[2];
        double v1_y = trianit->RealArray(node_cart_coords_tag)[3];
        double v2_x = trianit->RealArray(node_cart_coords_tag)[4];
        double v2_y = trianit->RealArray(node_cart_coords_tag)[5];

        std::vector<std::vector<double>> local_node_coords = {{v0_x, v0_y},
                                                              {v1_x, v1_y},
                                                              {v2_x, v2_y}};
        // find global ids of local nodes
        std::vector<int> nums_global = { local_nodes[0].Integer(id),
                                         local_nodes[1].Integer(id),
                                         local_nodes[2].Integer(id) }; 
        // assemble local 3x3 stiffness matrix
        std::vector<std::vector<double>> localLHS = LocalStiffnessMatrixAssembling(local_node_coords);
        if (params.is_fct)
        {
            if (trianit->GetStatus()!= Element::Ghost)
            {
                M_C_minus_M_L.push_back(localLHS);
            } 
        }
        // assemble global LHS mass matrix
        for (int i = 0; i < 3; ++i)
        {
            if(local_nodes[i].GetStatus() != Element::Ghost)
            {
                for (int j = 0; j < 3; ++j)
                {
                    LHS[nums_global[i]][nums_global[j]] += localLHS[i][j];
                }
            }
        }
        if (params.is_fct)
        {
            for (unsigned int i = 0; i < 3; ++i)
            {
                if(local_nodes[i].GetStatus() != Element::Ghost)
                {
                    for (unsigned int j = 0; j < 3; ++j)
                    {
                        LHS_low[nums_global[i]][nums_global[i]] += localLHS[i][j];
                        local_nodes[i]->Real(M_L_entire) += localLHS[i][j];
                    }
                }
            }
            if (trianit->GetStatus() != Element::Ghost)
            {
                for (unsigned int i = 0; i < 3; ++i)
                {
                    for (unsigned int j = 0; j < 3; ++j)
                    {
                        M_C_minus_M_L.back()[i][i] -= localLHS[i][j];
                    }
                }
            }
        } 
    }   
    BARRIER
    if (params.is_fct)
    {
        n.GetMesh()->ExchangeData(M_L_entire, NODE, 0);
        BARRIER
    }

    // multiply MC_minus_ML by ML^-1 for FCT
    if (params.is_fct)
    {
        int t_num = 0;
        for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
                trianit != n.GetMesh()->EndCell();
                ++trianit) 
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                INMOST::ElementArray<Node> local_nodes = trianit->getNodes(); 
                for (unsigned int i = 0; i < 3; ++i)
                {
                    for (unsigned int j = 0; j < 3; ++j)
                    {
                        M_C_minus_M_L[t_num][i][j] /= - local_nodes[i]->Real(M_L_entire);
                    }
                }
                ++t_num;
            }
        }
    }
    BARRIER
    n.GetMesh()->DeleteTag(M_L_entire, NODE);
}
*/