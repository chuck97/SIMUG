#include "advection.hpp"

using namespace INMOST;

namespace SIMUG
{

CgridAdvectionSolver::CgridAdvectionSolver(SIMUG::IceMesh *mesh_,
                                           double time_step_,
                                           velocity_tag vel_tag_,
                                           SIMUG::adv::timeScheme adv_time_scheme_,
                                           SIMUG::adv::spaceScheme adv_space_scheme_,
                                           SIMUG::adv::advFilter adv_filter_,
                                           const std::vector<double> &params_) : AdvectionSolver(mesh_, time_step_, vel_tag_, adv_time_scheme_, adv_space_scheme_, adv_filter_),
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

    if (mesh->GetMesh()->GetProcessorRank() == 0)
    {
        adv_log.Log("Maximal Courant number: " + std::to_string(max_courant) + " (" + std::to_string(duration) + " ms)\n");
        if (max_courant > 1.0)
            adv_log.Log("Courant number is greater than 1.0 - advection scheme may be unstable!\n");
    }
    BARRIER

    // create all auxilary tags anf perform preprocessing
    triangle_rhs_tag = mesh->GetMesh()->CreateTag("triangle_rhs_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
    temp_scal_tag = mesh->GetMesh()->CreateTag("temp_scal_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
    //mesh->GetMesh()->SetFileOption("Tag:triangle_rhs_tag", "nosave");
    mesh->GetMesh()->SetFileOption("Tag:temp_scal_tag", "nosave");

    if (adv_space_scheme == adv::spaceScheme::MUST)
    {
        node_scal_tag = mesh->GetMesh()->CreateTag("node_scal_tag", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
        trian_rev_dist_tags = mesh->GetMesh()->CreateTag("trian_rev_dist_tags", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        node_sum_rev_dist_tag = mesh->GetMesh()->CreateTag("node_sum_rev_dist_tag", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
        opposite_node_for_edge_tags = mesh->GetMesh()->CreateTag("opposite_node_for_edge_tags", INMOST::DATA_INTEGER, INMOST::CELL, INMOST::NONE, 3);
        phi_tags = mesh->GetMesh()->CreateTag("phi_tags", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);

        mesh->GetMesh()->SetFileOption("Tag:node_scal_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:trian_rev_dist_tags", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:node_sum_rev_dist_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:opposite_node_for_edge_tags", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:phi_tags", "nosave");

        // compute triangles reversed distances
        adv_timer.Launch();
        ComputeTrianDistances();
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        adv_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            adv_log.Log("Computation of triangle distances: OK (" + std::to_string(duration) + " ms)\n");
        }

        // compute oposite nodes for edges in trians
        adv_timer.Launch();
        ComputeOppositeNodes();
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        adv_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            adv_log.Log("Computation of oposite nodes for edges in every trian: OK (" + std::to_string(duration) + " ms)\n");
        }
    }
    else if (adv_space_scheme == adv::spaceScheme::MUSCL)
    {
        node_scal_tag = mesh->GetMesh()->CreateTag("node_scal_tag", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
        trian_rev_dist_tags = mesh->GetMesh()->CreateTag("trian_rev_dist_tags", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        node_sum_rev_dist_tag = mesh->GetMesh()->CreateTag("node_sum_rev_dist_tag", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
        gradient_trian_tag = mesh->GetMesh()->CreateTag("gradient_trian_tag", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 2);
        edge_distance_vector_tags = mesh->GetMesh()->CreateTag("edge_distance_vector_tags", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 6);
        adj_trian_baric_dist_vec_tags = mesh->GetMesh()->CreateTag("adj_trian_baric_dist_vec_tags", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 6);
        r_factor_tags = mesh->GetMesh()->CreateTag("r_factor_tags", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        phi_tags = mesh->GetMesh()->CreateTag("phi_tags", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
        trian_sc_diff = mesh->GetMesh()->CreateTag("trian_sc_diff", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);

        //mesh->GetMesh()->SetFileOption("Tag:node_scal_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:trian_rev_dist_tags", "nosave");
        //mesh->GetMesh()->SetFileOption("Tag:gradient_trian_tag", "nosave");
        mesh->GetMesh()->SetFileOption("Tag:edge_distance_vector_tags", "nosave");
        //mesh->GetMesh()->SetFileOption("Tag:adj_trian_baric_dist_vec_tags", "nosave");
        //mesh->GetMesh()->SetFileOption("Tag:r_factor_tags", "nosave");
        //mesh->GetMesh()->SetFileOption("Tag:phi_tags", "nosave");
        //mesh->GetMesh()->SetFileOption("Tag:trian_sc_diff", "nosave");

        // compute triangles reversed distances
        adv_timer.Launch();
        ComputeTrianDistances();
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        adv_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            adv_log.Log("Computation of triangle distances: OK (" + std::to_string(duration) + " ms)\n");
        }

        // compute distance vectors from baricenter to edges for every trian
        adv_timer.Launch();
        ComputeEdgeDistanceVectors();
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        adv_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            adv_log.Log("Computation of distance vectors from baricenter to edge center for triangles: OK (" + std::to_string(duration) + " ms)\n");
        }

        // compute adj trians baricenter distances
        adv_timer.Launch();
        ComputeAdjTrianBaricenterDistanceVectors();
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        adv_timer.Reset();

        if (mesh->GetMesh()->GetProcessorRank() == 0)
        {
            adv_log.Log("Computation of adjacent triangles baricenter distance vectors: OK (" + std::to_string(duration) + " ms)\n");
        }
    }
    else
    {
        if (adv_space_scheme != adv::spaceScheme::FVupwind)
        {
            SIMUG_ERR("unknown type of space scheme for C grid");
        }
    }

    if (mesh->GetMesh()->GetProcessorRank() == 0)
        adv_log.Log("=====================================================================\n");
    BARRIER
}

double CgridAdvectionSolver::GetMaxCourant()
{
    // get cart coord tag for node
    INMOST::Tag node_cart_coords_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->coords[coord::coordType::cart];

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
                    adj_nodes[0].RealArray(node_cart_coords_tag)[2]};

            std::vector<double> first_node_coords =
                {
                    adj_nodes[1].RealArray(node_cart_coords_tag)[0],
                    adj_nodes[1].RealArray(node_cart_coords_tag)[1],
                    adj_nodes[1].RealArray(node_cart_coords_tag)[2]};

            std::vector<double> second_node_coords =
                {
                    adj_nodes[2].RealArray(node_cart_coords_tag)[0],
                    adj_nodes[2].RealArray(node_cart_coords_tag)[1],
                    adj_nodes[2].RealArray(node_cart_coords_tag)[2]};

            double a0 = L2_norm_vec(first_node_coords - zero_node_coords);
            double a1 = L2_norm_vec(second_node_coords - zero_node_coords);
            double a2 = L2_norm_vec(second_node_coords - first_node_coords);

            double amin = std::min(std::min(a0, a1), a2);

            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

            double u0 = L2_norm_vec(std::vector<double>{adj_edges[0].RealArray(vel_tag)[0], adj_edges[0].RealArray(vel_tag)[1]});
            double u1 = L2_norm_vec(std::vector<double>{adj_edges[1].RealArray(vel_tag)[0], adj_edges[1].RealArray(vel_tag)[1]});
            double u2 = L2_norm_vec(std::vector<double>{adj_edges[2].RealArray(vel_tag)[0], adj_edges[2].RealArray(vel_tag)[1]});

            double umax = std::max(std::max(u0, u1), u2);
            double courant = umax * time_step / amin;

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
    INMOST::Tag node_cart_coords_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->coords[coord::coordType::cart];

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
                    adj_nodes[0].RealArray(node_cart_coords_tag)[2]};

            std::vector<double> first_node_coords =
                {
                    adj_nodes[1].RealArray(node_cart_coords_tag)[0],
                    adj_nodes[1].RealArray(node_cart_coords_tag)[1],
                    adj_nodes[1].RealArray(node_cart_coords_tag)[2]};

            std::vector<double> second_node_coords =
                {
                    adj_nodes[2].RealArray(node_cart_coords_tag)[0],
                    adj_nodes[2].RealArray(node_cart_coords_tag)[1],
                    adj_nodes[2].RealArray(node_cart_coords_tag)[2]};

            double a0 = L2_norm_vec(first_node_coords - zero_node_coords);
            double a1 = L2_norm_vec(second_node_coords - zero_node_coords);
            double a2 = L2_norm_vec(second_node_coords - first_node_coords);

            double amin = std::min(std::min(a0, a1), a2);

            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

            double u0 = L2_norm_vec(std::vector<double>{adj_edges[0].RealArray(vel_tag)[0], adj_edges[0].RealArray(vel_tag)[1]});
            double u1 = L2_norm_vec(std::vector<double>{adj_edges[1].RealArray(vel_tag)[0], adj_edges[1].RealArray(vel_tag)[1]});
            double u2 = L2_norm_vec(std::vector<double>{adj_edges[2].RealArray(vel_tag)[0], adj_edges[2].RealArray(vel_tag)[1]});

            double umax = std::max(std::max(u0, u1), u2);
            double u_div_dx = umax / amin;

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

void CgridAdvectionSolver::ComputeOppositeNodes()
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
            for (int ed_num = 0; ed_num < 3; ++ed_num)
            {
                int opposite_node_num_for_edge = 0;

                // iterate over nodes and find opposite edge
                for (int node_num = 0; node_num < 3; ++node_num)
                {
                    if ((adj_nodes[node_num].Integer(node_id_tag) != adj_edges[ed_num].getNodes()[0].Integer(node_id_tag)) &&
                        (adj_nodes[node_num].Integer(node_id_tag) != adj_edges[ed_num].getNodes()[1].Integer(node_id_tag)))
                    {
                        opposite_node_num_for_edge = node_num;
                        break;
                    }
                }

                // store opposite node num
                trianit->IntegerArray(opposite_node_for_edge_tags)[ed_num] = opposite_node_num_for_edge;
            }
        }
    }
    mesh->GetMesh()->ExchangeData(opposite_node_for_edge_tags, INMOST::CELL, 0);
    BARRIER
}

void CgridAdvectionSolver::ComputeTrianDistances()
{
    // tag for node coordinates (in trian basis)
    std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();

    for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
    {
        // get the adj nodes
        ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

        // get coords of nodes in trian basis
        std::vector<double> node0 =
        {
            trianit->RealArray(node_coords_in_trian_basis_tags[0])[0],
            trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]
        };

        std::vector<double> node1 =
        {
            trianit->RealArray(node_coords_in_trian_basis_tags[1])[0],
            trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]
        };

        std::vector<double> node2 =
        {
            trianit->RealArray(node_coords_in_trian_basis_tags[2])[0],
            trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]
        };

        // compute coordinates of baricenter of trian
        std::vector<double> baricenter = (node0 + node1 + node2) * (1.0 / 3.0);

        // store the dinstances
        trianit->RealArray(trian_rev_dist_tags)[0] = 1.0 / L2_norm_vec(node0 - baricenter);
        trianit->RealArray(trian_rev_dist_tags)[1] = 1.0 / L2_norm_vec(node1 - baricenter);
        trianit->RealArray(trian_rev_dist_tags)[2] = 1.0 / L2_norm_vec(node2 - baricenter);

        // add the part of rev dist sum to correspondent node
        if (adj_nodes[0].GetStatus() != Element::Ghost)
            adj_nodes[0].Real(node_sum_rev_dist_tag) += trianit->RealArray(trian_rev_dist_tags)[0];

        if (adj_nodes[1].GetStatus() != Element::Ghost)
            adj_nodes[1].Real(node_sum_rev_dist_tag) += trianit->RealArray(trian_rev_dist_tags)[1];

        if (adj_nodes[2].GetStatus() != Element::Ghost)
            adj_nodes[2].Real(node_sum_rev_dist_tag) += trianit->RealArray(trian_rev_dist_tags)[2];
    }
    // exchange data
    mesh->GetMesh()->ExchangeData(trian_rev_dist_tags, INMOST::CELL, 0);
    mesh->GetMesh()->ExchangeData(node_sum_rev_dist_tag, INMOST::NODE, 0);
    BARRIER
}

void CgridAdvectionSolver::InterpolateScalarNodes(INMOST::Tag trian_scalar_tag, INMOST::Tag node_scalar_tag)
{
    // assign zero node scalar
    for (auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
            nodeit->Real(node_scalar_tag) = 0.0;
    }
    mesh->GetMesh()->ExchangeData(node_scalar_tag, INMOST::NODE, 0);

    // compute the sum of scal*rev_dist in every node
    for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
    {
        ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

        for (int i = 0; i < 3; ++i)
        {
            if (adj_nodes[i].GetStatus() != Element::Ghost)
            {
                adj_nodes[i].Real(node_scalar_tag) = adj_nodes[i].Real(node_scalar_tag) +
                                                     trianit->Real(trian_scalar_tag) *
                                                     trianit->RealArray(trian_rev_dist_tags)[i];
            }
        }
    }
    mesh->GetMesh()->ExchangeData(node_scalar_tag, INMOST::NODE, 0);

    // divide node value by the sum of rev_dist
    for (auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
            nodeit->Real(node_scalar_tag) = nodeit->Real(node_scalar_tag) / nodeit->Real(node_sum_rev_dist_tag);
    }
    mesh->GetMesh()->ExchangeData(node_scalar_tag, INMOST::NODE, 0);
    BARRIER
}

void CgridAdvectionSolver::ComputeTrianGradients(INMOST::Tag node_scalar_tag, INMOST::Tag trian_grad_tag)
{
    // iterate over triangles
    for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            std::vector<double> sum_scaled_normal_with_weights = {0.0, 0.0};

            // iterate over edges
            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();

            for (int ed_num = 0; ed_num < 3; ++ed_num)
            {
                // get the normal of the edge in edge basis
                std::vector<double> normal_coords_in_edge_basis = {1.0, 0.0};

                 if (trianit->Integer(mesh->GetGridInfo(mesh::gridElemType::Trian)->GetIsXedgeBasisIsNormal()[ed_num]) == 0)
                        normal_coords_in_edge_basis = normal_coords_in_edge_basis * (-1.0);

                // move vec coords to trian basis
                std::vector<double> normal_coords_in_trian_basis = mesh->VecTransition(normal_coords_in_edge_basis, adj_edges[ed_num], trianit->getCells()[0]);

                // get the length of the edge
                double edge_len = adj_edges[ed_num]->Real(mesh->GetGridInfo(mesh::gridElemType::Edge)->GetCartesianSize());

                // get the scalar values on the edge
                std::vector<double> edge_scals = 
                {
                    adj_edges[ed_num].getNodes()[0].Real(node_scalar_tag),
                    adj_edges[ed_num].getNodes()[1].Real(node_scalar_tag)
                };

                // add the gradient part
                sum_scaled_normal_with_weights = sum_scaled_normal_with_weights + normal_coords_in_trian_basis * 0.5*(edge_scals[0] + edge_scals[1])*edge_len;
            }

            // compute gradient using Gauss theorem
            double trian_area = trianit->Real(mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize());
            std::vector<double> grad_trian = sum_scaled_normal_with_weights*(1.0/trian_area);
            trianit->RealArray(trian_grad_tag)[0] = grad_trian[0]; 
            trianit->RealArray(trian_grad_tag)[1] = grad_trian[1];
        }
    }
    mesh->GetMesh()->ExchangeData(trian_grad_tag, INMOST::CELL, 0);
}

void CgridAdvectionSolver::ComputeEdgeDistanceVectors()
{
    // tag for node coordinates (in trian basis)
    std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();

    // id tag for node
    INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

    // iterate over triangles
    for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            // get node coords of triangle in trian basis
            std::vector<std::vector<double>> node_coords = 
            {
                {trianit->RealArray(node_coords_in_trian_basis_tags[0])[0], trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]},
                {trianit->RealArray(node_coords_in_trian_basis_tags[1])[0], trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]},
                {trianit->RealArray(node_coords_in_trian_basis_tags[2])[0], trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]}
            };

            // compute coordinates of trian baricenter
            std::vector<double> baricenter_coords = (node_coords[0] + node_coords[1] + node_coords[2]) * (1.0 / 3.0);

            // get adjacent edges and nodes for current triangle
            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();
            ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

            // find edge center coordinates and store vector of difference
            for (int ed_num = 0; ed_num < 3; ++ed_num)
            {
                int wrong_node_index = 0;

                for (int i = 0; i < 3; ++i)
                {
                    if ((adj_nodes[i].Integer(node_id_tag) != adj_edges[ed_num].getNodes()[0].Integer(node_id_tag)) &&
                        (adj_nodes[i].Integer(node_id_tag) != adj_edges[ed_num].getNodes()[1].Integer(node_id_tag)))
                    {
                        wrong_node_index = i;
                        break;
                    }
                }

                std::pair<int, int> edge_node_indexes = {0, 0};
                
                if (wrong_node_index == 0)
                {
                    edge_node_indexes = {1 , 2};
                }
                else if (wrong_node_index == 1)
                {
                    edge_node_indexes = {0, 2};
                }
                else
                {
                    edge_node_indexes = {0, 1};
                }

                std::vector<double> edge_center_coords = (node_coords[edge_node_indexes.first] + node_coords[edge_node_indexes.second])*0.5;

                trianit->RealArray(edge_distance_vector_tags)[2*ed_num + 0] = (edge_center_coords - baricenter_coords)[0];
                trianit->RealArray(edge_distance_vector_tags)[2*ed_num + 1] = (edge_center_coords - baricenter_coords)[1];
            }
        }
    }
    mesh->GetMesh()->ExchangeData(edge_distance_vector_tags, INMOST::CELL, 0);
}

double CgridAdvectionSolver::ApplyFilter(adv::advFilter adv_filt, double r_factor)
{
    if (adv_filt == adv::advFilter::none)
    {
        return 1.0;
    }
    else if (adv_filt == adv::advFilter::Minmod)
    {
        return std::max(0.0, std::min(1.0, r_factor));
    }
    else if (adv_filt == adv::advFilter::VanLeer)
    {
        return (r_factor + std::abs(r_factor))/(1.0 + std::abs(r_factor));
    }
    else if (adv_filt == adv::advFilter::Superbee)
    {
        return std::max(std::max(0.0, std::min(1.0, 2.0*r_factor)), std::min(2.0, r_factor));
    }
    else if (adv_filt == adv::advFilter::BarthJesperson)
    {
        return 0.5*(r_factor + 1.0)*std::min(std::min(1.0, (4.0*r_factor)/(r_factor+1.0)), std::min(1.0, 4.0/(r_factor + 1.0)));
    }
    else 
    {
        SIMUG_ERR("unknown type of advection filter!");
    }
    return 1.0;
}

void CgridAdvectionSolver::ComputeAdjTrianBaricenterDistanceVectors()
{
    // get triangle id tag
    INMOST::Tag trian_id_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->id;
    
    // get node id tag
    INMOST::Tag node_id_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->id;

    // get edge id tag
    INMOST::Tag edge_id_tag = mesh->GetGridInfo(mesh::gridElemType::Edge)->id;

    // tag for node coordinates (in trian basis)
    std::vector<INMOST::Tag> node_coords_in_trian_basis_tags = mesh->GetGridInfo(mesh::gridElemType::Trian)->GetNodeCoordsInTrianBasis();

    // iterate over triangles
    for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            // get node coords of triangle in trian basis
            std::vector<std::vector<double>> node_coords = 
            {
                {trianit->RealArray(node_coords_in_trian_basis_tags[0])[0], trianit->RealArray(node_coords_in_trian_basis_tags[0])[1]},
                {trianit->RealArray(node_coords_in_trian_basis_tags[1])[0], trianit->RealArray(node_coords_in_trian_basis_tags[1])[1]},
                {trianit->RealArray(node_coords_in_trian_basis_tags[2])[0], trianit->RealArray(node_coords_in_trian_basis_tags[2])[1]}
            };

            // compute baricenter coords in trian basis
            std::vector<double> baricenter_coords = (node_coords[0] + node_coords[1] + node_coords[2])*(1.0/3.0);

            // get edges and nodes for current triangle
            ElementArray<INMOST::Face> adj_edges = trianit->getFaces();
            ElementArray<INMOST::Node> adj_nodes = trianit->getNodes();

            // iterate over adj edges
            for (int ed_num = 0; ed_num < 3; ++ ed_num)
            {
                // find current edge center coords in trian basis
                int wrong_node_index = 0;

                for (int i = 0; i < 3; ++i)
                {
                    if ((adj_nodes[i].Integer(node_id_tag) != adj_edges[ed_num].getNodes()[0].Integer(node_id_tag)) &&
                        (adj_nodes[i].Integer(node_id_tag) != adj_edges[ed_num].getNodes()[1].Integer(node_id_tag)))
                    {
                        wrong_node_index = i;
                        break;
                    }
                }

                std::pair<int, int> edge_node_indexes = {0, 0};
                
                if (wrong_node_index == 0)
                {
                    edge_node_indexes = {1, 2};
                }
                else if (wrong_node_index == 1)
                {
                    edge_node_indexes = {0, 2};
                }
                else
                {
                    edge_node_indexes = {0, 1};
                }

                std::vector<double> edge_center_coords = (node_coords[edge_node_indexes.first] + node_coords[edge_node_indexes.second])*0.5;

                // assemble vector from edge center to baricenter in trian basis
                std::vector<double> vec_edge_to_bari_trian = baricenter_coords - edge_center_coords;

                // move vector coords to edge basis
                std::vector<double> vec_edge_to_bari_edge = mesh->VecTransition(vec_edge_to_bari_trian, trianit->getCells()[0], adj_edges[ed_num]);

                // get adjacent triangle for current edge if it exists
                INMOST::Cell adj_trian;

                if (adj_edges[ed_num].getCells().size() != 2)
                {
                    continue;
                }
                else
                {
                    // get adjacent triangle
                    adj_trian = (trianit->Integer(trian_id_tag) == adj_edges[ed_num].getCells()[0]->Integer(trian_id_tag))
                                ? adj_edges[ed_num].getCells()[1]
                                : adj_edges[ed_num].getCells()[0];
                }

                // get nodes of adjacent triangle
                ElementArray<INMOST::Node> adj_nodes_adj =  adj_trian.getNodes();

                // get the node coordinates of adjacent triangle (in basis of this triangle)
                std::vector<std::vector<double>> node_coords_adj_trian = 
                {
                    {adj_trian.RealArray(node_coords_in_trian_basis_tags[0])[0], adj_trian.RealArray(node_coords_in_trian_basis_tags[0])[1]},
                    {adj_trian.RealArray(node_coords_in_trian_basis_tags[1])[0], adj_trian.RealArray(node_coords_in_trian_basis_tags[1])[1]},
                    {adj_trian.RealArray(node_coords_in_trian_basis_tags[2])[0], adj_trian.RealArray(node_coords_in_trian_basis_tags[2])[1]}
                };

                // calculate the baricenter of adjacent triangle
                std::vector<double> baricenter_coords_adj = (node_coords_adj_trian[0] + node_coords_adj_trian[1] + node_coords_adj_trian[2])*(1.0/3.0);

                // find current edge center coords in adj trian basis
                int wrong_node_index_adj = 0;

                for (int i = 0; i < 3; ++i)
                {
                    if ((adj_nodes_adj[i].Integer(node_id_tag) != adj_edges[ed_num].getNodes()[0].Integer(node_id_tag)) &&
                        (adj_nodes_adj[i].Integer(node_id_tag) != adj_edges[ed_num].getNodes()[1].Integer(node_id_tag)))
                    {
                        wrong_node_index_adj = i;
                        break;
                    }
                }

                std::pair<int, int> edge_node_indexes_adj = {0, 0};
                
                if (wrong_node_index_adj == 0)
                {
                    edge_node_indexes_adj = {1 , 2};
                }
                else if (wrong_node_index_adj == 1)
                {
                    edge_node_indexes_adj = {0, 2};
                }
                else
                {
                    edge_node_indexes_adj = {0, 1};
                }

                std::vector<double> edge_center_coords_adj = (node_coords_adj_trian[edge_node_indexes_adj.first] + node_coords_adj_trian[edge_node_indexes_adj.second])*0.5;

                // assemble vector from edge center to baricenter of adj trian in adj trian basis
                std::vector<double> vec_edge_to_bari_trian_adj = baricenter_coords_adj - edge_center_coords_adj;

                // move vector coords to edge basis
                std::vector<double> vec_edge_to_bari_edge_adj = mesh->VecTransition(vec_edge_to_bari_trian_adj, adj_trian, adj_edges[ed_num]);

                // finally compute the baricenter difference vector in edge basis
                std::vector<double> bari_diff_vec = vec_edge_to_bari_edge_adj - vec_edge_to_bari_edge;

                // store computed values
                trianit->RealArray(adj_trian_baric_dist_vec_tags)[2*ed_num + 0] = bari_diff_vec[0];
                trianit->RealArray(adj_trian_baric_dist_vec_tags)[2*ed_num + 1] = bari_diff_vec[1];
            }
        }
    }
    mesh->GetMesh()->ExchangeData(adj_trian_baric_dist_vec_tags, INMOST::CELL, 0);
}


void CgridAdvectionSolver::ComputeMUSCLrfactors(adv::advFilter adv_filt, scalar_tag scal_tag)
{
    if (adv_filt != adv::advFilter::none)
    {
        // get triangle id tag
        INMOST::Tag trian_id_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->id;

        // iterate over triangles
        for (auto trianit = mesh->GetMesh()->BeginCell(); trianit != mesh->GetMesh()->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                // get the gradient vector on trian
                std::vector<double> trian_grad = 
                {
                    trianit->RealArray(gradient_trian_tag)[0],
                    trianit->RealArray(gradient_trian_tag)[1]
                };

                // get edges of trian
                ElementArray<INMOST::Face> adj_edges = trianit->getFaces();
                
                // iterate over edges
                for (int ed_num = 0; ed_num < 3; ++ed_num)
                {
                    // get the baricenter difference vector for current edge
                    std::vector<double> bari_diff_vec = 
                    {
                        trianit->RealArray(adj_trian_baric_dist_vec_tags)[2*ed_num + 0],
                        trianit->RealArray(adj_trian_baric_dist_vec_tags)[2*ed_num + 1]
                    };

                    // move gradient vector to current edge basis
                    std::vector<double> edge_grad = mesh->VecTransition(trian_grad, trianit->getCells()[0], adj_edges[ed_num]);

                    // get the adjacent trian
                    INMOST::Cell adj_trian;

                    if (adj_edges[ed_num].getCells().size() != 2)
                    {
                        continue;
                    }
                    else
                    {
                        // get adjacent triangle
                        adj_trian = (trianit->Integer(trian_id_tag) == adj_edges[ed_num].getCells()[0]->Integer(trian_id_tag))
                                    ? adj_edges[ed_num].getCells()[1]
                                    : adj_edges[ed_num].getCells()[0];
                    }
                    
                    // compute and store r_value for current edge
                    if (std::abs(adj_trian.Real(scal_tag) - trianit->Real(scal_tag)) < 1e-5)
                    {
                        trianit->RealArray(r_factor_tags)[ed_num] = 1000.0;
                    }
                    else
                    {
                        trianit->RealArray(r_factor_tags)[ed_num] = 2.0*(edge_grad*bari_diff_vec)/(adj_trian.Real(scal_tag) - trianit->Real(scal_tag)) - 1.0;
                    }
                    trianit->RealArray(trian_sc_diff)[ed_num] = adj_trian.Real(scal_tag) - trianit->Real(scal_tag);
                }
            }
        }
        mesh->GetMesh()->ExchangeData(r_factor_tags, INMOST::CELL, 0);
        mesh->GetMesh()->ExchangeData(trian_sc_diff, INMOST::CELL, 0);
    }
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
                        normal_edge_velocity_component *= -1.0;

                    // get the length of the edge
                    double edge_len = adj_edges[ed_num]->Real(mesh->GetGridInfo(mesh::gridElemType::Edge)->GetCartesianSize());

                    // compute simple FV upwind rhs for triangle
                    current_trian_rhs += 0.5*normal_edge_velocity_component*edge_len*(trianit->Real(scalar_tag) + adj_trian->Real(scalar_tag)) -
                                         0.5*std::abs(normal_edge_velocity_component)*edge_len*(adj_trian->Real(scalar_tag) - trianit->Real(scalar_tag));
                }
                // store the value of computed rhs
                trianit->Real(triangle_rhs_tag) = current_trian_rhs;
            }
        }
    }
    else if (adv_space_scheme == adv::spaceScheme::MUST)
    {
        // get edge id tag 
        INMOST::Tag edge_id_tag = mesh->GetGridInfo(mesh::gridElemType::Edge)->id;

        // get triangle id tag
        INMOST::Tag trian_id_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->id;

        // interpolate the scalar values on nodes
        InterpolateScalarNodes(scalar_tag, node_scal_tag);

        // iterate over triangles
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
                    // get adj triangle if it exists  
                    INMOST::Cell adj_trian;

                    if (adj_edges[ed_num].getCells().size() != 2)
                    {
                        continue;
                    }
                    else
                    {
                        // get adjacent triangle
                        adj_trian = (trianit->Integer(trian_id_tag) == adj_edges[ed_num].getCells()[0]->Integer(trian_id_tag))
                                     ? adj_edges[ed_num].getCells()[1]
                                     : adj_edges[ed_num].getCells()[0];
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
                        normal_edge_velocity_component *= -1.0;

                    // get the length of the edge
                    double edge_len = adj_edges[ed_num]->Real(mesh->GetGridInfo(mesh::gridElemType::Edge)->GetCartesianSize());

                    // find calc trian according to upwind scheme
                    INMOST::Cell calc_trian = (normal_edge_velocity_component > 0.0) ? trianit->getCells()[0] : adj_trian;
                    
                    // find current edge index in calc trian
                    int edg_index = 0;

                    if (normal_edge_velocity_component > 0.0)
                    {
                        edg_index = ed_num;
                    }
                    else
                    {
                        for (int i = 0; i < 3; ++i)
                        {
                            if (calc_trian.getFaces()[i].Integer(edge_id_tag) == adj_edges[ed_num].Integer(edge_id_tag))
                            {
                                edg_index = i;
                                break;
                            } 
                        }
                    }

                    // get precalculated gradient values
                    std::vector<double> prec_grad = {0.5, 0.5};

                    int opposite_node_num = calc_trian.IntegerArray(opposite_node_for_edge_tags)[edg_index];

                    // assemble scalar difference vector
                    std::vector<double> dscal_vec = 
                    {
                        calc_trian.Real(scalar_tag) - calc_trian.getNodes()[opposite_node_num].Real(node_scal_tag),
                        calc_trian.Real(scalar_tag) - calc_trian.getNodes()[opposite_node_num].Real(node_scal_tag)
                    };

                    double increment = prec_grad*dscal_vec;

                    double phi = 1.0;

                    // compute r-factor if limiter is used and apply filter
                    if (adv_filter != adv::advFilter::none)
                    {
                        // get opposite triangle
                        INMOST::Cell opposite_trian = (normal_edge_velocity_component > 0.0) ? adj_trian : trianit->getCells()[0];

                        double r_factor = 0.0;

                        if (std::abs(opposite_trian->Real(scalar_tag) - calc_trian->Real(scalar_tag)) < 1e-5)
                        {
                            phi = 1.0;
                            calc_trian.RealArray(phi_tags)[ed_num] = phi;
                        }
                        else if (std::abs(increment) < 1e-5)
                        {
                            phi = 1.0;
                            calc_trian.RealArray(phi_tags)[ed_num] = phi;
                        }
                        else
                        {
                            // compute r-factor for current edge
                            r_factor = (opposite_trian->Real(scalar_tag) - calc_trian->Real(scalar_tag))/increment;
                            phi = ApplyFilter(adv_filter, r_factor);
                            calc_trian.RealArray(phi_tags)[ed_num] = phi;
                        }
                    } 
                    
                    // find scalar value at the center of the edge
                    double scal_edge = calc_trian.Real(scalar_tag) + 0.5*phi*increment;

                    // update rhs value according to MUST scheme (Upwind scheme with edge scalar)
                    current_trian_rhs += scal_edge*normal_edge_velocity_component*edge_len;
                }
                // store the value of computed rhs
                trianit->Real(triangle_rhs_tag) = current_trian_rhs;
            }
        }
    }
    else if (adv_space_scheme == adv::spaceScheme::MUSCL)
    {
        // get edge id tag 
        INMOST::Tag edge_id_tag = mesh->GetGridInfo(mesh::gridElemType::Edge)->id;

        // get triangle id tag
        INMOST::Tag trian_id_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->id;

        // interpolate the scalar values on nodes
        InterpolateScalarNodes(scalar_tag, node_scal_tag);

        // compute trian gradients (in trian basis) using Gauss theorem
        ComputeTrianGradients(node_scal_tag, gradient_trian_tag);
        
        // compute MUSCL rfactors
        ComputeMUSCLrfactors(adv_filter, scalar_tag);

        // iterate over triangles
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
                    
                    // get the adjacent triangle if it exists
                    INMOST::Cell adj_trian;

                    if (adj_edges[ed_num].getCells().size() != 2)
                    {
                        continue;
                    }
                    else
                    {
                        // get adjacent triangle
                        adj_trian = (trianit->Integer(trian_id_tag) == adj_edges[ed_num].getCells()[0]->Integer(trian_id_tag))
                                        ? adj_edges[ed_num].getCells()[1]
                                        : adj_edges[ed_num].getCells()[0];
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
                        normal_edge_velocity_component *= -1.0;

                    // get the length of the edge
                    double edge_len = adj_edges[ed_num]->Real(mesh->GetGridInfo(mesh::gridElemType::Edge)->GetCartesianSize());

                    // find calc trian according to upwind scheme
                    INMOST::Cell calc_trian = (normal_edge_velocity_component > 0.0) ? trianit->getCells()[0] : adj_trian;

                    // get the gradient value on calc trian
                    std::vector<double> grad_trian = 
                    {
                        calc_trian.RealArray(gradient_trian_tag)[0],
                        calc_trian.RealArray(gradient_trian_tag)[1]
                    };

                    // find current edge index in calc trian numeration
                    int edg_index = 0;

                    if (normal_edge_velocity_component > 0.0)
                    {
                        edg_index = ed_num;
                    }
                    else
                    {
                        for (int i = 0; i < 3; ++i)
                        {
                            if (calc_trian.getFaces()[i].Integer(edge_id_tag) == adj_edges[ed_num].Integer(edge_id_tag))
                            {
                                edg_index = i;
                                break;
                            } 
                        }
                    }

                    // get the distance vector between baricenter to the current edge in calc trian
                    std::vector<double> vec_to_edge = 
                    {
                        calc_trian.RealArray(edge_distance_vector_tags)[edg_index*2 + 0],
                        calc_trian.RealArray(edge_distance_vector_tags)[edg_index*2 + 1]
                    };

                    double phi = 1.0;

                    // apply filter
                    if (adv_filter != adv::advFilter::none)
                    {
                        // get r value for current edge in calc trian and compute filter value
                        double r_factor = calc_trian.RealArray(r_factor_tags)[edg_index];   
                        phi = ApplyFilter(adv_filter, r_factor);
                        trianit->RealArray(phi_tags)[ed_num] = phi;
                    } 

                    // compute edge scalar value for current edge
                    double edge_scalar = calc_trian.Real(scalar_tag) + 0.5*phi*(grad_trian*vec_to_edge);

                    // add simple upwind flux
                    //current_trian_rhs += 0.5*normal_edge_velocity_component*edge_len*(trianit->Real(scalar_tag) + adj_trian->Real(scalar_tag)) -
                    //                     0.5*std::abs(normal_edge_velocity_component)*edge_len*(adj_trian->Real(scalar_tag) - trianit->Real(scalar_tag));
                    
                    // add antidiffusive flux 
                    //current_trian_rhs += 0.5*phi*normal_edge_velocity_component*edge_len*(adj_trian->Real(scalar_tag) - trianit->Real(scalar_tag));
                    current_trian_rhs += edge_scalar*normal_edge_velocity_component*edge_len;
                }
                // store the value of computed rhs
                trianit->Real(triangle_rhs_tag) = current_trian_rhs;
            }
        }
    }
    else
    {
        SIMUG_ERR("Available advaction space schemes for C grid: FVupwind, MUST, MUSCL!");
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
                trianit->Real(scal_tag) = trianit->Real(scal_tag) - (time_step / trian_area) * trianit->Real(triangle_rhs_tag);
            }
        }
        // exchange scalar value
        mesh->GetMesh()->ExchangeData(scal_tag, INMOST::CELL, 0);

        // update step computation time
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        step_computation_time += duration;
        adv_timer.Reset();
    }
    else if (adv_time_scheme == adv::timeScheme::TRK2)
    {
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
                // get area of triangle
                double trian_area = trianit->Real(mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize());

                // compute temp scalar value
                trianit->Real(temp_scal_tag) = trianit->Real(scal_tag) - (time_step / (2.0 * trian_area)) * trianit->Real(triangle_rhs_tag);
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
                // get area of triangle
                double trian_area = trianit->Real(mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize());

                // compute temp scalar value
                trianit->Real(scal_tag) = trianit->Real(scal_tag) - (time_step / trian_area) * trianit->Real(triangle_rhs_tag);
            }
        }
        // exchange scalar value
        mesh->GetMesh()->ExchangeData(scal_tag, INMOST::CELL, 0);

        // update step computation time
        adv_timer.Stop();
        duration = adv_timer.GetMaxTime();
        step_computation_time += duration;
        adv_timer.Reset();
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
        adv_log.Log("## Advection profiling ##\n");
        adv_log.Log("Flux assembling time: " + std::to_string(flux_computation_time) + " ms\n");
        adv_log.Log("Step computation time: " + std::to_string(step_computation_time) + " ms\n");
    }
    flux_computation_time = 0.0;
    step_computation_time = 0.0;
    BARRIER
}

}