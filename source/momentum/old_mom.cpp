#include "momentum.h"

using namespace std;
using namespace INMOST;



MomentumSolver::MomentumSolver(IceMesh& im,
                               MomentumParams& mp,
                               ModelParams& modp,
                               MeshParams& mep,
                               INMOST::Solver& solver,
                               bool is_verbose_momentum_):
    ice_mesh(im),
    momentum_params(mp),
    model_params(modp),
    mesh_params(mep),
    Sol(solver),
    is_verbose_momentum(is_verbose_momentum_)
{
    // create temporal tags
    u_old_tag = ice_mesh.GetMesh()->CreateTag("u_old", DATA_REAL, NODE, NONE, 2);
    u_n_tag = ice_mesh.GetMesh()->CreateTag("u_n", DATA_REAL, NODE, NONE, 2);
    u_new_tag = ice_mesh.GetMesh()->CreateTag("u_new", DATA_REAL, NODE, NONE, 2);

    sigma_tag = ice_mesh.GetMesh()->CreateTag("sigma", DATA_REAL, CELL, NONE, 3);
    delta_tag = ice_mesh.GetMesh()->CreateTag("delta", DATA_REAL, CELL, NONE, 1);
    P_0_tag = ice_mesh.GetMesh()->CreateTag("P_0", DATA_REAL, CELL, NONE, 1);
    varepsilon_tag = ice_mesh.GetMesh()->CreateTag("varepsilon", DATA_REAL, CELL, NONE, 3);
    sigma_diff_tag = ice_mesh.GetMesh()->CreateTag("sigma_diff", DATA_REAL, CELL, NONE, 3);
    u_diff_tag = ice_mesh.GetMesh()->CreateTag("u_diff", DATA_REAL, NODE, NONE, 3);
    lumped_mass_matrix_entry_tag = ice_mesh.GetMesh()->CreateTag("lumped_mass_matrix_entry_tag", DATA_REAL, NODE, NONE, 1);
    Nx_matrix_entries_tag = ice_mesh.GetMesh()->CreateTag("Nx_matrix_entries_tag", DATA_REAL, CELL, NONE, 3);
    Ny_matrix_entries_tag = ice_mesh.GetMesh()->CreateTag("Ny_matrix_entries_tag", DATA_REAL, CELL, NONE, 3);
    Force_vector_tag = ice_mesh.GetMesh()->CreateTag("Force_vector_tag", DATA_REAL, NODE, NONE, 2);
    Level_vector_tag = ice_mesh.GetMesh()->CreateTag("Level_vector_tag", DATA_REAL, NODE, NONE, 2);
    shear_deformation_tag = ice_mesh.GetMesh()->CreateTag("shear_deformation", DATA_REAL, CELL, NONE, 1);

    // do not print this tags
    ice_mesh.GetMesh()->SetFileOption("Tag:u_old", "nosave");
    ice_mesh.GetMesh()->SetFileOption("Tag:u_n", "nosave");
    ice_mesh.GetMesh()->SetFileOption("Tag:u_new", "nosave");
    ice_mesh.GetMesh()->SetFileOption("Tag:sigma", "nosave");
    //ice_mesh.GetMesh()->SetFileOption("Tag:delta", "nosave");
    //ice_mesh.GetMesh()->SetFileOption("Tag:P_0", "nosave");
    ice_mesh.GetMesh()->SetFileOption("Tag:sigma_diff", "nosave");
    ice_mesh.GetMesh()->SetFileOption("Tag:u_diff", "nosave");
    ice_mesh.GetMesh()->SetFileOption("Tag:lumped_mass_matrix_entry_tag", "nosave");
    ice_mesh.GetMesh()->SetFileOption("Tag:Nx_matrix_entries_tag", "nosave");
    ice_mesh.GetMesh()->SetFileOption("Tag:Ny_matrix_entries_tag", "nosave");
    ice_mesh.GetMesh()->SetFileOption("Tag:Force_vector_tag", "nosave");
    ice_mesh.GetMesh()->SetFileOption("Tag:Level_vector_tag", "nosave");

    BARRIER

    // assemble Mass matrix and N matrix entries
    AssembleLHS();
    LogInit();
}

void MomentumSolver::AssembleLHS()
{
    // no-slip boundary LHS matrix
    if (momentum_params.GetMomentumBC() == MomentumBC::no_slip)
    {
        for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
                trianit != ice_mesh.GetMesh()->EndCell();
                ++trianit) 
        {
            INMOST::ElementArray<Node> local_nodes = trianit->getNodes();

            // first node
            double v0_x_0 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[0];
            double v0_y_0 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[1];
            double v1_x_0 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[2];
            double v1_y_0 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[3];
            double v2_x_0 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[4];
            double v2_y_0 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[5];
            std::vector<std::vector<double>> local_node_coords_0 = {{v0_x_0, v0_y_0},
                                                                    {v1_x_0, v1_y_0},
                                                                    {v2_x_0, v2_y_0}};

            vector<vector<double>> localLHS_0 = LocalStiffnessMatrixAssembling(local_node_coords_0);

            // second node
            double v0_x_1 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[6];
            double v0_y_1 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[7];
            double v1_x_1 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[8];
            double v1_y_1 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[9];
            double v2_x_1 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[10];
            double v2_y_1 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[11];
            std::vector<std::vector<double>> local_node_coords_1 = {{v0_x_1, v0_y_1},
                                                                    {v1_x_1, v1_y_1},
                                                                    {v2_x_1, v2_y_1}};

            vector<vector<double>> localLHS_1 = LocalStiffnessMatrixAssembling(local_node_coords_1);

            // third node
            double v0_x_2 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[12];
            double v0_y_2 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[13];
            double v1_x_2 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[14];
            double v1_y_2 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[15];
            double v2_x_2 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[16];
            double v2_y_2 = trianit->RealArray(ice_mesh.GetLocalBasisData().nodal_data.LocalNodeCoords)[17];
            
            std::vector<std::vector<double>> local_node_coords_2 = {{v0_x_2, v0_y_2},
                                                                    {v1_x_2, v1_y_2},
                                                                    {v2_x_2, v2_y_2}};
            // assemble local 3x3 mass matrix for first node 
            vector<vector<double>> localLHS_2 = LocalStiffnessMatrixAssembling(local_node_coords_2);

            vector<vector<vector<double>>> localLHS = {localLHS_0, localLHS_1, localLHS_2};

            // assemble lumped mass matrix entries
            for (int i = 0; i < 3; ++i)
            {
                if((local_nodes[i].GetStatus() != Element::Ghost) and
                   (local_nodes[i].Integer(ice_mesh.GetData().IsNodeBnd) == 0))
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        if (local_nodes[j].Integer(ice_mesh.GetData().IsNodeBnd) == 0)
                        {
                            local_nodes[i]->Real(lumped_mass_matrix_entry_tag) += localLHS[i][i][j];
                        }
                    }
                }
            }


            // assemble local N vector
            vector<double> local_Nx_0 = LocalNVectorAssembling(local_node_coords_0, 0);
            vector<double> local_Ny_0 = LocalNVectorAssembling(local_node_coords_0, 1);

            vector<double> local_Nx_1 = LocalNVectorAssembling(local_node_coords_1, 0);
            vector<double> local_Ny_1 = LocalNVectorAssembling(local_node_coords_1, 1);

            vector<double> local_Nx_2 = LocalNVectorAssembling(local_node_coords_2, 0);
            vector<double> local_Ny_2 = LocalNVectorAssembling(local_node_coords_2, 1);

            vector<vector<double>> local_Nx = {local_Nx_0, local_Nx_1, local_Nx_2};
            vector<vector<double>> local_Ny = {local_Ny_0, local_Ny_1, local_Ny_2};
            
            // assemble global N matrix etries
            for (int i = 0; i < 3; ++i)
            {
                if(local_nodes[i].Integer(ice_mesh.GetData().IsNodeBnd) == 0)
                {
                    trianit->RealArray(Nx_matrix_entries_tag)[i] = local_Nx[i][i];
                    trianit->RealArray(Ny_matrix_entries_tag)[i] = local_Ny[i][i];
                }
            }
        }
        BARRIER
        ice_mesh.GetMesh()->ExchangeData(lumped_mass_matrix_entry_tag, NODE, 0);
        BARRIER
        ice_mesh.GetMesh()->ExchangeData(Nx_matrix_entries_tag, CELL, 0);
        BARRIER
        ice_mesh.GetMesh()->ExchangeData(Ny_matrix_entries_tag, CELL, 0);
        BARRIER

        // assemble local 3x3 N matrix
    }
    else
    {
        INMOST_ICE_ERR("current momentum BC is not available")
    }    
};

void MomentumSolver::LogInit()
{
    if (ice_mesh.GetMesh()->GetProcessorRank() == 0)
    {
        cout << "=============================================" << endl;
        cout << "Momentum Solver Info" << endl;
        cout << "Momentum solver type: " << MomentumSolverTypeToName[momentum_params.GetMomentumSolverType()] << endl;
        if (momentum_params.GetMomentumSolverType() == MomentumSolverType::mEVP)
        {
            cout << "alpha: " << momentum_params.GetSolverParams()["alpha"].Double << endl;
            cout << "beta: " << momentum_params.GetSolverParams()["beta"].Double << endl;
            if (momentum_params.GetSolverParams()["number of iterations"].Int != 0)
            {
                cout << "number of iterations: " << momentum_params.GetSolverParams()["number of iterations"].Int << endl;
            }
            else
            {
                cout << "relative residual: " << momentum_params.GetSolverParams()["relative residual"].Double << endl;
            }   
        }
        cout << "Is Coriolis: " << momentum_params.GetIsCoriolis() << endl;
        cout << "Is water drag: " << momentum_params.GetIsWaterDrag() << endl;
        cout << "Is air drag: " << momentum_params.GetIsAirDrag() << endl;
        cout << "Boudary conditions: " << MomentumBCToName[momentum_params.GetMomentumBC()] << endl;
        cout << "=============================================" << endl;
    }
    BARRIER
};

void MomentumSolver::Evaluate()
{
    if (momentum_params.GetMomentumSolverType() == MomentumSolverType::mEVP)
    {
        mEVP_evaluation();
    }
    else
    {
        INMOST_ICE_ERR("currently onle mEVP solver is available");
    }
};

void MomentumSolver::UpdateH()
{
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                nodeit != ice_mesh.GetMesh()->EndNode();
                ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            double m = nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::m]);
            double a = nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]);
            double a_min = model_params.GetMinConcentration();
            double h_min = model_params.GetMinHeight();
            double rho = model_params.GetIceDensity();
            double h = 0.0;

            if ((a < a_min) or (m < h_min*a_min*rho)) 
            {
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]) = a_min;
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h]) = h_min;
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::m]) = rho*h_min*a_min;
                continue;
            }

            if (a > 1.0)
            {
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]) = 1.0;
                a = 1.0;
            }

            h = m/(rho*a);

            if (h < h_min)
            {
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]) = a_min;
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h]) = h_min;
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::m]) = rho*h_min*a_min;
            }
            else
            {
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h]) = h;
            }
        }
    } 
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::h], NODE, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::m], NODE, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::a], NODE, 0);
    BARRIER
};

void MomentumSolver::UpdateM()
{
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                nodeit != ice_mesh.GetMesh()->EndNode();
                ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            double h = nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h]);
            double a = nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]);
            double a_min = model_params.GetMinConcentration();
            double rho = model_params.GetIceDensity();
            double m = 0.0;

            if (a < a_min)
            {
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]) = 0.0;
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h]) = 0.0;
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::m]) = 0.0;
                continue;
            }

            if (m < 0.0)
            {
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::m]) = 0.0;
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]) = 0.0;
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h]) = 0.0;
                continue;
            }


            if (a > 1.0)
            {
                nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]) = 1.0;
                a = 1.0;
            }

            m = h*rho*a;

            nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::m]) = m;
        }
    } 
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::h], NODE, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::m], NODE, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::a], NODE, 0);
    BARRIER
};

void MomentumSolver::AssignSigma()
{
    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
                trianit != ice_mesh.GetMesh()->EndCell();
                ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            trianit->RealArray(sigma_tag)[0] = trianit->Real(ice_mesh.GetData().TriangleData[ModelVariableNotation::sig1]);
            trianit->RealArray(sigma_tag)[1] = trianit->Real(ice_mesh.GetData().TriangleData[ModelVariableNotation::sig2]);
            trianit->RealArray(sigma_tag)[2] = trianit->Real(ice_mesh.GetData().TriangleData[ModelVariableNotation::sig12]);
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(sigma_tag, CELL, 0);
    BARRIER
};

void MomentumSolver::AssignVelocity()
{
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                nodeit != ice_mesh.GetMesh()->EndNode();
                ++nodeit)
    {
        if ((nodeit->GetStatus() != Element::Ghost) and
            nodeit->Integer(ice_mesh.GetData().IsNodeBnd) == 0)
        {
            nodeit->RealArray(u_old_tag)[0] = nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_ice])[0];
            nodeit->RealArray(u_old_tag)[1] = nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_ice])[1];
            nodeit->RealArray(u_n_tag)[0] = nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_ice])[0];
            nodeit->RealArray(u_n_tag)[1] = nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_ice])[1];
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(u_old_tag, NODE, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(u_n_tag, NODE, 0);
    BARRIER
};

void MomentumSolver::CalculateP_0()
{
    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
                trianit != ice_mesh.GetMesh()->EndCell();
                ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            ElementArray<Node> nodes = trianit->getNodes();
            
            double h = (1/3.0)*(
            nodes[0]->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h]) +
            nodes[1]->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h]) +
            nodes[2]->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h]));

            double a = (1/3.0)*(
            nodes[0]->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]) +
            nodes[1]->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]) +
            nodes[2]->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]));

            double p_str = model_params.GetPressureStar();
            double p_coeff = model_params.GetPressureCoeff();


            trianit->Real(P_0_tag) = p_str*h*exp(-p_coeff*(1.0 - a));
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(P_0_tag, CELL, 0);
    BARRIER
};

void MomentumSolver::CalculateVarepsilonDelta()
{
    double e = model_params.GetEccentricity();
    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
                trianit != ice_mesh.GetMesh()->EndCell();
                ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            // node coords
            ElementArray<Node> local_nodes = trianit->getNodes();
            double v0_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[0];
            double v0_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[1];
            double v1_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[2];
            double v1_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[3];
            double v2_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[4];
            double v2_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[5];
            vector<std::vector<double>> local_node_coords = {{v0_x, v0_y},
                                                             {v1_x, v1_y},
                                                             {v2_x, v2_y}};

            // vec coords
            double u0 = local_nodes[0]->RealArray(u_old_tag)[0];
            double v0 = local_nodes[0]->RealArray(u_old_tag)[1];
            double u1 = local_nodes[1]->RealArray(u_old_tag)[0];
            double v1 = local_nodes[1]->RealArray(u_old_tag)[1];
            double u2 = local_nodes[2]->RealArray(u_old_tag)[0];
            double v2 = local_nodes[2]->RealArray(u_old_tag)[1];

            vector<std::vector<double>> velocity_coords = {{u0, v0},
                                                           {u1, v1},
                                                           {u2, v2}
                                                          };
            
            std::vector<double> u0_local = 
            {
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[0]*velocity_coords[0][0] +
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[1]*velocity_coords[0][1],
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[2]*velocity_coords[0][0] +
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[3]*velocity_coords[0][1]
            };

            std::vector<double> u1_local = 
            {
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[4]*velocity_coords[1][0] +
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[5]*velocity_coords[1][1],
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[6]*velocity_coords[1][0] +
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[7]*velocity_coords[1][1]
            };
            std::vector<double> u2_local = 
            {
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[8]*velocity_coords[2][0] +
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[9]*velocity_coords[2][1],
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[10]*velocity_coords[2][0] +
                trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[11]*velocity_coords[2][1]
            };
            std::vector<std::vector<double>> local_u =  {u0_local, u1_local, u2_local};
            
            ScalarFunction u_loc = {{local_u[0][0], local_u[1][0], local_u[2][0]},
                                    local_node_coords};
            
            ScalarFunction v_loc = {{local_u[0][1], local_u[1][1], local_u[2][1]},
                                    local_node_coords};

            ScalarFunction du_dx_local = d_dx(u_loc);
            ScalarFunction du_dy_local = d_dy(u_loc);
            ScalarFunction dv_dx_local = d_dx(v_loc);
            ScalarFunction dv_dy_local = d_dy(v_loc);

            double du_dx = du_dx_local.GetParams()[0];
            double du_dy = du_dy_local.GetParams()[0];
            double dv_dx = dv_dx_local.GetParams()[0];
            double dv_dy = dv_dy_local.GetParams()[0];

            double eps11 = du_dx;
            double eps22 = dv_dy;
            double eps12 = 0.5*(du_dy + dv_dx);

            double delta = sqrt( (eps11*eps11 + eps22*eps22)*(1.0 + 1.0/(e*e)) +
                                  eps12*eps12*(4.0/(e*e)) +
                                  eps11*eps22*2.0*(1.0 - 1.0/(e*e))
                               );

            trianit->RealArray(varepsilon_tag)[0] = eps11 + eps22;
            trianit->RealArray(varepsilon_tag)[1] = eps11 - eps22;
            trianit->RealArray(varepsilon_tag)[2] = eps12;
            trianit->Real(delta_tag) = delta;
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(varepsilon_tag, CELL, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(delta_tag, CELL, 0);
    BARRIER
};

void MomentumSolver::UpdateTemporalSigma()
{ 
    double e = model_params.GetEccentricity();
    double delta_min = model_params.GetDeltaMin();
    double alpha = momentum_params.GetSolverParams()["alpha"].Double;

    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
                trianit != ice_mesh.GetMesh()->EndCell();
                ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            double P_0 = trianit->Real(P_0_tag);
            double delta = trianit->Real(delta_tag);
            double sig1 = trianit->RealArray(sigma_tag)[0];
            double sig2 = trianit->RealArray(sigma_tag)[1];
            double sig12 = trianit->RealArray(sigma_tag)[2];

            double eps1 = trianit->RealArray(varepsilon_tag)[0];
            double eps2 = trianit->RealArray(varepsilon_tag)[1];
            double eps12 = trianit->RealArray(varepsilon_tag)[2];

            trianit->RealArray(sigma_tag)[0] = (1.0 - 1.0/alpha)*sig1 + 
            (1.0/alpha)*(P_0/(delta + delta_min))*(eps1 - delta);
            trianit->RealArray(sigma_tag)[1] = (1.0 - 1.0/alpha)*sig2 + 
            (1.0/(alpha*e*e))*(P_0/(delta + delta_min))*(eps2);
            trianit->RealArray(sigma_tag)[2] = (1.0 - 1.0/alpha)*sig12 + 
            (1.0/(alpha*e*e))*(P_0/(delta + delta_min))*(eps12);

            // calculate sigma diff
            trianit->RealArray(sigma_diff_tag)[0] = trianit->RealArray(sigma_tag)[0] - sig1;
            trianit->RealArray(sigma_diff_tag)[1] = trianit->RealArray(sigma_tag)[1] - sig2;
            trianit->RealArray(sigma_diff_tag)[2] = trianit->RealArray(sigma_tag)[2] - sig12;
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(sigma_tag, CELL, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(sigma_diff_tag, CELL, 0);
    BARRIER
};

void MomentumSolver::AssembleLevelVector()
{
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                nodeit != ice_mesh.GetMesh()->EndNode();
                ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            nodeit->RealArray(Level_vector_tag)[0] = -model_params.GetCoriolisParam()*
                                                      nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_water])[1]/
                                                      model_params.GetGravity();
            nodeit->RealArray(Level_vector_tag)[1] =  model_params.GetCoriolisParam()*
                                                      nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_water])[0]/
                                                      model_params.GetGravity();
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(Level_vector_tag, NODE, 0);
    BARRIER
};

void MomentumSolver::AssembleForceVector()
{
    // clean Force_vector_tag
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                nodeit != ice_mesh.GetMesh()->EndNode();
                ++nodeit)
    {
        nodeit->RealArray(Force_vector_tag)[0] = 0.0;
        nodeit->RealArray(Force_vector_tag)[1] = 0.0;
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(Force_vector_tag, NODE, 0);
    BARRIER

    // calculate new Force_vector_tag
    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
                trianit != ice_mesh.GetMesh()->EndCell();
                ++trianit)
    {
        INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
        
        double sigma1_tr = trianit->RealArray(sigma_tag)[0];
        double sigma2_tr = trianit->RealArray(sigma_tag)[1];
        double sigma12_tr = trianit->RealArray(sigma_tag)[2];
        double sigma11_tr = 0.5*(sigma1_tr + sigma2_tr);
        double sigma22_tr = 0.5*(sigma1_tr - sigma2_tr);

        for (int i = 0; i < 3; ++i)
        {
            if ((local_nodes[i]->GetStatus() != Element::Ghost) and 
                (local_nodes[i].Integer(ice_mesh.GetData().IsNodeBnd) == 0))
            { 
                double q0 = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[i*4 + 0];
                double q1 = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[i*4 + 1];
                double q2 = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[i*4 + 2];
                double q3 = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[i*4 + 3];

                double sigma11_n = q0*q0*sigma11_tr + 2.0*q0*q2*sigma12_tr + q2*q2*sigma22_tr;
                double sigma22_n = q1*q1*sigma11_tr + 2.0*q1*q3*sigma12_tr + q3*q3*sigma22_tr;
                double sigma12_n = q0*q1*sigma11_tr + (q1*q2 + q0*q3)*sigma12_tr + q2*q3*sigma22_tr;

                double sigma1_n = sigma11_n + sigma22_n;
                double sigma2_n = sigma11_n - sigma22_n; 

                local_nodes[i]->RealArray(Force_vector_tag)[0] += 
                -0.5*(trianit->RealArray(Nx_matrix_entries_tag)[i])*sigma1_n
                -0.5*(trianit->RealArray(Nx_matrix_entries_tag)[i])*sigma2_n
                -(trianit->RealArray(Ny_matrix_entries_tag)[i])*sigma12_n;
                local_nodes[i]->RealArray(Force_vector_tag)[1] += 
                -0.5*(trianit->RealArray(Ny_matrix_entries_tag)[i])*sigma1_n
                +0.5*(trianit->RealArray(Ny_matrix_entries_tag)[i])*sigma2_n
                -(trianit->RealArray(Nx_matrix_entries_tag)[i])*sigma12_n;
            }
        }

    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(Force_vector_tag, NODE, 0);
    BARRIER
    for (Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                nodeit != ice_mesh.GetMesh()->EndNode();
                ++nodeit)
    {
        if ((nodeit->GetStatus() != Element::Ghost) and
            (nodeit->Integer(ice_mesh.GetData().IsNodeBnd) == 0))
        {
            nodeit->RealArray(Force_vector_tag)[0] /= nodeit->Real(lumped_mass_matrix_entry_tag);
            nodeit->RealArray(Force_vector_tag)[1] /= nodeit->Real(lumped_mass_matrix_entry_tag); 
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(Force_vector_tag, NODE, 0);
    BARRIER
};

void MomentumSolver::UpdateTemporalVelocity()
{
    double rho_air = model_params.GetAirDensity();
    double rho_water = model_params.GetWaterDensity();
    double C_air = model_params.GetAirDragCoeff();
    double C_water = model_params.GetWaterDragCoeff();
    double del_t = model_params.GetTimeStepHours()*3600.0;
    double f_cor = model_params.GetCoriolisParam();
    double beta = momentum_params.GetSolverParams()["beta"].Double;
    double g = model_params.GetGravity();

    if (!momentum_params.GetIsCoriolis())
    {
        f_cor = 0.0;
    }

    if (!momentum_params.GetIsWaterDrag())
    {
        C_water = 0.0;
    }

    if (!momentum_params.GetIsAirDrag())
    {
        C_air = 0.0;
    }

    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                nodeit != ice_mesh.GetMesh()->EndNode();
                ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            double u = nodeit->RealArray(u_old_tag)[0];
            double v = nodeit->RealArray(u_old_tag)[1];
            double u_n = nodeit->RealArray(u_n_tag)[0];
            double v_n = nodeit->RealArray(u_n_tag)[1];
            double u_air = nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_air])[0];
            double v_air = nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_air])[1];
            double u_water = nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_water])[0];
            double v_water = nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_water])[1];
            double m = nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::m]);
            double a = nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]);
            double f_1 = nodeit->RealArray(Force_vector_tag)[0];
            double f_2 = nodeit->RealArray(Force_vector_tag)[1];
            double lev_1 = nodeit->RealArray(Level_vector_tag)[0]; 
            double lev_2 = nodeit->RealArray(Level_vector_tag)[1];
            double abs_u_air = sqrt(u_air*u_air + v_air*v_air);
            double abs_u_min_u_water = sqrt((u-u_water)*(u-u_water) + (v-v_water)*(v-v_water));


            double K = (beta + 1.0) + (del_t/m)*C_water*rho_water*a*abs_u_min_u_water;
            double F = del_t*f_cor;


            double rhs1 = u_n + beta*u + (del_t/m)*(f_1 + C_air*rho_air*a*u_air*abs_u_air + C_water*rho_water*a*u_water*abs_u_min_u_water - m*g*lev_1);
            double rhs2 = v_n + beta*v + (del_t/m)*(f_2 + C_air*rho_air*a*v_air*abs_u_air + C_water*rho_water*a*v_water*abs_u_min_u_water - m*g*lev_2);

            nodeit->RealArray(u_new_tag)[0] = (rhs2*F + rhs1*K)/(K*K + F*F);
            nodeit->RealArray(u_new_tag)[1] = (rhs2*K - rhs1*F)/(K*K + F*F);

            nodeit->RealArray(u_diff_tag)[0] = nodeit->RealArray(u_new_tag)[0] - u;
            nodeit->RealArray(u_diff_tag)[1] = nodeit->RealArray(u_new_tag)[1] - v;

            nodeit->RealArray(u_old_tag)[0] = nodeit->RealArray(u_new_tag)[0];
            nodeit->RealArray(u_old_tag)[1] = nodeit->RealArray(u_new_tag)[1];
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(u_old_tag, NODE, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(u_new_tag, NODE, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(u_diff_tag, NODE, 0);
    BARRIER
};

void MomentumSolver::UpdateSigma()
{
    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
                trianit != ice_mesh.GetMesh()->EndCell();
                ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            trianit->Real(ice_mesh.GetData().TriangleData[ModelVariableNotation::sig1]) = trianit->RealArray(sigma_tag)[0];
            trianit->Real(ice_mesh.GetData().TriangleData[ModelVariableNotation::sig2]) = trianit->RealArray(sigma_tag)[1];
            trianit->Real(ice_mesh.GetData().TriangleData[ModelVariableNotation::sig12]) = trianit->RealArray(sigma_tag)[2];
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().TriangleData[ModelVariableNotation::sig1], CELL, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().TriangleData[ModelVariableNotation::sig2], CELL, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().TriangleData[ModelVariableNotation::sig12], CELL, 0);
    BARRIER
};

void MomentumSolver::UpdateVelocity()
{
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                nodeit != ice_mesh.GetMesh()->EndNode();
                ++nodeit)
    {
        if ((nodeit->GetStatus() != Element::Ghost) and
            (nodeit->Integer(ice_mesh.GetData().IsNodeBnd) == 0))
        {
            nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_ice])[0] = nodeit->RealArray(u_old_tag)[0];
            nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_ice])[1] = nodeit->RealArray(u_old_tag)[1];
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::u_ice], NODE, 0);
    BARRIER
}

void MomentumSolver::UpdateShearDeformation()
{
    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
                trianit != ice_mesh.GetMesh()->EndCell();
                ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            //double q0 = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[0];
            //double q1 = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[1];
            //double q2 = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[2];
            //double q3 = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[3];

            double eps1 = trianit->RealArray(varepsilon_tag)[0];
            double eps2 = trianit->RealArray(varepsilon_tag)[1];
            double eps12 = trianit->RealArray(varepsilon_tag)[2];

            double eps11 = 0.5*(eps1 + eps2);
            double eps22 = 0.5*(eps1 - eps2);

            //double eps11_cart = q0*q0*eps11 + 2.0*q0*q2*eps12 + q2*q2*eps22;
            //double eps22_cart = q1*q1*eps11 + 2.0*q1*q3*eps12 + q3*q3*eps22;
            //double eps12_cart = q0*q1*eps11 + (q1*q2 + q0*q3)*eps12 + q2*q3*eps22;

            trianit->Real(shear_deformation_tag) = sqrt((eps11 - eps22)*(eps11 - eps22) + 4.0*eps12*eps12);
        }
    }
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(shear_deformation_tag, CELL, 0);
    BARRIER
}

 void MomentumSolver::LogMevpError(int stepn)
 {
     // find max abs diff u
    double u_res = std::numeric_limits<double>::min();
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                nodeit != ice_mesh.GetMesh()->EndNode();
                ++nodeit)
    {
        if ((nodeit->GetStatus() != Element::Ghost) and
            (nodeit->Integer(ice_mesh.GetData().IsNodeBnd) == 0))
        {
            double u_val =  sqrt((nodeit->RealArray(u_diff_tag)[0]*nodeit->RealArray(u_diff_tag)[0]) +
                                 (nodeit->RealArray(u_diff_tag)[1]*nodeit->RealArray(u_diff_tag)[1]));
            if (u_val > u_res)
            {
                u_res = u_val;
            }
        }
    }
    BARRIER

    double diff_u_max = 0.0;
    BARRIER
#if defined(USE_MPI)
   MPI_Allreduce(&u_res, &diff_u_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    // find max abs diff sigma
    double res = std::numeric_limits<double>::min();
    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
        trianit != ice_mesh.GetMesh()->EndCell();
        ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            double val = sqrt(trianit->RealArray(sigma_diff_tag)[0]*trianit->RealArray(sigma_diff_tag)[0] +
                              trianit->RealArray(sigma_diff_tag)[1]*trianit->RealArray(sigma_diff_tag)[1] +
                              trianit->RealArray(sigma_diff_tag)[2]*trianit->RealArray(sigma_diff_tag)[2]);
            if (val > res)
            {
                res = val;
            }
        }
    }
    // init_m max
    double diff_sig_max = 0.0;
    BARRIER
#if defined(USE_MPI)
   MPI_Allreduce(&res, &diff_sig_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    BARRIER

    if (ice_mesh.GetMesh()->GetProcessorRank() == 0)
    {
        cout << "iter " << stepn << ": sig diff = " << diff_sig_max << ", u diff = " << diff_u_max << ";" << endl;     
    }

    BARRIER
};


void MomentumSolver::mEVP_evaluation()
{
    UpdateH();
    AssignSigma();
    AssignVelocity();
    CalculateP_0();
    AssembleLevelVector();

    if (momentum_params.GetSolverParams()["number of iterations"].Int != 0)
    {
        for (int pseudostep = 0;
             pseudostep < momentum_params.GetSolverParams()["number of iterations"].Int;
             pseudostep++)
        {
            CalculateVarepsilonDelta();
            UpdateTemporalSigma();
            AssembleForceVector();
            UpdateTemporalVelocity();
            BARRIER
            if (is_verbose_momentum)
            {
                if ((pseudostep%50) == 0)
                {
                    LogMevpError(pseudostep);
                    //ice_mesh.PrintPVTU("PreudoStep"+ to_string(pseudostep) + ".pvtu");
                }
            }
        }
    }
    else
    {
        INMOST_ICE_ERR("currently available only fixed nit mEVP");
    }
    BARRIER
    UpdateSigma();
    UpdateVelocity();
    UpdateShearDeformation();
    BARRIER
};