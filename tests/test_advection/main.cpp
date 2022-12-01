#include "simug.hpp"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <memory>
#include <vector>
#include <string>
#include <iostream>

#define PLANE_PATH "../../../../SIMUG_v0/MESHES/pmf/Box_low_res.pmf"

using namespace SIMUG;
using SIMUG::operator*;

bool test_advection_assembling_basal()
{
    SIMUG::Timer timer;
    // create local nodes coordinates
    std::vector<std::vector<double>> local_node_coords = 
    {
        {5.0, 6.0},
        {10.0, 2.0},
        {7.0, 8.0}
    };

    std::vector<std::vector<double>> Jac =
    {
        {(local_node_coords[2][0] - local_node_coords[0][0]), (local_node_coords[1][0] - local_node_coords[0][0])},
        {(local_node_coords[2][1] - local_node_coords[0][1]), (local_node_coords[1][1] - local_node_coords[0][1])}
    };

    std::vector<std::vector<double>> Jac_inv = inv(Jac);

    double Jacobian = std::abs(det(Jac));

    // create vector coordinates
    std::vector<std::vector<double>> local_vel_components = 
    {
        {0.05, 0.06},
        {-0.07, -0.08},
        {0.09, -0.1}
    };

    std::vector<double> u_comp = 
    {
        local_vel_components[0][0],
        local_vel_components[1][0],
        local_vel_components[2][0]
    };

    std::vector<double> v_comp = 
    {
        local_vel_components[0][1],
        local_vel_components[1][1],
        local_vel_components[2][1]
    };

    // create local scalar vector
    std::vector<double> local_scal = 
    {
        0.005, 
        -0.006,
        0.007
    };

    double time_step = 5.0*3600.0;

    timer.Launch();
    std::vector<double> local_rhs = LocalTG2RhsAssembling(local_node_coords, local_vel_components, local_scal, time_step);
    timer.Stop();

    double duration = timer.GetMaxTime();
    timer.Reset();

    std::cout << local_rhs[0] << " " << local_rhs[1] << " " << local_rhs[2] << " (" << duration << " ms)" << std::endl;

    auto M_matr = LocaReferenceMassMatrixAssembling();

    auto D_matr = LocaReferenceFirstDerivativesMatrixAssembling();
    auto Dx = D_matr.first;
    auto Dy = D_matr.second;

    auto DD_matr = LocaReferenceSecondDerivativesMatrixAssembling();

    std::vector<double> v1d = {0.0, 0.0, 0.0};
    
    std::vector<std::vector<double>> first_matr = {v1d, v1d, v1d};
    std::vector<std::vector<double>> second_matr = {v1d, v1d, v1d};
    std::vector<std::vector<double>> third_matr = {v1d, v1d, v1d};

    timer.Launch();

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            first_matr[i][j] += M_matr[i][j]*Jacobian;
        }
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                second_matr[i][j] += Dx[i][j][k]*(u_comp[k]*Jac_inv[0][0] + v_comp[k]*Jac_inv[0][1])*Jacobian +
                                     Dy[i][j][k]*(u_comp[k]*Jac_inv[1][0] + v_comp[k]*Jac_inv[1][1])*Jacobian;
            }
        }
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                for (int l = 0; l < 3; ++l)
                {
                    third_matr[i][j] += DD_matr[0][i][j][k][l]*(u_comp[k]*u_comp[l]*Jac_inv[0][0]*Jac_inv[0][0] + u_comp[k]*v_comp[l]*Jac_inv[0][0]*Jac_inv[0][1] + v_comp[k]*u_comp[l]*Jac_inv[0][1]*Jac_inv[0][0] + v_comp[k]*v_comp[l]*Jac_inv[0][1]*Jac_inv[0][1])*Jacobian + 
                                        DD_matr[1][i][j][k][l]*(u_comp[k]*u_comp[l]*Jac_inv[0][0]*Jac_inv[1][0] + u_comp[k]*v_comp[l]*Jac_inv[0][0]*Jac_inv[1][1] + v_comp[k]*u_comp[l]*Jac_inv[0][1]*Jac_inv[1][0] + v_comp[k]*v_comp[l]*Jac_inv[0][1]*Jac_inv[1][1])*Jacobian + 
                                        DD_matr[2][i][j][k][l]*(u_comp[k]*u_comp[l]*Jac_inv[1][0]*Jac_inv[0][0] + u_comp[k]*v_comp[l]*Jac_inv[1][0]*Jac_inv[0][1] + v_comp[k]*u_comp[l]*Jac_inv[1][1]*Jac_inv[0][0] + v_comp[k]*v_comp[l]*Jac_inv[1][1]*Jac_inv[0][1])*Jacobian +
                                        DD_matr[3][i][j][k][l]*(u_comp[k]*u_comp[l]*Jac_inv[1][0]*Jac_inv[1][0] + u_comp[k]*v_comp[l]*Jac_inv[1][0]*Jac_inv[1][1] + v_comp[k]*u_comp[l]*Jac_inv[1][1]*Jac_inv[1][0] + v_comp[k]*v_comp[l]*Jac_inv[1][1]*Jac_inv[1][1])*Jacobian; 
                                       
                }
            }
        }
    } 



    std::vector<double> res_rhs = (first_matr + time_step*second_matr - (time_step*time_step*0.5)*third_matr)*local_scal;
    timer.Stop();
    duration = timer.GetMaxTime();
    timer.Reset();
    std::cout << res_rhs[0] << " " << res_rhs[1] << " " << res_rhs[2] << " (" << duration << " ms)" << std::endl;

    return true;
}


int main()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (test_advection_assembling_basal())
    {
        if (rank == 0)
            std::cout << "Adv assembling test: OK!\n";
    }
    else
        SIMUG_ERR("Adv assembling test: FAILED!\n");

    BARRIER

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}