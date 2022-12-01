#pragma once

#include "defines.hpp"
#include "vecmath.hpp"
#include "function.hpp"
#include "numerical.hpp"

namespace SIMUG
{
    std::vector<std::vector<double>> LocalMassMatrixAssembling(const std::vector<std::vector<double>>& adj_nodes);

    double L_ij_entry(const std::vector<std::vector<double>>& node_coords,
                      const ScalarFunction& phi_i,
                      const ScalarFunction& phi_j,
                      const ScalarFunction& u0,
                      const ScalarFunction& u1);


    std::vector<double> LocalTG2RhsAssembling(const std::vector<std::vector<double>>& adj_nodes,
                                                        const std::vector<std::vector<double>>& uvalues,
                                                        const std::vector<double>& localmass,
                                                        double time_step);

    std::vector<double> LocalCG2RhsAssembling(const std::vector<std::vector<double>>& adj_nodes,
                                                        const std::vector<std::vector<double>>& uvalues,
                                                        const std::vector<double>& localmass,
                                                        double time_step);

    std::vector<double> LocalTTG2RhsAssembling(const std::vector<std::vector<double>>& adj_nodes,
                                                        const std::vector<std::vector<double>>& uvalues,
                                                        const std::vector<double>& localmass,
                                                        const std::vector<double>& localmass_half,
                                                        double time_step,
                                                        int step);

    std::vector<double> LocalTTG3RhsAssembling(const std::vector<std::vector<double>>& adj_nodes,
                                                        const std::vector<std::vector<double>>& uvalues,
                                                        const std::vector<double>& localmass,
                                                        const std::vector<double>& localmass_half,
                                                        double time_step,
                                                        int step);

    std::vector<double> LocalTTG4RhsAssembling(const std::vector<std::vector<double>>& adj_nodes,
                                                        const std::vector<std::vector<double>>& uvalues,
                                                        const std::vector<double>& localmass,
                                                        const std::vector<double>& localmass_half,
                                                        double time_step,
                                                        int step);

    using Vec2d = std::vector<std::vector<double>>;
    using Vec3d = std::vector<std::vector<std::vector<double>>>;
    using Vec4d = std::vector<std::vector<std::vector<std::vector<double>>>>;

    Vec2d LocaReferenceMassMatrixAssembling();
    std::pair<Vec3d, Vec3d> LocaReferenceFirstDerivativesMatrixAssembling();
    std::vector<Vec4d> LocaReferenceSecondDerivativesMatrixAssembling();

    Vec2d FastLocalMassMatrixAssembling(const Vec2d& mass_tensor,
                                        const std::vector<double>& Jacobi_info_vec);

    Vec2d FastLocalFirstDerivMatrixAssembling(const std::pair<Vec3d, Vec3d>& first_deriv_tensors,
                                              const std::vector<double>& u_components,
                                              const std::vector<double>& v_components,
                                              const std::vector<double>& Jacobi_info_vec);

    Vec2d FastLocalSecondDerivMatrixAssembling(const std::vector<Vec4d>& second_deriv_tensors,
                                               const std::vector<double>& u_components,
                                               const std::vector<double>& v_components,
                                               const std::vector<double>& Jacobi_info_vec);
    

}