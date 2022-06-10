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
}