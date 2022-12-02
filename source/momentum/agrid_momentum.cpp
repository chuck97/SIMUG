#include "momentum.hpp"

namespace SIMUG
{
    AgridMomentumSolver::AgridMomentumSolver(SIMUG::IceMesh* mesh_,
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
                           SIMUG::dyn::spaceScheme::CFE,
                           mom_press_param_,
                           mom_bc_type_),
            real_params(real_params_),
            integer_params(integer_params_)
    {

    }

    Agrid_mEVP_Solver::Agrid_mEVP_Solver(SIMUG::IceMesh* mesh_,
                                         double time_step_,
                                         velocity_tag vel_tag_,
                                         scalar_tag mass_tag_,
                                         scalar_tag conc_tag_,
                                         scalar_tag thick_tag_,
                                         SIMUG::dyn::pressParam mom_press_param_,
                                         SIMUG::dyn::bcType mom_bc_type_,
                                         const std::vector<double>& real_params_,
                                         const std::vector<int>& integer_params_) :
            AgridMomentumSolver(mesh_,
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

    }

    void Agrid_mEVP_Solver::ComputeVelocity()
    {
        std::cout << "calculations" << std::endl;
    }

    void Agrid_mEVP_Solver::PrintProfiling()
    {
        std::cout << "profiling" << std::endl;
    }
}