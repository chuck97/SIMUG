#include "momentum.hpp"

namespace SIMUG
{
    MomentumSolver::MomentumSolver(SIMUG::IceMesh* mesh_,
                                   double time_step_,
                                   velocity_tag vel_tag_,
                                   scalar_tag mass_tag_,
                                   scalar_tag conc_tag_,
                                   scalar_tag thick_tag_,
                                   SIMUG::dyn::timeScheme mom_time_scheme_,
                                   SIMUG::dyn::spaceScheme mom_space_scheme_,
                                   SIMUG::dyn::pressParam mom_press_param_,
                                   SIMUG::dyn::bcType mom_bc_type_):
            mesh(mesh_),
            time_step(time_step_),
            vel_tag(vel_tag_),
            mass_tag(mass_tag_),
            conc_tag(conc_tag_),
            thick_tag(thick_tag_),
            mom_time_scheme(mom_time_scheme_),
            mom_space_scheme(mom_space_scheme_),
            mom_press_param(mom_press_param_),
            mom_bc_type(mom_bc_type_)
    {
        // initialize timer and logger
        SIMUG::Logger mom_log(std::cout);
        SIMUG::Timer mom_timer;

        // log constructor
        if (mesh->GetMesh()->GetProcessorRank()==0)
        {
            mom_log.Log("================== Momentum solver initialization ==================\n");
            mom_log.Log("Momentum time scheme: " + dyn::momTimeSchemeName.at(mom_time_scheme) + "\n");
            mom_log.Log("Momentum space scheme: " + dyn::momSpaceSchemeName.at(mom_space_scheme) + "\n");
            mom_log.Log("Momentum pressure parameterization: " + dyn::momPressParamName.at(mom_press_param) + "\n");
            mom_log.Log("Momentum boundary conditions: " + dyn::momBcTypeName.at(mom_bc_type) + "\n");
        }
        BARRIER
    }
}