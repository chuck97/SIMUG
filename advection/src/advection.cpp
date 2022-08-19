#include "advection.hpp"

using namespace std;
using namespace INMOST;
using namespace SIMUG;

AdvectionSolver::AdvectionSolver(SIMUG::IceMesh* mesh_,
                                 double time_step_,
                                 velocity_tag vel_tag_,
                                 SIMUG::adv::timeScheme adv_time_scheme_,
                                 SIMUG::adv::spaceScheme adv_space_scheme_,
                                 SIMUG::adv::advFilter adv_filter_):
        mesh(mesh_),
        time_step(time_step_),
        vel_tag(vel_tag_),
        adv_time_scheme(adv_time_scheme_),
        adv_space_scheme(adv_space_scheme_),
        adv_filter(adv_filter_)
{
    // initialize timer and logger
    SIMUG::Logger adv_log(std::cout);
    SIMUG::Timer adv_timer;

    // log constructor
    if (mesh->GetMesh()->GetProcessorRank()==0)
    {
        adv_log.Log("================== Advection solver initialization ==================\n");
        adv_log.Log("Advection time scheme: " + adv::advTimeSchemeName.at(adv_time_scheme) + "\n");
        adv_log.Log("Advection space scheme: " + adv::advSpaceSchemeName.at(adv_space_scheme) + "\n");
        adv_log.Log("Advection filter: " + adv::advFilterName.at(adv_filter) + "\n");
    }
    BARRIER
}


void AdvectionSolver::TransportScalars()
{
    SIMUG::Logger adv_log(std::cout);
    SIMUG::Timer adv_timer;
    double duration;

    adv_timer.Launch();
    for (auto& scal_tag: scal_tags)
    {
        Evaluate(vel_tag, scal_tag);
    } 
    adv_timer.Stop();
    duration = adv_timer.GetMaxTime();
    adv_timer.Reset();

    if (mesh->GetMesh()->GetProcessorRank() == 0)
        adv_log.Log("\nSuccessfull transport of " + std::to_string(scal_tags.size()) + " scalars! (" + std::to_string(duration) + " ms)\n");
    
    PrintProfiling();
    BARRIER
};