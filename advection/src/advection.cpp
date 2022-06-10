#include "advection.hpp"

using namespace std;
using namespace INMOST;
using namespace SIMUG;

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