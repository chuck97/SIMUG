#pragma once
#include <initializer_list>
#include <string>
#include <map>

// all enum variables for advection time discretezatiion
namespace SIMUG::adv
{
    // advection time schemes
    enum timeScheme
    {
        TG2,
        TTG2,
        TTG3,
        TTG4,
        Euler,
        TRK2
    };

    // advection space schemes
    enum spaceScheme
    {
        CFE,
        FVupwind,
        MUST
    };

    // advection filter type
    enum advFilter
    {
        none,
        Zalesak,
        MUSTfilter
    };

    // is advection time scheme single step
    static std::map<timeScheme, bool> is_single_step =
    {
        {TG2,   true  },
        {TTG2,  false },
        {TTG3,  false },
        {TTG4,  false },
        {Euler, true  },
        {TRK2,  false }
    };

    // advection time scheme type -> name
    static std::map<timeScheme, std::string> advTimeSchemeName =
    {
        {TG2,    "One-step Taylor-Galerkin of 2nd order on time (TG2)"      },
        {TTG2,   "Two-step step Taylor-Galerkin of 2nd order on time (TTG2)"},
        {TTG3,   "Two-step step Taylor-Galerkin of 3rd order on time (TTG3)"},
        {TTG4,   "Two-step step Taylor-Galerkin of 4th order on time (TTG4)"},
        {Euler,  "One-step Euler scheme of 1st order on time (Euler)"       },
        {TRK2,   "Two-step Runge-Kutta scheme of 2nd order on time (TRK2)"  }
    };

    // advection space scheme type -> name
    static std::map<spaceScheme, std::string> advSpaceSchemeName =
    {
        {CFE,      "Continous finite element scheme of 2nd order in space (CFE)" },
        {FVupwind, "Finite volume upwind scheme of 1st order in space (FVupwind)"},
        {MUST,     "Finite volume monotonic upwind scheme of 2nd order in space (MUST)"},
    };

    // advection scheme filter -> name
    static std::map<advFilter, std::string> advFilterName =
    {
        {none,       "none"                                     },
        {Zalesak,    "Zalesak flux-correction filter (Zalesak)" },
        {MUSTfilter, "MUST flux filter (MUSTfilter)"            }
    };
};  