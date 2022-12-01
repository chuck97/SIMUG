#pragma once
#include <initializer_list>
#include <map>

// all physical variables
namespace SIMUG::phys
{
    //  variables 
    enum physVar
    {
        g,      // gravity acceleration
        rhow,   // water density
        rhoa,   // air density
        rhoi,   // ice density
        Cw,     // water-ice drag coefficient
        Ca,     // air-ice drag coefficient
        C,      // ice concentration parameter
        pstr,   // p* coefficient
        e,      // elipse eccentricity
        delmin, // minimal value of delta
        f,      // Coriolis parameter
        Reth    // Earth radius
    };

    // model variable list
    constexpr static std::initializer_list<physVar> vars = 
    {
        g, rhow, rhoa, rhoi, Cw, Ca, C, pstr, e, delmin, f, Reth 
    };


    // variable notation -> variable name
    static std::map<physVar, std::string> name =
    {
        {g,      "gravity acceleration"       },
        {rhow,   "water density"              },
        {rhoa,   "air density"                },
        {rhoi,   "ice density"                },
        {Cw,     "water-ice drag coefficient" },
        {Ca,     "air-ice drag coefficient"   },
        {C,      "concentration coefficient"  },
        {pstr,   "pressure star coefficient"  },
        {e,      "ellipse aspect ratio"       },
        {delmin, "minimal delta value"        },
        {f,      "Coriolis parameter"         },
        {Reth,   "Earth radius"               }
    };
}