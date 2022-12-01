#pragma once
#include<map>

// all model variables
namespace SIMUG::model
{
    // model variables
    enum modelVar
    {
        dt,     // time step
        T,      // total time
        amin,   // minimal ice concentration  
        hmin,   // minimal ice thickness
        mesh    // mesh filename
    };

    // model variable notation -> variable name
    static const std::map<modelVar, std::string> modelVarName =
    {
        {dt,     "time step"                 },
        {T,      "total time"                },
        {amin,   "minimal ice concentration" },
        {hmin,   "minimal ice thickness"     },
        {mesh,   "mesh file"                 }
    };

    // model variable name -> variable notation 
    static const std::map<std::string, modelVar> modelVarNotation =
    {
        {"time step"                , dt   },
        {"total time"               , T    },
        {"minimal ice concentration", amin },
        {"minimal ice thickness"    , hmin },
        {"mesh file"                , mesh }
    };
}

