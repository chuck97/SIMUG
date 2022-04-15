#pragma once

// all model variables
namespace SIMUG::model
{
    // model variables
    enum var
    {
        dt,     // time step
        T,      // total time
        amin,   // minimal ice concentration  
        hmin,   // minimal ice thickness
        mesh    // mesh filename
    };

    // model variables list
    constexpr static std::initializer_list<var> vars = 
    {
        dt, T, amin, hmin, mesh
    };

    // variable notation -> variable name
    static std::map<var, std::string> name =
    {
        {dt,     "time step"                  },
        {T,      "total time"                 },
        {amin,   "minimal ice concentration"  },
        {hmin,   "minimal ice thickness"      },
        {mesh,   "mesh file"                  }
    };
}

