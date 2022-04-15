#pragma once
#include <initializer_list>
#include <map>

// all mesh data variables
namespace SIMUG::mesh
{
    enum meshVar
    {
        mi,      // ice mass
        hi,      // ice thickness
        ai,      // ice concentration
        ui,      // ice velocity
        sig,     // ice stress tensor
        eps,     // ice strain rate tensor
        del,     // delta function of strain rate invariants
        P0,      // ice pressure
        ua,      // air velocity
        uw,      // water velocity
        hw       // water level  
    }; 

    enum meshDim
    {
        scalar,
        vector,
        tensor
    }; 

    // mesh variable notation -> mesh variable name
    static std::map<meshVar, std::string> meshVarName =
    {
        {mi,  "ice mass"               },
        {hi,  "ice thickness"          },
        {ai,  "ice concentration"      },
        {ui,  "ice velocity"           },
        {sig, "ice stress tensor"      },
        {eps, "ice strain rate tensor" },
        {del, "delta function"         },
        {P0,  "ice pressure"           },
        {ua,  "air velocity"           },
        {uw,  "water velocity"         },
        {hw,  "water level"            }
    };

    // mesh variable name -> mesh variable name
    static std::map<std::string, meshVar> meshVarNotation =
    {
        { "ice mass"               , mi  },
        { "ice thickness"          , hi  },
        { "ice concentration"      , ai  },
        { "ice velocity"           , ui  },
        { "ice stress tensor"      , sig },
        { "ice strain rate tensor" , eps },
        { "delta function"         , del },
        { "ice pressure"           , P0  },
        { "air velocity"           , ua  },
        { "water velocity"         , uw  },
        { "water level"            , hw  }
    };

    // variable dimension notation -> variable dimension name
    static std::map<meshDim, std::string> meshDimName =
    {
        {scalar,  "scalar" },
        {vector,  "vector" },
        {tensor,  "tensor" }
    };

    // variable dimension name -> variable dimension notation
    static std::map<std::string, meshDim> meshDimNotation =
    {
        {"scalar", scalar },
        {"vector", vector },
        {"tensor", tensor }
    };
} 