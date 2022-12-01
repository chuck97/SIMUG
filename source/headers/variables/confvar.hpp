#pragma once
#include <initializer_list>
#include <string>
#include <map>

// all configuration variables 
namespace SIMUG::conf
{
    // configuration class types
    enum type
    {
        notype, // no-type config class
        phys,   // physical constants config class
        model,  // model parameters config class
        adv,    // advection parameters config class
        dyn,    // momentum solver config class
        mesh,   // mesh config params
        forc    // forcing config params
    };

    // list of all conf types
    constexpr static std::initializer_list<type> types = 
    {
        notype, phys, model, adv, dyn, mesh, forc
    };

    // configuration class type -> name
    static std::map<type, std::string> name =
    {
        {notype, "No type configuration class"            },
        {phys,   "Physical constants configuration class" },
        {model,  "Model parameters configuration class"   },
        {adv,    "Advection scheme configuration class"   },
        {dyn,    "Momentum scheme configuration class"    },
        {mesh,   "Mesh info configuration class"          },
        {forc,   "Forcing configuration class"            }
    };

    // configuration class type -> json key
    static std::map<type, std::string> jsonkey =
    {
        {phys,   "Physical"      },
        {model,  "Model"         },
        {adv,    "Advection"     },
        {dyn,    "Momentum"      },
        {mesh,   "Mesh"          },
        {forc,   "Forcing"       }
    };
};  