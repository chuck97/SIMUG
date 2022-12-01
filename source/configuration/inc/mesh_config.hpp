#pragma once
#include "config.hpp"
#include "meshvar.hpp"
#include <nlohmann/json.hpp>
#include <map>
#include <iostream>
#include <fstream>

namespace SIMUG
{
    // model params
    template <typename RT>
    class MeshConfig: public Configuration<RT>
    {
    public:
        // constructors
        MeshConfig();

        MeshConfig(const std::string& json_path_);

        // manually set surface type
        void SetSurf(mesh::surf surf_);

        // get surface type
        mesh::surf GetSurf() const; 

        // manually set variable 
        void SetVar(mesh::var var_notation_,
                    const Var<std::string>& var_);
        
        // manually set map of variables
        void SetVars(const std::map<mesh::var, Var<std::string>>& vars_);

        // get variable
        Var<std::string> GetVar(mesh::var var_notation_) const;

        // get current configuration type
        conf::type GetConfigType() const override;
        
        // log model variables
        void Log(std::ostream& os) const override;

    private:

        // default surface
        mesh::surf surface = mesh::surf::plane;

        // default model variables
        std::map<mesh::var, Var<std::string>> MeshVars = 
        {
            {mesh::var::mi,      {"ice mass",               "kg m-2" }},
            {mesh::var::hi,      {"ice thickness",          "m"      }},
            {mesh::var::ai,      {"ice concentration",      ""       }},
            {mesh::var::ui,      {"ice velocity",           "m s-1"  }},
            {mesh::var::sig,     {"ice stress tensor",      "N m-2"  }},
            {mesh::var::eps,     {"ice strain rate tensor", "s-1"    }},
            {mesh::var::del,     {"delta function",         "s-1"    }},
            {mesh::var::P0,      {"ice pressure",           "N m-1"  }},
            {mesh::var::ua,      {"air velocity",           "m s-1"  }},
            {mesh::var::uw,      {"water velocity",         "m s-1"  }},
            {mesh::var::hw,      {"water level",            "m"      }}
        };
    };
}