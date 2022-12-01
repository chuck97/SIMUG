#pragma once
#include "config.hpp"
#include "physvar.hpp"
#include <nlohmann/json.hpp>
#include <map>
#include <iostream>
#include <fstream>

namespace SIMUG
{
    // model params
    template <typename RT>
    class PhysConfig: public Configuration<RT>
    {
    public:
        // constructors
        PhysConfig();

        PhysConfig(const std::string& json_path_);

        // manually set variable 
        void SetVar(phys::var var_notation_,
                    const Var<RT>& var_);
        
        // manually set map of variables
        void SetVars(const std::map<phys::var, Var<RT>>& vars_);

        // get variable
        Var<RT> GetVar(phys::var var_notation_) const;

        // get current configuration type
        conf::type GetConfigType() const override;
        
        // log model variables
        void Log(std::ostream& os) const override;

    private:

        // default model variables
        std::map<phys::var, Var<RT>> PhysVars = 
        {
            {phys::var::g,      {9.8,     "m s-2"  }},
            {phys::var::rhow,   {1026.0,  "kg m-3" }},
            {phys::var::rhoa,   {1.3,     "kg m-3" }},
            {phys::var::rhoi,   {900.0,   "kg m-3" }},
            {phys::var::Cw,     {55e-4,   ""       }},
            {phys::var::Ca,     {12e-4,   ""       }},
            {phys::var::C,      {20.0,    ""       }},
            {phys::var::pstr,   {275e2,   "Pa"     }},
            {phys::var::e,      {2.0,     ""       }},
            {phys::var::delmin, {2e-9,    "s-1"    }},
            {phys::var::f,      {1.46e-4, "s-1"    }},
            {phys::var::Reth,   {6371e3,  "m"      }}
        };
    };
}