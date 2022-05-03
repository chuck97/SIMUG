#pragma once
#include "config.hpp"
#include "forcvar.hpp"
#include "meshvar.hpp"
#include <nlohmann/json.hpp>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>

namespace SIMUG
{
    // model params
    template <typename RT, forc::forc ftype>
    class ForcConfig: public Configuration<RT>
    {
    public:
        // constructors
        ForcConfig();

        ForcConfig(const std::string& json_path_);

        // manually set forcing variable 
        void SetVar(mesh::var mesh_var_notation_,
                   const std::pair<Var<std::string>, forc::setup>& mesh_var_);

        // get forcing variable 
        std::pair<Var<std::string>, forc::setup> 
        GetVar(mesh::var mesh_var_notation_);

        // get current configuration type
        conf::type GetConfigType() const override;

        // get current forcing type
        forc::forc GetForcingType() const;
        
        // log model variables
        void Log(std::ostream& os) const override;

    private:

        // data
        forc::forc forctype;
        std::map<mesh::var, std::pair<Var<std::string>, forc::setup>> ForcVars;
    };
}