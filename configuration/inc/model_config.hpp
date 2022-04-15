#pragma once
#include "config.hpp"
#include "modvar.hpp"
#include <nlohmann/json.hpp>
#include <map>
#include <iostream>
#include <fstream>

namespace SIMUG
{
    // model params
    template <typename RT>
    class ModelConfig: public Configuration<RT>
    {
    public:
        // constructors
        ModelConfig();

        ModelConfig(const Var<RT>& time_step_,
                    const Var<RT>& total_time_,
                    const Var<RT>& mesh_);
        
        ModelConfig(const Var<RT>& time_step_,
                    const Var<RT>& total_time_,
                    const Var<RT>& amin_,
                    const Var<RT>& hmin_,
                    const Var<RT>& mesh_);

        ModelConfig(const std::string& json_path_);

        // manually set variable 
        void SetVar(model::var var_notation_,
                    const Var<RT>& var_);
        
        // manually set map of variables
        void SetVars(const std::map<model::var, Var<RT>>& vars_);

        // get variable
        Var<RT> GetVar(model::var var_notation_) const;

        // get current configuration type
        conf::type GetConfigType() const override;
        
        // log model variables
        void Log(std::ostream& os) const override;

    private:

        // default model variables
        std::map<model::var, Var<RT>> ModelVars = 
        {
            {model::var::dt,   {60.0  , "s"}},
            {model::var::T,    {3600.0, "s"}},
            {model::var::amin, {1e-6  , "" }},
            {model::var::hmin, {1e-2  , "m"}},
            {model::var::mesh, {0.0   , "" }}
        };
    };
}