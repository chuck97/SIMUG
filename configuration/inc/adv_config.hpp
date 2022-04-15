#pragma once
#include "config.hpp"
#include "advvar.hpp"
#include <nlohmann/json.hpp>
#include <map>
#include <iostream>
#include <fstream>

namespace SIMUG
{
    // model params
    template <typename RT>
    class AdvConfig: public Configuration<RT>
    {
    public:
        // constructors
        AdvConfig();

        AdvConfig(const std::string& json_path_);

        // manually set variable 
        void SetScheme(adv::scheme scheme_notation_);
        
        // manually set map of variables
        void SetFilter(adv::filter filter_notation_);

        // manually set real advection parameter
        void SetRealParam(const std::string& param_name_,
                                         RT  param_value_);

        // manually set int advection parameter
        void SetIntParam(const std::string& param_name_,
                                       int  param_value_);

        // get advection scheme
        adv::scheme GetScheme() const;

        // get advection filter
        adv::filter GetFilter() const;

        // get real advection parameter
        RT GetRealParam(const std::string param_name_);

        // get int advection parameter
        int GetIntParam(const std::string param_name_);

        // get current configuration type
        conf::type GetConfigType() const override;
        
        // log model variables
        void Log(std::ostream& os) const override;

    private:

        // advection scheme and filter
        adv::scheme scheme = adv::scheme::TG2;
        adv::filter filter = adv::filter::Zalesak;
        std::map<std::string, RT> realparameters = {{"Zalesak diffusion parameter", 0.5}};
        std::map<std::string, RT> intparameters;
    };
}