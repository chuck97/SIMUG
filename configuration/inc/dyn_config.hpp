#pragma once
#include "config.hpp"
#include "dynvar.hpp"
#include <nlohmann/json.hpp>
#include <map>
#include <iostream>
#include <fstream>

namespace SIMUG
{
    // dynamics params
    template <typename RT>
    class DynConfig: public Configuration<RT>
    {
    public:
        // constructors
        DynConfig();

        DynConfig(const std::string& json_path_);

        // manually set scheme 
        void SetScheme(dyn::scheme scheme_notation_);

        // manually set pressure type 
        void SetPress(dyn::press press_notation_);

        // manually set boundary condition type 
        void SetBC(dyn::bc bc_notation_);
        
        // manually set real dynamics parameter
        void SetRealParam(const std::string& param_name_,
                                          RT param_value_);

        // manually set int dynamics parameter
        void SetIntParam(const std::string& param_name_,
                                        int param_value_);

        // get current scheme
        dyn::scheme GetScheme() const;

        // get current pressure type
        dyn::press GetPress() const;

        // get current bc type
        dyn::bc GetBC() const;

        // get real dynamics parameter
        RT GetRealParam(const std::string& param_name_);

        // get int dynamics parameter
        int GetIntParam(const std::string& param_name_);

        // get current configuration type
        conf::type GetConfigType() const override;
        
        // log model variables
        void Log(std::ostream& os) const override;

    private:

        dyn::scheme scheme = dyn::scheme::mEVP;
        dyn::press pressure = dyn::press::clas;
        dyn::bc boundary_conditions = dyn::bc::noslip;

        // default dynamics parameters
        std::map<std::string, RT> realparameters = 
        {
            {"alpha mEVP", 500.0},
            {"beta mEVP", 500.0}
        };

        std::map<std::string, int> intparameters = 
        {
            {"iters mEVP", 500}
        };
    };
}