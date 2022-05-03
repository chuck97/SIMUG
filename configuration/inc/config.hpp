#pragma once
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include "defines.hpp"
#include "stringtrim.hpp"
#include "allvars.hpp"

// physical, model and mesh variable classes
namespace SIMUG
{
    template <typename T>
    struct Var
    {
        T value;
        std::string unit;
    };

    // configuration base class
    template <typename RT>
    class Configuration
    {
    public:
        Configuration();
        Configuration(const std::string& json_path_);
        std::string GetJsonFilename() const;
        virtual conf::type GetConfigType() const = 0;
        virtual void Log(std::ostream& os) const = 0;

    protected:
        std::string json_filename;
        conf::type config_type;
    };
}