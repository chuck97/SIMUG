#pragma once
#include <string>

namespace SIMUG
{
    struct NcFileInfo
    {
        std::string filename;
        std::string xname;
        std::string yname;
        std::string lonname;
        std::string latname;
        std::string scale_factor_name;
        std::string invalid_value_name;
        std::string offset_name;
        bool is_depth;
    };
}