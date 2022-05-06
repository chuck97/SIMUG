#pragma once
#include <initializer_list>
#include <string>
#include <map>

// all enum variables for coordinates
namespace SIMUG::coord
{
    // coordinate names
    enum coordType
    {
        model, // model coordinates
        geo,   // geographical coordinates
        cart   // cartesian coordinates
    };

    constexpr static std::initializer_list<coordType> coords = 
    {
        model, geo, cart
    };


    // coord type -> name
    const static std::map<coordType, std::string> coordname =
    {
        {model,     "model"        },
        {geo,       "geographical" },
        {cart,      "Cartesian"    }
    };

    // coord name -> type
    const static std::map<std::string, coordType> coordnotation =
    {
        {"model",        model },
        {"geographical", geo   },
        {"Cartesian",    cart  }
    };
};  