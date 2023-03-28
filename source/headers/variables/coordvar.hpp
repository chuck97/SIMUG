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
        cart,  // cartesian coordinates
        topaz  // TOPAZ coordinates 
    };

    constexpr static std::initializer_list<coordType> coords = 
    {
        model, geo, cart, topaz
    };


    // coord type -> name
    const static std::map<coordType, std::string> coordname =
    {
        {model,     "model"        },
        {geo,       "geographical" },
        {cart,      "Cartesian"    },
        {topaz,     "TOPAZ"}
    };

    // coord name -> type
    const static std::map<std::string, coordType> coordnotation =
    {
        {"model",        model },
        {"geographical", geo   },
        {"Cartesian",    cart  },
        {"TOPAZ",        topaz }
    };
};  