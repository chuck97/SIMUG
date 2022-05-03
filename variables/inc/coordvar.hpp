#pragma once
#include <initializer_list>
#include <string>
#include <map>

// all enum variables for coordinates
namespace SIMUG::coords
{
    // coordinate names
    enum coord
    {
        model, // model coordinates
        geo,   // geographical coordinates
        cart,  // cartesian coordinates
        temp1, // first temporal coordinates
        temp2, // second temporal coordinates
        temp3  // third temporal coordinates 
    };

    constexpr static std::initializer_list<coord> coords = 
    {
        model, geo, cart, temp1, temp2, temp3
    };


    // coord type -> name
    static std::map<coord, std::string> coordname =
    {
        {model,     "model"        },
        {geo,       "geographical" },
        {cart,      "cartesian"    },
        {temp1,     "temporal1"    },
        {temp2,     "temporal2"    },
        {temp3,     "temporal3"    },
    };

    // coord name -> type
    static std::map<std::string, coord> coordnotation =
    {
        {"model",        model },
        {"geographical", geo   },
        {"cartesian",    cart  },
        {"temporal 1",   temp1 },
        {"temporal 2",   temp2 },
        {"temporal 3",   temp3 }
    };
};  