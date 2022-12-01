#pragma once
#include <initializer_list>
#include <map>

// all physical variables
namespace SIMUG::forc
{
    // types of forcing
    enum forc
    {
        ice,
        air,
        water
    };

    // names of forcing type
    static std::map<forc, std::string> forcname =
    {
        {ice,    "Ice"   },
        {air,    "Air"   },
        {water,  "Water" }
    };

    // types of forcing setup 
    enum setup
    {
        analyt,  // internal analytical setup
        file,    // external file setup 
        online   // external online setup
    };

    // names of forcing setup
    static std::map<setup, std::string> forcsetupname =
    {
        {analyt, "analytical" },
        {file,   "file"       },
        {online, "online"     }
    };

    // notations of forcing setup
    static std::map<std::string, setup> forcsetupnotation =
    {
        {"analytical", analyt },
        {"file",       file   },
        {"online",     online }
    };



}