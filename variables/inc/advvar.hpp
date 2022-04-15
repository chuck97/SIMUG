#pragma once
#include <initializer_list>
#include <string>
#include <map>

// all enum variables for advection time discretezatiion
namespace SIMUG::adv
{
    // advection time discretization schemes
    enum scheme
    {
        TG2,
        TTG2,
        TTG3,
        TTG4
    };

    constexpr static std::initializer_list<scheme> schemes = 
    {
        TG2, TTG2, TTG3, TTG4
    };

    // advection scheme filter type
    enum filter
    {
        none,
        Zalesak
    };

    constexpr static std::initializer_list<filter> filters = 
    {
        none, Zalesak
    };

    // is advection scheme single step
    static std::map<scheme, bool> is_fct =
    {
        {TG2,   true  },
        {TTG2,  false },
        {TTG3,  false },
        {TTG4,  false }
    };

    // advection scheme type -> name
    static std::map<scheme, std::string> advschname =
    {
        {TG2,   "Single-step Taylor-Galerkin of 2nd order on time (TG2)"   },
        {TTG2,  "Two-step step Taylor-Galerkin of 2nd order on time (TTG2)"},
        {TTG3,  "Two-step step Taylor-Galerkin of 3rd order on time (TTG3)"},
        {TTG4,  "Two-step step Taylor-Galerkin of 4th order on time (TTG4)"}
    };

    // advection scheme type -> short name
    static std::map<scheme, std::string> advschshortname =
    {
        {TG2,   "TG2" },
        {TTG2,  "TTG2"},
        {TTG3,  "TTG3"},
        {TTG4,  "TTG4"}
    };

    // advection scheme short name -> type
    static std::map<std::string, scheme> advschnotation =
    {
        {"TG2",  TG2 },
        {"TTG2", TTG2},
        {"TTG3", TTG3},
        {"TTG4", TTG4}
    };

    // advection scheme filter -> name
    static std::map<filter, std::string> advfiltname =
    {
        {none,    "none"                    },
        {Zalesak, "Zalesak flux-correction" }
    };

    // advection scheme name -> type
    static std::map<std::string, filter> advfiltnotation =
    {
        {"none",                    none   },
        {"Zalesak flux-correction", Zalesak}
    };

    //
};  