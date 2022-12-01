#pragma once
#include <initializer_list>
#include <string>
#include <map>

// all dynamic variables
namespace SIMUG::dyn
{
    // momentum solver schemes
    enum scheme
    {
        EVP,      // classical EVP solver 
        mEVP,     // modified EVP solver
        aEVP,     // adaptive EVP solver
        mEVPopt,  // optimized mEVP solver
        Picard,   // Picard solver
        Newton    // Jacobian free Newton-Krylov solver
    };

    // list of momentum schemes
    constexpr static std::initializer_list<scheme> schemes = 
    {
        EVP, mEVP, aEVP, mEVPopt, Picard, Newton 
    };

    // pressure parametrization methods
    enum press
    {
        clas, // classical pressure parametrization
        repl  // replaced pressure parametrization
    };

    // list of pressure parametrizations
    constexpr static std::initializer_list<press> press_list = 
    {
        clas, repl 
    };

    // boundary condition types
    enum bc
    {
        noslip,   // no-slip boundary conditions  (u = 0 on BND)
        slip,     // slip boundary conditions     (u|_n = 0 on BND)
        fric      // friction boundary conditions (+ friction in momentum equ. for BND nodes)
    };

    // list of boundary condition types
    constexpr static std::initializer_list<bc> bc_list = 
    {
        noslip, slip, fric
    };

    // momentum solver scheme type -> name
    static std::map<scheme, std::string> momschname =
    {
        {EVP,     "EVP momentum solver"     },
        {mEVP,    "mEVP momentum solver"    },
        {aEVP,    "aEVP momentum solver"    },
        {mEVPopt, "mEVPopt momentum solver" },
        {Picard,  "Picard momentum solver"  },
        {Newton,  "JFNKS momentum solver"   }
    };

    // momentum solver short name -> type
    static std::map<std::string, scheme> momschnotation =
    {
        {"EVP",     EVP     },
        {"mEVP",    mEVP    },
        {"aEVP",    aEVP    },
        {"mEVPopt", mEVPopt },
        {"Picard",  Picard  },
        {"Newton",  Newton  }
    };

    // pressure parametrization -> name
    static std::map<press, std::string> pressname =
    {
        {clas,  "classical"  },
        {repl,  "replacement"}
    };

    // pressure name -> parametrization
    static std::map<std::string, press> pressnotation =
    {
        {"classical",   clas},
        {"replacement", repl}
    };

    // bc type -> name
    static std::map<bc, std::string> bcname =
    {
        {noslip,  "no-slip (u = 0 on BND)"                         },
        {slip,    "slip (u|_n = 0 on BND)"                         },
        {fric,    "friction (+ friction in mom. eq. for BND nodes)"}
    };

    // bc name -> type
    static std::map<std::string, bc> bcnotation =
    {
        {"no-slip" , noslip },
        {"slip"    , slip   },
        {"friction", fric   }
    };
}