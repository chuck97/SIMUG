#pragma once
#include <initializer_list>
#include <string>
#include <map>

// all dynamic variables
namespace SIMUG::dyn
{
    // momentum solver schemes
    enum timeScheme
    {
        EVP,      // classical EVP solver 
        mEVP,     // modified EVP solver
        aEVP,     // adaptive EVP solver
        mEVPopt,  // optimized mEVP solver
        Picard,   // Picard solver
        Newton    // Jacobian free Newton-Krylov solver
    };

    enum spaceScheme
    {
        CFE,      // Continous Finite Elements
        stabCR    // Stabilized nonconforming Crouzeix-Raviart
    };

    enum pressParam
    {
        clas,    // classical definition
        repl     // pressure with replacement
    };

    enum bcType
    {
        noslip, // no-slip boundary conditions
        slip,   // slip boundary conditions
        fric    // friction boundary conditions
    };

    // momentum solver time scheme type -> name
    static std::map<timeScheme, std::string> momTimeSchemeName =
    {
        {EVP,     "EVP momentum solver"     },
        {mEVP,    "mEVP momentum solver"    },
        {aEVP,    "aEVP momentum solver"    },
        {mEVPopt, "mEVPopt momentum solver" },
        {Picard,  "Picard momentum solver"  },
        {Newton,  "JFNKS momentum solver"   }
    };

    // momentum solver time scheme short name -> type
    static std::map<std::string, timeScheme> momTimeSchemeNotation =
    {
        {"EVP",     EVP     },
        {"mEVP",    mEVP    },
        {"aEVP",    aEVP    },
        {"mEVPopt", mEVPopt },
        {"Picard",  Picard  },
        {"Newton",  Newton  }
    };

    // momentum solver space scheme type -> name
    static std::map<spaceScheme, std::string> momSpaceSchemeName =
    {
        {CFE,     "continous finite elements"                                 },
        {stabCR,  "nonconforming stabilized Crouzeix-Raviart finite elements" }
    };

    // momentum solver space scheme type -> name
    static std::map<std::string, spaceScheme> momSpaceSchemeNotation =
    {
        {"continous finite elements", CFE,                                    },
        {"nonconforming stabilized Crouzeix-Raviart finite elements", stabCR  }
    };

    // pressure parametrization -> name
    static std::map<pressParam, std::string> momPressParamName =
    {
        {clas,  "classical"  },
        {repl,  "replacement"}
    };

    // pressure name -> parametrization
    static std::map<std::string, pressParam> momPressParamNotation =
    {
        {"classical",   clas},
        {"replacement", repl}
    };

    // bc type -> name
    static std::map<bcType, std::string> momBcTypeName =
    {
        {noslip,  "no-slip (u = 0 on BND)"                             },
        {slip,    "slip (u|_n = 0 on BND)"                             },
        {fric,    "friction (+ friction in mom. eq. for BND elements)" }
    };

    // bc name -> type
    static std::map<std::string, bcType> momBcTypeNotation =
    {
        {"no-slip" , noslip },
        {"slip"    , slip   },
        {"friction", fric   }
    };
}