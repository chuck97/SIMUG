#pragma once
#include <initializer_list>
#include <map>

// all mesh data variables
namespace SIMUG::mesh
{
    // types of surface
    enum surfType
    {
        plane,
        sphere,
        basin
    };

    // types of grid
    enum gridType
    {
        Agrid,
        Bgrid,
        Cgrid
    };

    // surface type notation -> name
    static const std::map<surfType, std::string> surfTypeName =
    {
        {plane,  "plane"      },
        {sphere, "sphere"     },
        {basin,  "water area" }
    };

    // surface type name -> notation
    static const std::map<std::string, surfType> surfTypeNotation =
    {
        {"plane",      plane  },
        {"sphere",     sphere },
        {"water area", basin  }
    };

    // grid type notation -> name
    static const std::map<gridType, std::string> gridTypeName =
    {
        {Agrid, "A-grid (scalars and vectors at nodes)"           },
        {Bgrid, "B-grid (scalars at nodes, vectors on triangles)" },
        {Cgrid, "C-grid (scalars on triangles, vectors on edges)" }
    };

    // surface type name -> notation
    static const std::map<std::string, gridType> gridTypeNotation =
    {
        {"A-grid", Agrid },
        {"B-grid", Bgrid },
        {"C-grid", Cgrid }
    };
} 