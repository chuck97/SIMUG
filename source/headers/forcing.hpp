#pragma once

#include "mesh.hpp"

namespace SIMUG
{
    class Forcing
    {
    public:
        typedef std::vector<double>(*FuncPtr)(std::pair<double, double>, double);
    
    public:
        inline Forcing(SIMUG::IceMesh* mesh_):
            mesh(mesh_)
        {};
        
        // assign analytical expression for mesh variable
        void SetAnalytical(mesh::meshVar mesh_var,
                           FuncPtr func_ptr
                           );

        // assign analytical expression for mesh variable (with layer)
        void SetAnalytical(mesh::meshVar mesh_var,
                           int layer,
                           FuncPtr func_ptr
                           );
        
        // assign analytical expression for arbitrary variable
        void SetAnalytical(const std::string& varname, 
                           mesh::gridElemType elem_type,
                           FuncPtr func_ptr
                          );

        // assign analytical expression for arbitrary variable (with layer)
        void SetAnalytical(const std::string& varname, 
                           mesh::gridElemType elem_type,
                           int layer,
                           FuncPtr func_ptr
                          );

        // update mesh variable according to analytical expression
        void Update(mesh::meshVar mesh_var,
                    coord::coordType coord_type,
                    double time);

        void Update(mesh::meshVar mesh_var,
                    int layer,
                    coord::coordType coord_type,
                    double time);

        void Update(const std::string& varname,
                    mesh::gridElemType elem_type,
                    coord::coordType coord_type,
                    double time);
        
        void Update(const std::string& varname,
                    mesh::gridElemType elem_type,
                    int layer,
                    coord::coordType coord_type,
                    double time);

    protected:
        SIMUG::IceMesh* mesh;

        void UpdateVariable(INMOST::Tag var_tag,
                            mesh::gridElemType elem_type,
                            coord::coordType coord_type,
                            FuncPtr func_ptr,
                            double time);
    };

    class IceForcing : public Forcing
    {
    public:
        inline IceForcing(IceMesh* mesh_):
            Forcing(mesh_)
        {};

    };

    class AirForcing : public Forcing
    {
        inline AirForcing(IceMesh* mesh_):
            Forcing(mesh_)
        {};

    };

    class WaterForcing : public Forcing
    {
        inline WaterForcing(IceMesh* mesh_):
            Forcing(mesh_)
        {};

    };
}