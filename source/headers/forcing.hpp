#pragma once

#include "mesh.hpp"
#include "netcdf.h"
#include "ncfile_info.hpp"
#include "mesh_data.hpp"
#include "topaz_interpolation.hpp"
#include "cams_interpolation.hpp"

// Log netcdf error message 
#define NC_ERR(e) {std::cerr << "Error:" <<  nc_strerror(e) << std::endl; exit(1);}

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

        // assign file for mesh variable
        void SetFile(mesh::meshVar mesh_var,
                     NcFileInfo file_info
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
        
        // update mesh variable from file interpolation
        void UpdateTOPAZ(mesh::meshVar mesh_var,
                         const std::vector<std::string>& nc_variable_name,
                         double max_abs_value,
                         double invalid_value_fill,
                         double no_extrapolation_fill,
                         int index);

        void UpdateTOPAZ(mesh::meshVar mesh_var,
                         int layer,
                         const std::vector<std::string>& nc_variable_name,
                         double max_abs_value,
                         double invalid_value_fill,
                         double no_extrapolation_fill,
                         int index);
        
        void UpdateCAMS(mesh::meshVar mesh_var,
                        const std::vector<std::string>& nc_variable_names,
                        double max_abs_value,
                        double invalid_value_fill,
                        double no_extrapolation_fill,
                        int index);

    private:
        SIMUG::IceMesh* mesh;

        void UpdateVariable(INMOST::Tag var_tag,
                            mesh::gridElemType elem_type,
                            coord::coordType coord_type,
                            FuncPtr func_ptr,
                            double time);
        
        std::map<mesh::meshVar, NcFileInfo> file_info_map;
    };
}