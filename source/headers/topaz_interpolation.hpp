#pragma once
#include "forcing.hpp"
#include "interpolation2d.hpp"
#include "ncfile_info.hpp"

namespace SIMUG
{
void TopazScalarInterpolation(SIMUG::IceMesh* mesh,
                              const NcFileInfo& file_info,
                              int netcdf_index,
                              const std::string& nc_variable_name,
                              double max_abs_value,
                              double invalid_value_fill,
                              double no_extrapolation_fill,
                              INMOST::Tag scalar_tag,
                              INMOST::Tag netcdf_coords_tag);
}