#pragma once
#include "forcing.hpp"
#include "interpolation2d.hpp"
#include "ncfile_info.hpp"

namespace SIMUG
{
void CamsVectorInterpolation(SIMUG::IceMesh* mesh,
                              const NcFileInfo& file_info,
                              int netcdf_index,
                              const std::string& nc_variable_u,
                              const std::string& nc_variable_v,
                              double max_abs_value,
                              double invalid_value_fill,
                              double no_extrapolation_fill,
                              INMOST::Tag vector_tag,
                              INMOST::Tag geo_coords_tag);

}