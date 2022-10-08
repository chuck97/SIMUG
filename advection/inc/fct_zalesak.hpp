#pragma once

#include "inmost.h"
#include "defines.hpp"
#include "vecmath.hpp"

namespace SIMUG
{
    void FCT_Procedure(INMOST::Mesh* mesh,
                       const std::vector<std::vector<std::vector<double>>>& M_C_minus_M_L,
                       INMOST::Tag scal_tag,
                       INMOST::Tag scal_low_tag,
                       INMOST::Tag scal_high_tag,
                       INMOST::Tag node_id_tag,
                       double fct_cd);

    void CalculateAlpha(INMOST::Mesh* mesh,
                        INMOST::Tag Fe_trian_tag,
                        INMOST::Tag alpha_trian_tag,
                        INMOST::Tag scal_tag,
                        INMOST::Tag scal_low_tag,
                        INMOST::Tag node_id_tag);
}