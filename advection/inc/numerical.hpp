#pragma once
#include "defines.hpp"
#include "function.hpp"
#include "vecmath.hpp"

namespace SIMUG
{
    std::vector<double> FromBariToOrdinary(const std::vector<std::vector<double>>& trnodes,
                                           const std::vector<double>& bcoords);

    double integral_over_triangle(const std::vector<std::vector<double>>& trnodes,
                                  ScalarFunction f);
}