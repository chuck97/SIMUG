#pragma once

#include "defines.hpp"

namespace SIMUG
{

double bilinear_interpolation(double x1, double x2, 
							  double y1, double y2, 
							  double x , double y,
							  double f_x1_y1,
							  double f_x1_y2,
							  double f_x2_y1,
							  double f_x2_y2);
}