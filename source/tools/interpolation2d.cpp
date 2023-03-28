#include "interpolation2d.hpp"

namespace SIMUG
{
double bilinear_interpolation(double x1, double x2, 
							  double y1, double y2, 
							  double x , double y,
							  double f_x1_y1,
							  double f_x1_y2,
							  double f_x2_y1,
							  double f_x2_y2)
{
	if ((x < x1) or (x > x2) or (y < y1) or (y > y2))
	{
		SIMUG_ERR("Bilinear interpolation failed - point is outside rectangle!");
	}
	double f_x_y1 = ((x2 - x)/(x2 - x1))*f_x1_y1 + ((x - x1)/(x2 - x1))*f_x2_y1; 
	double f_x_y2 = ((x2 - x)/(x2 - x1))*f_x1_y2 + ((x - x1)/(x2 - x1))*f_x2_y2; 
	double f_x_y  = ((y2 - y)/(y2 - y1))*f_x_y1 + ((y - y1)/(y2 - y1))*f_x_y2;   
	return f_x_y;
}
}