#pragma once
#include "simug.hpp"

#define VELOCITY_SCALE_FACTOR 0.00000625 // 40km/R_eth
#define FINAL_TIME 1080000.0             // 300 hours
#define EARTH_RADIUS 6400000.0           // 6400 km

std::vector<double> init_mass_gaussian(std::pair<double, double> coords, double time)
{
    double lon = coords.first;
    double lat = coords.second;

    double lon0 = 3.0*M_PI/4.0;
    double lat0 = 0.0;

    double lon1 = 5.0*M_PI/4.0;
    double lat1 = 0.0;

    double x = std::cos(lat)*std::cos(lon);
    double y = std::cos(lat)*std::sin(lon);
    double z = std::sin(lat);

    double x0 = std::cos(lat0)*std::cos(lon0);
    double y0 = std::cos(lat0)*std::sin(lon0);
    double z0 = std::sin(lat0);

    double x1 = std::cos(lat1)*std::cos(lon1);
    double y1 = std::cos(lat1)*std::sin(lon1);
    double z1 = std::sin(lat1);

    double width = 5.0;
    double ampl = 1.0;

    double h0 = ampl*std::exp( -width*
                                ((x - x0)*(x - x0) +
                                 (y - y0)*(y - y0) +
                                 (z - z0)*(z - z0)));

    double h1 = ampl*std::exp( -width*
                                ((x - x1)*(x - x1) +
                                 (y - y1)*(y - y1) +
                                 (z - z1)*(z - z1)));
    return {h0 + h1};
}

std::vector<double> init_mass_slotted_cylinders(std::pair<double, double> coords, double time)
{
    double lon = coords.first;
    double lat = coords.second;

    double first_center_lon = 3*M_PI/4.0;
    double first_center_lat = 0.0;

    double second_center_lon = 5*M_PI/4.0;
    double second_center_lat = 0.0;

    double background = 0.0;
    double scale_factor = 1.0;
    double radius = 0.5;

    double r1 = std::acos(std::sin(first_center_lat)*std::sin(lat) +
                          std::cos(first_center_lat)*std::cos(lat)*
                          std::cos(lon - first_center_lon));
        
    double r2 = std::acos(std::sin(second_center_lat)*std::sin(lat) +
                          std::cos(second_center_lat)*std::cos(lat)*
                          std::cos(lon - second_center_lon));

                          

    if (((r1 <= radius) and
         (std::abs(lon - first_center_lon) >= radius/6.0)) or
         ((r2 <= radius) and
         (std::abs(lon - second_center_lon) >= radius/6.0)))
    {
        return {scale_factor};
    }
    else if (((r1 <= radius) and
             ((std::abs(lon - first_center_lon) < radius/6.0))
             and ((lat - first_center_lat) < -(5.0/12.0)*radius)) or
             ((r2 <= radius) and
             ((std::abs(lon - second_center_lon) < radius/6.0))
             and ((lat - second_center_lat) > (5.0/12.0)*radius)))
    {
        return {scale_factor};
    }
    else
    {
        return {background};
    }
    return {0.0};
}

std::vector<double> non_div_velocity_1(std::pair<double, double> coords, double time)
{
    double lon = coords.first;
    double lat = coords.second;

    double u = VELOCITY_SCALE_FACTOR*
               std::sin(lon)*std::sin(lon)* 
               std::sin(lat*2.0)*
               std::cos(M_PI*time/(FINAL_TIME));

    double v = (VELOCITY_SCALE_FACTOR)*
               std::sin(2.0*lon)* 
               std::cos(lat)*
               std::cos(M_PI*time/(FINAL_TIME));

    return {u, v};
}

std::vector<double> non_div_velocity_2(std::pair<double, double> coords, double time)
{
    double lon = coords.first;
    double lat = coords.second;

    double lonn = lon - 2.0*M_PI*time/FINAL_TIME;

    double u = VELOCITY_SCALE_FACTOR*
               std::sin(lonn)*std::sin(lonn)* 
               std::sin(lat*2.0)*
               std::cos(M_PI*time/FINAL_TIME) +
               2*M_PI*std::cos(lat)/(FINAL_TIME);

    double v = VELOCITY_SCALE_FACTOR*
               std::sin(2.0*lonn)* 
               std::cos(lat)*
               std::cos(M_PI*time/FINAL_TIME);
    
    return {u, v};
}

std::vector<double> div_velocity(std::pair<double, double> coords, double time)
{
    double lon = coords.first;
    double lat = coords.second;

    double u = -VELOCITY_SCALE_FACTOR*
                std::sin(lon/2.0)*std::sin(lon/2.0)* 
                std::sin(lat*2.0)*
                std::cos(lat)*std::cos(lat)*
                std::cos(M_PI*time/FINAL_TIME);

    double v = (VELOCITY_SCALE_FACTOR/2.0)*
               std::sin(lon)* 
               std::cos(lat)*std::cos(lat)*std::cos(lat)*
               std::cos(M_PI*time/FINAL_TIME);

    return {u, v};
}

