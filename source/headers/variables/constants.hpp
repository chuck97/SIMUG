#pragma once

// all physical constants
namespace SIMUG
{
    struct GenConsts
    {
        static constexpr double g = 9.8;          // gravity acceleration value (m s-1)
        static constexpr double f = 1.46e-4;      // Coriolis parameter value (s-1)
        static constexpr double Reth = 6.371e6;   // Earth radius (m)
    };

    struct AirConsts
    {
        static constexpr double rhoa = 1.28;      // air density value (kg m-3)
    };

    struct WaterConsts
    {
        static constexpr double rhow = 1023.0;   // water density value (kg m-3)
    };

    struct IceConsts
    {
        static constexpr double rhoi = 900.0;     // ice density value (kg m-3)
        static constexpr double Cw = 5.5e-3;      // water-ice drag coefficient value (-)
        static constexpr double Ca = 1.2e-3;      // water-air drag coefficient value (-)
        static constexpr double C = 20.0;         // ice concentration parameter in pressure defenition value(-)
        static constexpr double pstr = 27.5e3;    // p* coefficient in pressure defenition value (N m-2)
        static constexpr double e = 2.0;          // eccentricity of VP elliptic yield curvee (-)
        static constexpr double delmin = 2e-9;    // minimal value of delta dunctional (s-1) 
        static constexpr double amin = 1e-6;      // minimal value of ice concentration (-)
        static constexpr double hmin = 1e-2;      // minimalvalue of mean ice thickness (m)
    };
}