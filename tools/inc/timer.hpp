#pragma once
#include <chrono>
#include "defines.hpp"
#include <variant>

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace SIMUG
{
    class Timer
    {
    public:
        typedef std::chrono::high_resolution_clock clock;
        typedef std::chrono::duration<double, std::milli> duration;

    public:
        std::variant<double, clock::time_point> start_time;
        std::variant<double, clock::time_point> stop_time;

    public:
        Timer();
        void Reset();
        void Launch();
        void Stop();
        double GetTime() const;
        double GetMinTime() const;
        double GetMaxTime() const;
        double GetAvgTime() const;
    };
}