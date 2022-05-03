#include "timer.hpp"

using namespace SIMUG;

Timer::Timer()
{
    Reset();
}

void Timer::Reset()
{
#ifdef USE_MPI
    double current_time = MPI_Wtime();
    start_time = current_time;
    stop_time = current_time;
#else      
    clock::time_point current_time = clock::now();
    start_time = current_time;
    stop_time = current_time;
#endif
};

void Timer::Launch()
{
#ifdef USE_MPI
    BARRIER
    start_time = MPI_Wtime();
#else
    start_time = clock::now();
#endif
};

void Timer::Stop()
{
#ifdef USE_MPI
    stop_time = MPI_Wtime();
#else
    stop_time = clock::now();
#endif
};

double Timer::GetTime() const
{
#ifdef USE_MPI
    return (std::get<double>(stop_time) - std::get<double>(start_time))*1e3;
#else
    duration elapsed = std::get<clock::time_point>(stop_time) - std::get<clock::time_point>(start_time);
    return elapsed.count();
#endif
};

double Timer::GetMinTime() const
{
#ifdef USE_MPI
    BARRIER
    double duration = GetTime();
    double duration_min = 0.0;

    MPI_Allreduce(&duration, &duration_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    return duration_min;

#else
    return GetTime();
#endif
};

double Timer::GetMaxTime() const
{
#ifdef USE_MPI
    BARRIER
    double duration = GetTime();
    double duration_max = 0.0;

    MPI_Allreduce(&duration, &duration_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return duration_max;

#else
    return GetTime();
#endif
};

double Timer::GetAvgTime() const 
{
#ifdef USE_MPI
    BARRIER
    double duration = GetTime();
    double duration_sum = 0.0;

    MPI_Allreduce(&duration, &duration_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    return duration_sum/nprocs;

#else
    return GetTime();
#endif
};