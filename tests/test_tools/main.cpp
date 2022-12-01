#include <iostream>
#include <iomanip>
#include <variant>
#include <unistd.h>

#include "simug.hpp"

#if defined(USE_MPI)
#include "mpi.h"
#endif

bool trim_test()
{
    std::string a = " abcd ";
    std::string b = " abcd ";
    std::string c = " abcd ";

    SIMUG::tools::ltrim(a);
    SIMUG::tools::rtrim(b);
    SIMUG::tools::trim(c);
    
    if (a != "abcd ")
        return false;

    if (b != " abcd")
        return false;

    if (c != "abcd")
        return false;
    
    return true;
}

bool log_test()
{
    int rank = 0;

#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if (rank == 0)
    {
        {
        SIMUG::Logger l1(std::cout);
        l1.Log("log console \n");
        }
        {
        SIMUG::Logger l2("test.txt", SIMUG::fwriteMode::rewrite);
        l2.Log("log file \n");
        }
        {
        SIMUG::Logger l3(std::cout, "test.txt", SIMUG::fwriteMode::append);
        l3.Log("log console and file \n");
        }
    }

    BARRIER
    return true;
};

bool timer_test()
{
    unsigned int microsecond = 1000000;
    int rank = 0;
    int size = 1;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    SIMUG::Timer test_timer;

    test_timer.Launch();

    for (int i = 0; i < size; ++i)
    {
        if (rank == i)
        {
            usleep((i+1)*microsecond);
        }
    }

    test_timer.Stop();

    double time = test_timer.GetTime();
    double min_time = test_timer.GetMinTime();
    double max_time = test_timer.GetMaxTime();
    double avg_time = test_timer.GetAvgTime();

    if (rank == 0)
    {
        std::cout << std::setprecision(15) << "time: "     << time    << " ms" << std::endl;
        std::cout << std::setprecision(15) << "min time: " << min_time << " ms" << std::endl;
        std::cout << std::setprecision(15) << "max time: " << max_time << " ms" << std::endl;
        std::cout << std::setprecision(15) << "avg time: " << avg_time << " ms" << std::endl;
    }

    test_timer.Reset();
    test_timer.Launch();

   for (int i = 0; i < size; ++i)
    {
        if (rank == i)
        {
            usleep((i+1)*microsecond);
        }
    }

    test_timer.Stop();

    time = test_timer.GetTime();
    min_time = test_timer.GetMinTime();
    max_time = test_timer.GetMaxTime();
    avg_time = test_timer.GetAvgTime();

    if (rank == 0)
    {
        std::cout << std::setprecision(15) << "time: "     << time     << " ms" << std::endl;
        std::cout << std::setprecision(15) << "min time: " << min_time << " ms" << std::endl;
        std::cout << std::setprecision(15) <<  "max time: " << max_time << " ms" << std::endl;
        std::cout << std::setprecision(15) << "avg time: " << avg_time << " ms" << std::endl;
    }

    BARRIER

    return true;
}

int main(int argc, char *argv[])
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (trim_test())
    {
        if (rank == 0)
            std::cout << "Trim test: OK!\n";
    }
    else
        SIMUG_ERR("Trim test: FAILED!\n")

    if (log_test())
    {
        if (rank == 0)
            std::cout << "Logger test: OK!\n";
    }
    else
        SIMUG_ERR("Logger test: FAILED!\n")

    if (timer_test())
    {
        if (rank == 0)
            std::cout << "Timer test: OK!\n";
    }
    else
        SIMUG_ERR("Timer test: FAILED!\n")
    return 0;

#ifdef USE_MPI
    MPI_Finalize();
#endif
}