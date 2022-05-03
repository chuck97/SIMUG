#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include "defines.hpp"

#if defined(USE_MPI)
#include "mpi.h"
#endif

namespace SIMUG
{
    enum fwriteMode
    {
        rewrite,
        append
    };

    class Logger
    {
    private:
        std::ostream* os = NULL;
        std::fstream* ofs = NULL;

    public:
        inline Logger(std::ostream& os_)
        {
            os = &os_;
        };

        inline Logger(const std::string& filename_, fwriteMode mode)
        {
            ofs = new std::fstream;

            if (mode == fwriteMode::rewrite)
                ofs->open(filename_, std::ios::out);
            else 
                ofs->open(filename_, std::ios::app);

            if (!ofs->is_open())
            {
                SIMUG_ERR("can't open file "+ filename_ + " for logging!");
            }
        };

        inline Logger(std::ostream& os_, const std::string& filename_, fwriteMode mode)
        {
            os = &os_;
            ofs = new std::fstream;

            if (mode == fwriteMode::rewrite)
                ofs->open(filename_, std::ios::out);
            else 
                ofs->open(filename_, std::ios::app);

            if (!ofs->is_open())
            {
                SIMUG_ERR("can't open file "+ filename_ + " for logging!");
            }
        };

        inline void Log(const std::string& message)
        {
            if (os)
                *os << message;

            if (ofs)
                *ofs << message;
        };

        inline ~Logger()
        {
            if (ofs)
            {
                ofs->close(); 
                delete ofs;
            }

            if (os)
                os = NULL;
        };
    };
}