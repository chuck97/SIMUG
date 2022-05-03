#pragma once

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

#if defined(USE_MPI)
#define SIMUG_ERR(message) {std::cerr << "Error: " << message  << std::endl; MPI_Finalize(); exit(1);}
#else
#define SIMUG_ERR(message) {std::cerr << "Error: " << message  << std::endl; exit(1);}
#endif