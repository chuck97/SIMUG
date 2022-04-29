#include "tests.hpp"

int main()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (test_data())
    {
        if (rank == 0)
            std::cout << "MeshDataVar test: OK!\n";
    }
    else
        SIMUG_ERR("MeshDataVar test: FAILED!\n");

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}
