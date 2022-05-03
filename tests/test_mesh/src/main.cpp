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

    if (test_mesh_load())
    {
        if (rank == 0)
            std::cout << "Mesh loading test: OK!\n";
    }
    else
        SIMUG_ERR("Mesh loading test: FAILED!\n");

    if (test_bnd_selection())
    {
        if (rank == 0)
            std::cout << "Mesh bnd elems selection test: OK!\n";
    }
    else
        SIMUG_ERR("Mesh bnd elems selection test: FAILED!\n");

    if (test_id())
    {
        if (rank == 0)
            std::cout << "Mesh ID test: OK!\n";
    }
    else
        SIMUG_ERR("Mesh ID test: FAILED!\n");

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}
