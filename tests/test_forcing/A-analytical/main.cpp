#include "defines.hpp"
#include "inmost.h"
#include "model_var.hpp"
#include "mesh_info.hpp"
#include "data.hpp"
#include "mesh.hpp"
#include "forcing.hpp"

#ifdef USE_MPI
#include "mpi.h"
#endif


#define PLANE_PATH "../../../../SIMUG_v0/MESHES/pmf/Box_low_res.pmf"

using namespace INMOST;
using namespace SIMUG;
using namespace std;

std::vector<double> sin_mass(std::pair<double, double> coords, double time)
{
    double x = coords.first;
    double y = coords.second;

    return {sin(x*2.0*M_PI)*sin(y*2.0*M_PI)};
}

std::vector<double> sin_velocity(std::pair<double, double> coords, double time)
{
    double x = coords.first;
    double y = coords.second;

    return {sin(x*2.0*M_PI)*sin(y*2.0*M_PI), cos(x*2.0*M_PI)*cos(y*2.0*M_PI), 0.0};
}

bool test_analytical_forcing()
{
    // make plane mesh
    IceMesh mesh_plane(PLANE_PATH, mesh::surfType::plane, mesh::gridType::Agrid);   

    // make simple forcing class
    Forcing test_forcing(&mesh_plane);

    // set analytical mass
    test_forcing.SetAnalytical(mesh::meshVar::mi, 0, sin_mass);

    // update mass value according to analytical function
    test_forcing.Update(mesh::meshVar::mi, 0, coord::coordType::model, 0.0);

    // create temporarily grid vector on triangles
    mesh_plane.GetDataSingle(mesh::gridElemType::Trian)->Create("velocity", 3, INMOST::DATA_REAL);

    // set analytical velocity
    test_forcing.SetAnalytical("velocity", mesh::gridElemType::Trian, sin_velocity);

    // update velocity value according to analytical function
    test_forcing.Update("velocity", mesh::gridElemType::Trian, coord::coordType::model, 0.0);

    // save mesh
    mesh_plane.SaveVTU("forcing_test");

    return true;
}

int main()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (test_analytical_forcing())
    {
        if (rank == 0)
            std::cout << "Analytical forcing test: OK!\n";
    }
    else
        SIMUG_ERR("Analytical forcing test: FAILED!\n");

    BARRIER

#ifdef USE_MPI
    MPI_Finalize();
#endif
	return 0;
}