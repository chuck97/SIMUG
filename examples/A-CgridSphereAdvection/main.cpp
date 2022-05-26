#include "defines.hpp"
#include "inmost.h"
#include "model_var.hpp"
#include "mesh_info.hpp"
#include "data.hpp"
#include "mesh.hpp"
#include "forcing.hpp"
#include "advection.hpp"

#include<sstream>
#include<iomanip>

#ifdef USE_MPI
#include "mpi.h"
#endif


#define SPHERE_PATH "../../../SIMUG_v0/MESHES/pmf/Sphere.pmf"
#define FINAL_TIME 300.0*3600.0
#define VELOCITY_SCALE_FACTOR 80.0/6400000.0

using namespace INMOST;
using namespace SIMUG;
using namespace std;

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

    double h1 = ampl*std::exp(-width*
                                ((x - x1)*(x - x1) +
                                 (y - y1)*(y - y1) +
                                 (z - z1)*(z - z1)));
    return {h0 + h1};
}

std::vector<double> reversible_velocity(std::pair<double, double> coords, double time)
{
    double lon = coords.first;
    double lat = coords.second;

    double u = VELOCITY_SCALE_FACTOR*
                   std::sin(lon/2.0)*std::sin(lon/2.0)* 
                   std::sin(lat*2.0)*
                   std::cos(M_PI*time/FINAL_TIME);

    double v = (VELOCITY_SCALE_FACTOR/2.0)*
               std::sin(lon)* 
               std::cos(lat)*
               std::cos(M_PI*time/FINAL_TIME);

    return {u, v};
}

int main()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    IceMesh mesh_sphere(SPHERE_PATH, mesh::surfType::sphere, mesh::gridType::Cgrid); 
 
    Forcing test_forcing(&mesh_sphere);

    
    test_forcing.SetAnalytical(mesh::meshVar::mi, init_mass_gaussian);
    test_forcing.Update(mesh::meshVar::mi, coord::coordType::geo, 0.0);
    test_forcing.SetAnalytical(mesh::meshVar::ui, reversible_velocity);

    double time_step = 0.2*3600.0;
    //double time_step = 50.0*3600.0;
    double current_time = 0.0;

    CgridAdvectionSolver advection(&mesh_sphere, time_step);

    // get velocity tag
    INMOST::Tag mass_tag = mesh_sphere.GetProgData(mesh::gridElemType::Trian, 0)->Get(mesh::meshVar::mi);
    INMOST::Tag vel_tag = mesh_sphere.GetForcData(mesh::gridElemType::Edge)->Get(mesh::meshVar::ui);

    advection.AddScalar(vel_tag, mass_tag);

    size_t n_steps = (size_t)(FINAL_TIME/time_step);

    for (size_t stepnum = 0; stepnum < n_steps; ++stepnum)
    {   
        test_forcing.Update(mesh::meshVar::ui, coord::coordType::geo, current_time);
        advection.TransportScalars();
        current_time += time_step;
        
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(5) << stepnum;

        mesh_sphere.SaveVTU("advection_test" + ss.str());
    }

	return 0;
}