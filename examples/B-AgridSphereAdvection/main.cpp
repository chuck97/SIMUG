#include "defines.hpp"
#include "inmost.h"
#include "model_var.hpp"
#include "mesh_info.hpp"
#include "data.hpp"
#include "mesh.hpp"
#include "forcing.hpp"
#include "advection.hpp"

#include <sstream>
#include <iomanip>

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

    double h1 = ampl*std::exp( -width*
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

int main(int argc, char* argv[])
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // initialize logger
    SIMUG::Logger logger(std::cout);

    // parse courant number from comand line
    double courant = 0.0;
    if (argc != 2)
    {
        SIMUG_ERR("One argument (desirable Courant number) should be given!");
    }
    else
    {
        istringstream ss(argv[1]);
        ss >> courant;
    }

    // mesh initialization
    IceMesh mesh_sphere(SPHERE_PATH, mesh::surfType::sphere, mesh::gridType::Agrid); 
 
    // forcing initialization
    Forcing test_forcing(&mesh_sphere);

 
    // assign initial scalar and velocity field 
    test_forcing.SetAnalytical(mesh::meshVar::mi, 0, init_mass_gaussian);
    test_forcing.Update(mesh::meshVar::mi, 0, coord::coordType::geo, 0.0);
    test_forcing.SetAnalytical(mesh::meshVar::ui, reversible_velocity);
    test_forcing.Update(mesh::meshVar::ui, coord::coordType::geo, 0.0);


    // initialize SLAE solver
    INMOST::Solver slae_solver("inner_ilu2"); 
    slae_solver.SetParameter("absolute_tolerance", "1e-9");

    // get mass and velocity tag
    INMOST::Tag mass_tag = mesh_sphere.GetDataMulti(mesh::gridElemType::Node, 0)->Get(mesh::meshVar::mi);
    INMOST::Tag vel_tag = mesh_sphere.GetDataSingle(mesh::gridElemType::Node)->Get(mesh::meshVar::ui);


    // initialize advection solver
    AgridAdvectionSolver advection(&mesh_sphere,
                                   1.0, 
                                   vel_tag,
                                   &slae_solver,
                                   adv::timeScheme::TG2,
                                   adv::spaceScheme::CFE,
                                   adv::advFilter::none,
                                   {});
    
    // add mass scalar for advection
    advection.AddScalar(mass_tag);
    
    // compute time step for desirable Courant number
    double time_step = courant/advection.GetMaxUdivDx();

    // setup computed time step
    advection.SetTimeStep(time_step); 

    // compute new courant number to check
    double check_courant = advection.GetMaxCourant();

    // compute number of steps
    size_t n_steps = (size_t)(FINAL_TIME/time_step);

    // log Courant number, time step and number of steps
    if (mesh_sphere.GetMesh()->GetProcessorRank() == 0)
    {
        logger.Log("Courant number: " + std::to_string(check_courant) + "\n");
        logger.Log("Time step: " + std::to_string(time_step) + "\n");
        logger.Log("Number of steps: " + std::to_string(n_steps) + "\n");
    }
    BARRIER

    // write initial state to file
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(5) << 0;
    mesh_sphere.SaveVTU("AgridAdvectionSphere" + ss.str());

/*
    // time stepping
    for (size_t stepnum = 0; stepnum < n_steps; ++stepnum)
    {   
        
        test_forcing.Update(mesh::meshVar::ui, coord::coordType::geo, current_time);
        advection.TransportScalars();
        current_time += time_step;
        
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(5) << stepnum;

        mesh_sphere.SaveVTU("advection_test" + ss.str());

    }
*/
	return 0;
}