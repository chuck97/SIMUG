#include "simug.hpp"

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

    double h1 = ampl*std::exp( -width*
                                ((x - x1)*(x - x1) +
                                 (y - y1)*(y - y1) +
                                 (z - z1)*(z - z1)));
    return {h0 + h1};
}

std::vector<double> init_mass_slotted_cylinders(std::pair<double, double> coords, double time)
{
    double lon = coords.first;
    double lat = coords.second;

    double first_center_lon = 3*M_PI/4.0;
    double first_center_lat = 0.0;

    double second_center_lon = 5*M_PI/4.0;
    double second_center_lat = 0.0;

    double background = 0.0;
    double scale_factor = 1.0;
    double radius = 0.5;

    double r1 = std::acos(std::sin(first_center_lat)*std::sin(lat) +
                          std::cos(first_center_lat)*std::cos(lat)*
                          std::cos(lon - first_center_lon));
        
    double r2 = std::acos(std::sin(second_center_lat)*std::sin(lat) +
                          std::cos(second_center_lat)*std::cos(lat)*
                          std::cos(lon - second_center_lon));

                          

    if (((r1 <= radius) and
         (std::abs(lon - first_center_lon) >= radius/6.0)) or
         ((r2 <= radius) and
         (std::abs(lon - second_center_lon) >= radius/6.0)))
    {
        return {scale_factor};
    }
    else if (((r1 <= radius) and
             ((std::abs(lon - first_center_lon) < radius/6.0))
             and ((lat - first_center_lat) < -(5.0/12.0)*radius)) or
             ((r2 <= radius) and
             ((std::abs(lon - second_center_lon) < radius/6.0))
             and ((lat - second_center_lat) > (5.0/12.0)*radius)))
    {
        return {scale_factor};
    }
    else
    {
        return {background};
    }
    return {0.0};
}

std::vector<double> reversible_velocity(std::pair<double, double> coords, double time)
{
    double lon = coords.first;
    double lat = coords.second;

    double u = VELOCITY_SCALE_FACTOR*
                   std::sin(lon/2.0)*std::sin(lon/2.0)* 
                   std::sin(lat*2.0)*
                   std::cos(M_PI*time/(FINAL_TIME));

    double v = (VELOCITY_SCALE_FACTOR/2.0)*
               std::sin(lon)* 
               std::cos(lat)*
               std::cos(M_PI*time/(FINAL_TIME));

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
    IceMesh mesh_sphere(SPHERE_PATH, mesh::surfType::sphere, mesh::gridType::Cgrid); 
 
    // forcing initialization
    Forcing test_forcing(&mesh_sphere);

    
    // assign initial scalar and velocity field 
    test_forcing.SetAnalytical(mesh::meshVar::mi, 0, init_mass_slotted_cylinders);
    test_forcing.Update(mesh::meshVar::mi, 0, coord::coordType::geo, 0.0);
    test_forcing.SetAnalytical(mesh::meshVar::ui, reversible_velocity);
    test_forcing.Update(mesh::meshVar::ui, coord::coordType::geo, 0.0);


    // get mass and velocity tag
    INMOST::Tag mass_tag = mesh_sphere.GetDataMulti(mesh::gridElemType::Trian, 0)->Get(mesh::meshVar::mi);
    INMOST::Tag vel_tag = mesh_sphere.GetDataSingle(mesh::gridElemType::Edge)->Get(mesh::meshVar::ui);

    // initialize advection solver
    CgridAdvectionSolver advection(&mesh_sphere,
                                   1.0,
                                   vel_tag,
                                   adv::timeScheme::Euler,
                                   adv::spaceScheme::MUSCL,
                                   adv::advFilter::Superbee,
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
        logger.Log("Courant number: " + std::to_string(check_courant) +"\n");
        logger.Log("Time step: " + std::to_string(time_step) + " s\n");
        logger.Log("Number of steps: " + std::to_string(n_steps) + "\n");
    }
    BARRIER

    // write initial state to file
    mesh_sphere.SaveVTU("CgridAdvectionSphere", 0);

    // time stepping
    double current_time = 0.0;

    for (size_t stepnum = 1; stepnum < n_steps; ++stepnum)
    {   
        // update time
        current_time += time_step;

        // update forcing
        test_forcing.Update(mesh::meshVar::ui, coord::coordType::geo, current_time);
        
        // transport scalars
        advection.TransportScalars();
        
        // write output to file
        if (stepnum % 1 == 0)
        {
            if (mesh_sphere.GetMesh()->GetProcessorRank() == 0)
                logger.Log("Step: " + std::to_string(stepnum) +"\n");
                
            mesh_sphere.SaveVTU("CgridAdvectionSphere", stepnum);
        }
        BARRIER
    }

	return 0;
}