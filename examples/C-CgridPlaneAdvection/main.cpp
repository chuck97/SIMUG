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


#define PLANE_PATH "../../../SIMUG_v0/MESHES/pmf/Box_low_res.pmf"
#define FINAL_TIME 60.0

using namespace INMOST;
using namespace SIMUG;
using namespace std;

std::vector<double> init_mass_discont(std::pair<double, double> coords, double time)
{
    double x = coords.first;
    double y = coords.second;

    double x0 = 0.25;
    double y0 = 0.5;
    double r0 = 0.15;

    double x1 = 0.75;
    double y1 = 0.5;
    double r1 = 0.15;

    double dist0 = sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0));
    double dist1 = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1));

    double scale_factor = 1.0;
    double background = 0.0;

    if ((dist0 <= r0) or (dist1 <= r1))
    {
        return {scale_factor};
    }
    else
    {
        return {background};
    }
    return {0.0};
}

std::vector<double> velocity(std::pair<double, double> coords, double time)
{
    double x = coords.first;
    double y = coords.second;

    double xc = 0.5;
    double yc = 0.5;

    std::vector<double> radius_vec = {(x - xc), (y - yc), 0.0};
    std::vector<double> unit_normal = {0.0, 0.0, 1.0};

    double radius = L2_norm_vec(radius_vec);

    double scale_factor = 2*M_PI*radius/FINAL_TIME;

    std::vector<double> vel_vec = (unit_normal%radius_vec)*(1.0/L2_norm_vec(unit_normal%radius_vec))*scale_factor;

    return {vel_vec[0], vel_vec[1]};
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
    IceMesh mesh_plane(PLANE_PATH, mesh::surfType::plane, mesh::gridType::Cgrid); 
 
    // forcing initialization
    Forcing test_forcing(&mesh_plane);

    
    // assign initial scalar and velocity field 
    test_forcing.SetAnalytical(mesh::meshVar::mi, 0, init_mass_discont);
    test_forcing.Update(mesh::meshVar::mi, 0, coord::coordType::cart, 0.0);
    test_forcing.SetAnalytical(mesh::meshVar::ui, velocity);
    test_forcing.Update(mesh::meshVar::ui, coord::coordType::cart, 0.0);


    // get mass and velocity tag
    INMOST::Tag mass_tag = mesh_plane.GetDataMulti(mesh::gridElemType::Trian, 0)->Get(mesh::meshVar::mi);
    INMOST::Tag vel_tag = mesh_plane.GetDataSingle(mesh::gridElemType::Edge)->Get(mesh::meshVar::ui);

    // initialize advection solver
    CgridAdvectionSolver advection(&mesh_plane,
                                   1.0,
                                   vel_tag,
                                   adv::timeScheme::Euler,
                                   adv::spaceScheme::MUSCL,
                                   adv::advFilter::Minmod,
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
    if (mesh_plane.GetMesh()->GetProcessorRank() == 0)
    {
        logger.Log("Courant number: " + std::to_string(check_courant) +"\n");
        logger.Log("Time step: " + std::to_string(time_step) + " s\n");
        logger.Log("Number of steps: " + std::to_string(n_steps) + "\n");
    }
    BARRIER

    // write initial state to file
    mesh_plane.SaveVTU("CgridAdvectionPlane", 0);

    // time stepping
    double current_time = 0.0;

    for (size_t stepnum = 1; stepnum < n_steps; ++stepnum)
    {   
        // update time
        current_time += time_step;

        // update forcing
        test_forcing.Update(mesh::meshVar::ui, coord::coordType::cart, current_time);
        
        // transport scalars
        advection.TransportScalars();
        
        // write output to file
        if (stepnum % 1 == 0)
        {
            if (mesh_plane.GetMesh()->GetProcessorRank() == 0)
                logger.Log("Step: " + std::to_string(stepnum) +"\n");
                
            mesh_plane.SaveVTU("CgridAdvectionPlane", stepnum);
        }
        BARRIER
    }
	return 0;
}