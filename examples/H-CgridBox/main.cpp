#include "simug.hpp"
#include "inmost.h"

#include<sstream>
#include<iomanip>

#ifdef USE_MPI
#include "mpi.h"
#endif

#define SQUARE_PATH "../../../SIMUG_v0/MESHES/pmf/square_1000km_15km.pmf" 

#define BOX_LEN_SCALE 1000e3
#define TIME_SCALE 345600.0

using namespace INMOST;
using namespace SIMUG;
using namespace std;


// function for initial ice thickness
std::vector<double> init_ice_thickness(std::pair<double, double> coords, double time)
{
    return {2.0};
}

// function for initial ice concentration
std::vector<double> init_ice_concentration(std::pair<double, double> coords, double time)
{
    double x = coords.first;
    double y = coords.second;

    return {x/BOX_LEN_SCALE};
}

// function for initial ice velocity
std::vector<double> init_ice_velocity(std::pair<double, double> coords, double time)
{
    return {0.0, 0.0};
} 

// function for wind velocity
std::vector<double> wind_velocity(std::pair<double, double> coords, double time)
{
    double x = coords.first;
    double y = coords.second;

    return 
    {
        5.0 + (std::sin(2.0*M_PI*time/TIME_SCALE) - 3.0)*std::sin(2.0*M_PI*x/BOX_LEN_SCALE)*std::sin(2.0*M_PI*y/BOX_LEN_SCALE),
        5.0 + (std::sin(2.0*M_PI*time/TIME_SCALE) - 3.0)*std::sin(2.0*M_PI*y/BOX_LEN_SCALE)*std::sin(2.0*M_PI*x/BOX_LEN_SCALE)
    };
}

// function for ocean velocity
std::vector<double> ocean_velocity(std::pair<double, double> coords, double time)
{
    double x = coords.first;
    double y = coords.second;

    return 
    {
        0.1*(2.0*y - BOX_LEN_SCALE)/BOX_LEN_SCALE,
        -0.1*(2.0*x - BOX_LEN_SCALE)/BOX_LEN_SCALE 
    };
}

// function for ocean level
std::vector<double> ocean_level(std::pair<double, double> coords, double time)
{
    return {0.0};
}

void run_model(double time_step,
               double total_time,
               const std::string& mesh_path,
               int output_frequency,
               adv::timeScheme advection_time_scheme,
               adv::spaceScheme advection_space_scheme,
               adv::advFilter advection_filter,
               const std::vector<double>& advection_params,
               const std::vector<double>& mevp_real_params,
               const std::vector<int>& mevp_integer_params,
               const std::string& output_prefix)
{
    // start MPI activity
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // initialize logger
    SIMUG::Logger logger(std::cout);

    // mesh initialization
    IceMesh* mesh_plane = new IceMesh(mesh_path, mesh::surfType::plane, mesh::gridType::Cgrid); 
 
    // forcing initialization
    Forcing forcing(mesh_plane);
    
    // assign initial ice thickness
    forcing.SetAnalytical
    (
        mesh::meshVar::hi,
        0,
        [](std::pair<double, double> coords, double time)->std::vector<double>
        {
            return init_ice_thickness(coords, time);
        }
    );
    forcing.Update(mesh::meshVar::hi, 0, coord::coordType::cart, 0.0);

    // assign initial ice concentration
    forcing.SetAnalytical
    (
        mesh::meshVar::ai,
        0,
        [](std::pair<double, double> coords, double time)->std::vector<double>
        {
            return init_ice_concentration(coords, time);
        }
    );
    forcing.Update(mesh::meshVar::ai, 0, coord::coordType::cart, 0.0);


    // assign initial ice velocity
    forcing.SetAnalytical
    (
        mesh::meshVar::ui,
        [](std::pair<double, double> coords, double time)->std::vector<double>
        {
            return init_ice_velocity(coords, time);
        }
    );
    forcing.Update(mesh::meshVar::ui, coord::coordType::cart, 0.0);

    // set analytical air velocity
    forcing.SetAnalytical
    (
        mesh::meshVar::ua,
        [](std::pair<double, double> coords, double time)->std::vector<double>
        {
            return wind_velocity(coords, time);
        }
    );
    forcing.Update(mesh::meshVar::ua, coord::coordType::cart, 0.0);

    // set analytical ocean velocity
    forcing.SetAnalytical
    (
        mesh::meshVar::uw,
        [](std::pair<double, double> coords, double time)->std::vector<double>
        {
            return ocean_velocity(coords, time);
        }
    );
    forcing.Update(mesh::meshVar::uw, coord::coordType::cart, 0.0);

    // set analytical ocean level
    forcing.SetAnalytical
    (
        mesh::meshVar::hw,
        [](std::pair<double, double> coords, double time)->std::vector<double>
        {
            return ocean_level(coords, time);
        }
    );
    forcing.Update(mesh::meshVar::hw, coord::coordType::cart, 0.0);


    // get tags for ice values
    INMOST::Tag conc_tag = mesh_plane->GetDataMulti(mesh::gridElemType::Trian, 0)->Get(mesh::meshVar::ai);
    INMOST::Tag thick_tag = mesh_plane->GetDataMulti(mesh::gridElemType::Trian, 0)->Get(mesh::meshVar::hi);
    INMOST::Tag vel_tag = mesh_plane->GetDataSingle(mesh::gridElemType::Edge)->Get(mesh::meshVar::ui);


    // initialize advection solver
    CgridAdvectionSolver advection(mesh_plane,
                                   time_step,
                                   vel_tag,
                                   advection_time_scheme,
                                   advection_space_scheme,
                                   advection_filter,
                                   advection_params);

    // add mass scalar for advection
    advection.AddScalar(thick_tag);
    advection.AddScalar(conc_tag);

    // compute number of steps
    size_t n_steps = (size_t)(total_time/time_step);

    // initialize momentum solver
    Cgrid_mEVP_Solver momentum(mesh_plane,
                               time_step,
                               vel_tag,
                               conc_tag,
                               thick_tag,
                               SIMUG::dyn::pressParam::clas,
                               SIMUG::dyn::bcType::noslip,
                               mevp_real_params,
                               mevp_integer_params
                               );


    // log time step and number of steps
    if (mesh_plane->GetMesh()->GetProcessorRank() == 0)
    {
        logger.Log("Time step: " + std::to_string(time_step) + " s\n");
        logger.Log("Number of steps: " + std::to_string(n_steps) + "\n");
    }
    BARRIER

    // write initial state to file
    mesh_plane->SaveVTU(output_prefix, 0);

    // time stepping
    double current_time = 0.0;

    for (size_t stepnum = 1; stepnum < (n_steps + 1); ++stepnum)
    {   
        // update time
        current_time += time_step;

        // update air velocty
        forcing.Update(mesh::meshVar::ua, coord::coordType::cart, current_time);

        // update ocean velocty
        forcing.Update(mesh::meshVar::uw, coord::coordType::cart, current_time);

        // update ocean level
        forcing.Update(mesh::meshVar::hw, coord::coordType::cart, current_time);
        
        // transport scalars
        advection.TransportScalars();

        // update ice velocity
        momentum.ComputeVelocity();
        
        // write output to file
        if ((stepnum % output_frequency == 0) or (stepnum == (n_steps-1)))
        {       
            mesh_plane->SaveVTU(output_prefix, stepnum);
        }
        BARRIER
    }
    
    //delete slae solver and mesh
    delete mesh_plane;

    // end MPI activity
#ifdef USE_MPI
    MPI_Finalize();
#endif

}

int main()
{
    run_model(600.0,                                // time step (seconds)
              604800,                               // total time (seconds)
              (std::string)SQUARE_PATH,             // path to mesh (LOW_RES_MESH_PATH, MIDDLE_RES_MESH_PATH, HIGH_RES_MESH_PATH) 
              1,                                    // output frequncy n (every n-th time will be written to file)
              adv::timeScheme::TRK2,                // advection time scheme (TG2, TTG2, TTG3, TTG4)
              adv::spaceScheme::FVupwind,           // advection space scheme (FVupwind, MUST, MUSCL)
              adv::advFilter::none,                 // advection filter (Minmod, VanLeer, Superbee, BarthJesperson, none)
              {},                                   // vector pf advection filter parameters ()
              {500.0, 500.0, 1.0},                  // vector of mEVP real params and stab param (alpha, beta, alpha_stab)
              {200},                                // vector of mEVP integer params (number of interations)
              "CgridBox"                            // output prefix
              );
}