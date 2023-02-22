#include "simug.hpp"
#include "inmost.h"
#include "parser.hpp"

#include<sstream>
#include<iomanip>

#ifdef USE_MPI
#include "mpi.h"
#endif

#define LKF_BOX_LEN_SCALE 512e3

using namespace INMOST;
using namespace SIMUG;
using namespace std;


// function for initial ice thickness
std::vector<double> init_ice_thickness(std::pair<double, double> coords, double time)
{
    return {0.3};
}

// function for initial ice concentration
std::vector<double> init_ice_concentration(std::pair<double, double> coords, double time)
{
    return {1.0};
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

    // current time in hours
    double t  = time/3600.0;

    // position of the center of cyclone
    double m_x = 0.5*LKF_BOX_LEN_SCALE + 0.1*LKF_BOX_LEN_SCALE*(t/24.0);
    double m_y = m_x;

    // angle between veclone velocity and radius
    double alpha = 72.0*(M_PI/180.0);

    // maximal air velocity
    double v_air_max = 15.0;

    // distance to cyclone center
    double r = std::sqrt((m_x-x)*(m_x-x) + (m_y-y)*(m_y-y));

    // reduction factor
    double s = (1.0/50.0)*std::exp(-(r*1e-5));

    
    return 
    {
        -s*v_air_max*(std::cos(alpha)*(x-m_x)*1e-3 + std::sin(alpha)*(y-m_y)*1e-3),
        -s*v_air_max*(-std::sin(alpha)*(x-m_x)*1e-3 + std::cos(alpha)*(y-m_y)*1e-3)
    };
}

// function for ocean velocity
std::vector<double> ocean_velocity(std::pair<double, double> coords, double time)
{
    double x = coords.first;
    double y = coords.second;

    double vel_scale = 0.01;

    std::vector<double> res = 
    {
        (2.0*y - 1.0*LKF_BOX_LEN_SCALE) / (1.0*LKF_BOX_LEN_SCALE),
       -(2.0*x - 1.0*LKF_BOX_LEN_SCALE) / (1.0*LKF_BOX_LEN_SCALE) 
    };

    return vel_scale*res;
}

// function for ocean level
std::vector<double> ocean_level(std::pair<double, double> coords, double time)
{
    return {0.0};
}

void run_model(double time_step,                             // time step (seconds)
               double total_time,                            // total time (seconds)
               const std::string& mesh_path,                 // path to mesh .pmf file
               int output_frequency,                         // output frequency
               bool is_advection,                            // turn on the advection?
               adv::timeScheme advection_time_scheme,        // advection time scheme
               adv::spaceScheme advection_space_scheme,      // advection space scheme
               adv::advFilter advection_filter,              // advection filter type
               const std::vector<double>& advection_params,  // vector of real advection params
               const std::vector<double>& mevp_real_params,  // vector of real momentum params
               const std::vector<int>& mevp_integer_params,  // vector of integer momentum params
               const std::string& output_prefix,             // output prefix
               const std::string& output_folder)             // output folder
{
    // initialize logger
    SIMUG::Logger logger(std::cout);

    // mesh initialization
    IceMesh* mesh_plane = new IceMesh(mesh_path, output_folder, mesh::surfType::plane, mesh::gridType::Agrid); 
 
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
    INMOST::Tag conc_tag = mesh_plane->GetDataMulti(mesh::gridElemType::Node, 0)->Get(mesh::meshVar::ai);
    INMOST::Tag thick_tag = mesh_plane->GetDataMulti(mesh::gridElemType::Node, 0)->Get(mesh::meshVar::hi);
    INMOST::Tag vel_tag = mesh_plane->GetDataSingle(mesh::gridElemType::Node)->Get(mesh::meshVar::ui);

    // initialize SLAE solver
    INMOST::Solver* slae_solver = new INMOST::Solver("inner_ilu2"); 
    slae_solver->SetParameter("absolute_tolerance", "1e-9");

    // initialize advection solver
    AgridAdvectionSolver advection(mesh_plane,
                                   time_step,
                                   vel_tag,
                                   slae_solver,
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
    Agrid_mEVP_Solver momentum(mesh_plane,
                               time_step,
                               vel_tag,
                               conc_tag,
                               thick_tag,
                               SIMUG::dyn::pressParam::clas,
                               SIMUG::dyn::bcType::noslip,
                               mevp_real_params,
                               mevp_integer_params
                               );

    // log Courant number, time step and number of steps
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

        if (mesh_plane->GetMesh()->GetProcessorRank() == 0)
        {
            logger.Log("\n!!! Step " + std::to_string(stepnum) + " out of " + std::to_string(n_steps) + " !!!\n");
        }

        // update air velocty
        forcing.Update(mesh::meshVar::ua, coord::coordType::cart, current_time);

        // update ocean velocty
        forcing.Update(mesh::meshVar::uw, coord::coordType::cart, current_time);

        // update ocean level
        forcing.Update(mesh::meshVar::hw, coord::coordType::cart, current_time);
        
        // transport scalars
        if (is_advection)
        {
            advection.TransportScalars();
        }

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
    delete slae_solver;
}

int main(int argc, char* argv[])
{
// start MPI activity
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // parse input
    std::string current_exec_name = argv[0]; 
    std::vector<std::string> all_args;

    std::string input_external_path;

    if (argc > 1) 
    {
        all_args.assign(argv + 1, argv + argc);
    }

    if (all_args.size() > 1)
    {
        SIMUG_ERR("should be only one configuration \".txt\" file on input");
    }
    else if (all_args.size() == 0)
    {
        SIMUG_ERR("configuration \".txt\" file should be given on input");
    }

    Parser config(all_args.back());

    run_model(config.time_step_seconds,
              config.total_time_seconds,
              config.grid_file,
              config.output_frequency,
              (config.is_advection == 1) ? true : false,
              adv::timeScheme::TTG2,    // advection time scheme (TG2, TTG2, TTG3, TTG4)
              adv::spaceScheme::CFE,    // advection space scheme (CFE)
              adv::advFilter::Zalesak,  // advection filter (Zalesak, None)
              {0.5},                    // vector of advection filter parameters (fct cd value)
              {
                  config.alpha_mEVP,
                  config.beta_mEVP
              },
              {config.Nits_mEVP},
              config.output_prefix,
              config.output_dir
              );

// end MPI activity
#ifdef USE_MPI
    MPI_Finalize();
#endif
}