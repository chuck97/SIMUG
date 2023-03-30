#include "simug.hpp"
#include "inmost.h"
//#include "parser.hpp"

#include<sstream>
#include<iomanip>

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace INMOST;
using namespace SIMUG;
using namespace std;

#define MESH_ARCTIC "/home/users/spetrov/SIMUG/SIMUG_v0/MESHES/pmf/Arctic.pmf"

/*
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
    IceMesh* mesh_sphere = new IceMesh(mesh_path, output_folder, mesh::surfType::sphere, mesh::gridType::Agrid); 
 
    // forcing initialization
    Forcing forcing(mesh_sphere);
    
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
    forcing.Update(mesh::meshVar::hi, 0, coord::coordType::geo, 0.0);

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
    forcing.Update(mesh::meshVar::ai, 0, coord::coordType::geo, 0.0);


    // assign initial ice velocity
    forcing.SetAnalytical
    (
        mesh::meshVar::ui,
        [](std::pair<double, double> coords, double time)->std::vector<double>
        {
            return init_ice_velocity(coords, time);
        }
    );
    forcing.Update(mesh::meshVar::ui, coord::coordType::geo, 0.0);

    // set analytical air velocity
    forcing.SetAnalytical
    (
        mesh::meshVar::ua,
        [](std::pair<double, double> coords, double time)->std::vector<double>
        {
            return wind_velocity(coords, time);
        }
    );
    forcing.Update(mesh::meshVar::ua, coord::coordType::geo, 0.0);

    // set analytical ocean velocity
    forcing.SetAnalytical
    (
        mesh::meshVar::uw,
        [](std::pair<double, double> coords, double time)->std::vector<double>
        {
            return ocean_velocity(coords, time);
        }
    );
    forcing.Update(mesh::meshVar::uw, coord::coordType::geo, 0.0);

    // set analytical ocean level
    forcing.SetAnalytical
    (
        mesh::meshVar::hw,
        [](std::pair<double, double> coords, double time)->std::vector<double>
        {
            return ocean_level(coords, time);
        }
    );
    forcing.Update(mesh::meshVar::hw, coord::coordType::geo, 0.0);


    // get tags for ice values
    INMOST::Tag conc_tag = mesh_sphere->GetDataMulti(mesh::gridElemType::Node, 0)->Get(mesh::meshVar::ai);
    INMOST::Tag thick_tag = mesh_sphere->GetDataMulti(mesh::gridElemType::Node, 0)->Get(mesh::meshVar::hi);
    INMOST::Tag vel_tag = mesh_sphere->GetDataSingle(mesh::gridElemType::Node)->Get(mesh::meshVar::ui);

    // initialize SLAE solver
    INMOST::Solver* slae_solver = new INMOST::Solver("inner_ilu2"); 
    slae_solver->SetParameter("absolute_tolerance", "1e-9");

    // initialize advection solver
    AgridAdvectionSolver advection(mesh_sphere,
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
    Agrid_mEVP_Solver momentum(mesh_sphere,
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
    if (mesh_sphere->GetMesh()->GetProcessorRank() == 0)
    {
        logger.Log("Time step: " + std::to_string(time_step) + " s\n");
        logger.Log("Number of steps: " + std::to_string(n_steps) + "\n");
    }
    BARRIER

    // write initial state to file
    mesh_sphere->SaveVTU(output_prefix, 0);

    // time stepping
    double current_time = 0.0;

    for (size_t stepnum = 1; stepnum < (n_steps + 1); ++stepnum)
    {   
        // update time
        current_time += time_step;

        if (mesh_sphere->GetMesh()->GetProcessorRank() == 0)
        {
            logger.Log("\n!!! Step " + std::to_string(stepnum) + " out of " + std::to_string(n_steps) + " !!!\n");
        }

        // update air velocty
        forcing.Update(mesh::meshVar::ua, coord::coordType::geo, current_time);

        // update ocean velocty
        forcing.Update(mesh::meshVar::uw, coord::coordType::geo, current_time);

        // update ocean level
        forcing.Update(mesh::meshVar::hw, coord::coordType::geo, current_time);
        
        // transport scalars
        if (is_advection)
        {
            advection.TransportScalars();
        }

        // update ice velocity
        momentum.ComputeVelocity();
        
        // write output to file
        if ((stepnum % output_frequency == 0) or (stepnum == (n_steps)))
        {       
            mesh_sphere->SaveVTU(output_prefix, stepnum);
        }
        BARRIER
    }
    
    //delete slae solver and mesh
    delete mesh_sphere;
    delete slae_solver;
}
*/

int main(int argc, char* argv[])
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
    IceMesh* mesh_arctic = new IceMesh(MESH_ARCTIC, mesh::surfType::basin, mesh::gridType::Cgrid); 

    // forcing initialization
    Forcing topaz_forcing(mesh_arctic);
    Forcing cams_forcing(mesh_arctic);

    // concentration file
    topaz_forcing.SetFile
    (
        mesh::meshVar::ai,
        {
            "/home/users/spetrov/nc_topaz/dataset-ice.nc",
            "x",
            "y",
            "longitude",
            "latitude",
            "scale_factor",
            "missing_value",
            "add_offset",
            false
        }
    );

    // thickness file
    topaz_forcing.SetFile
    (
        mesh::meshVar::hi,
        {
            "/home/users/spetrov/nc_topaz/dataset-ice.nc",
            "x",
            "y",
            "longitude",
            "latitude",
            "scale_factor",
            "missing_value",
            "add_offset",
            false
        }
    );

    // ssh file
    topaz_forcing.SetFile
    (
        mesh::meshVar::hw, 
        {
            "/home/users/spetrov/nc_topaz/dataset-ice.nc",
            "x",
            "y",
            "longitude",
            "latitude",
            "",
            "missing_value",
            "",
            false
        }
    );

    // ice velocity file
    topaz_forcing.SetFile
    (
        mesh::meshVar::uw, 
        {
            "/home/users/spetrov/nc_topaz/dataset-ice.nc",
            "x",
            "y",
            "longitude",
            "latitude",
            "",
            "missing_value",
            "",
            false
        }
    );

    // water velocity file
    topaz_forcing.SetFile
    (
        mesh::meshVar::ui, 
        {
            "/home/users/spetrov/nc_topaz/dataset-ice.nc",
            "x",
            "y",
            "longitude",
            "latitude",
            "scale_factor",
            "missing_value",
            "add_offset",
            false
        }
    );

    // air velocity file
    cams_forcing.SetFile
    (
        mesh::meshVar::ua, 
        {
            "/home/users/spetrov/nc_topaz/atm-forcing-2019-rev.nc",
            "",
            "",
            "longitude",
            "latitude",
            "scale_factor",
            "missing_value",
            "add_offset",
            false
        }
    );

    // initialize ice concentration using topaz
    topaz_forcing.UpdateTOPAZ(mesh::meshVar::ai, 0, {"fice"}, 1.1, 0.0, 0.0, 0);

    std::cout << "ai ok" << std::endl;

    // initialize ice thickness using topaz
    topaz_forcing.UpdateTOPAZ(mesh::meshVar::hi, 0, {"hice"}, 50.0, 0.0, 0.0, 0);

    std::cout << "hi ok" << std::endl;

    // initialize ice thickness using topaz
    topaz_forcing.UpdateTOPAZ(mesh::meshVar::hw, {"ssh"}, 10.0, 0.0, 0.0, 0);

    std::cout << "hw ok" << std::endl;

    // initialize ice velocity using topaz
    topaz_forcing.UpdateTOPAZ(mesh::meshVar::ui, {"uice", "vice"}, 10.0, 0.0, 0.0, 0);

    std::cout << "ui ok" << std::endl;

    // initialize ice velocity using topaz
    topaz_forcing.UpdateTOPAZ(mesh::meshVar::uw, {"u", "v"}, 10.0, 0.0, 0.0, 0);

    std::cout << "uw ok" << std::endl;

    // initialize atm velocity using cams
    cams_forcing.UpdateCAMS(mesh::meshVar::ua, {"u10", "v10"}, 100.0, 0.0, 0.0, 0);

    std::cout << "ua ok" << std::endl;

    // save mesh
    mesh_arctic->SaveVTU("CgridArctic", 0);

    // delete mesh
    delete mesh_arctic;

    // end MPI activity
#ifdef USE_MPI
    MPI_Finalize();
#endif
}