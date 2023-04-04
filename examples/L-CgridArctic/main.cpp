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


void run_model(double time_step,                             // time step (seconds)
               int output_frequency,                         // output frequency
               adv::timeScheme advection_time_scheme,        // advection time scheme
               adv::spaceScheme advection_space_scheme,      // advection space scheme
               adv::advFilter advection_filter,              // advection filter type
               const std::vector<double>& advection_params,  // vector of real advection params
               const std::vector<double>& mevp_real_params,  // vector of real momentum params
               const std::vector<int>& mevp_integer_params)  // vector of integer momentum params
{
    // initialize logger
    SIMUG::Logger logger(std::cout);

    // mesh initialization
    IceMesh* mesh_arctic = new IceMesh(MESH_ARCTIC, mesh::surfType::basin, mesh::gridType::Cgrid); 

    // forcing initialization
    Forcing topaz_forcing(mesh_arctic);
    Forcing cams_forcing(mesh_arctic);

    // prepare concentration file
    topaz_forcing.SetFile
    (
        mesh::meshVar::ai,
        {
            "/home/users/spetrov/Forcing/dataset-topaz4-arc_2020_03_01-2020_03_04_1hr.nc",
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

    // prepare thickness file
    topaz_forcing.SetFile
    (
        mesh::meshVar::hi,
        {
            "/home/users/spetrov/Forcing/dataset-topaz4-arc_2020_03_01-2020_03_04_1hr.nc",
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

    // prepare ssh file
    topaz_forcing.SetFile
    (
        mesh::meshVar::hw, 
        {
            "/home/users/spetrov/Forcing/dataset-topaz4-arc_2020_03_01-2020_03_04_1hr.nc",
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
        mesh::meshVar::uw, 
        {
            "/home/users/spetrov/Forcing/dataset-topaz4-arc_2020_03_01-2020_03_04_1hr.nc",
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


    // prepare air velocity file
    cams_forcing.SetFile
    (
        mesh::meshVar::ua, 
        {
            "/home/users/spetrov/Forcing/cams_2021_03_01-2021_03_04_3h_3hr.nc",
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

    // initialize ice thickness using topaz
    topaz_forcing.UpdateTOPAZ(mesh::meshVar::hi, 0, {"hice"}, 8.0, 0.0, 0.0, 0);

    // initialize ice thickness using topaz
    topaz_forcing.UpdateTOPAZ(mesh::meshVar::hw, {"ssh"}, 1.0, 0.0, 0.0, 0);

    // initialize water velocity using topaz
    topaz_forcing.UpdateTOPAZ(mesh::meshVar::uw, {"u", "v"}, 2.0, 0.0, 0.0, 0);

    // initialize atm velocity using cams
    cams_forcing.UpdateCAMS(mesh::meshVar::ua, {"u10", "v10"}, 100.0, 0.0, 0.0, 0);

    if (mesh_arctic->GetMesh()->GetProcessorRank() == 0)
    {
        logger.Log("\n Forcing initialization done! \n");
    }
    BARRIER

    // get tags for ice values
    INMOST::Tag conc_tag = mesh_arctic->GetDataMulti(mesh::gridElemType::Trian, 0)->Get(mesh::meshVar::ai);
    INMOST::Tag thick_tag = mesh_arctic->GetDataMulti(mesh::gridElemType::Trian, 0)->Get(mesh::meshVar::hi);
    INMOST::Tag vel_tag = mesh_arctic->GetDataSingle(mesh::gridElemType::Edge)->Get(mesh::meshVar::ui);

    // initialize advection solver
    CgridAdvectionSolver advection(mesh_arctic,
                                   time_step,
                                   vel_tag,
                                   advection_time_scheme,
                                   advection_space_scheme,
                                   advection_filter,
                                   advection_params);
    
    // add concentration and thickness scalar for advection
    advection.AddScalar(conc_tag);
    advection.AddScalar(thick_tag);

    if (mesh_arctic->GetMesh()->GetProcessorRank() == 0)
    {
        logger.Log("\n Advection initialization done! \n");
    }
    BARRIER

    // initialize momentum solver
    Cgrid_mEVP_Solver momentum(mesh_arctic,
                               time_step,
                               vel_tag,
                               conc_tag,
                               thick_tag,
                               SIMUG::dyn::pressParam::clas,
                               SIMUG::dyn::bcType::noslip,
                               mevp_real_params,
                               mevp_integer_params
                               );
    
    if (mesh_arctic->GetMesh()->GetProcessorRank() == 0)
    {
        logger.Log("\n Momentum initialization done! \n");
    }
    BARRIER
    
    // compute number of steps
    size_t n_steps = (size_t)(93.0*60.0*60.0/time_step);

    // compute atm and ocn forcing frequency
    size_t atm_freq = (size_t)(3.0*60.0*60.0/time_step);
    size_t ocn_freq = (size_t)(1.0*60.0*60.0/time_step);

    // log time step and number of steps
    if (mesh_arctic->GetMesh()->GetProcessorRank() == 0)
    {
        logger.Log("Time step: " + std::to_string(time_step) + " s\n");
        logger.Log("Number of steps: " + std::to_string(n_steps) + "\n");
    }
    BARRIER

    // save initial state
    mesh_arctic->SaveVTU("CgridArctic", 0);

    // time stepping
    double current_time = 0.0;
    size_t forcing_counter = 1;
    size_t atm_forcing_ind = 1;
    size_t ocn_forcing_ind = 1;

    for (size_t stepnum = 1; stepnum < (n_steps + 1); ++stepnum)
    {   
        // update time
        current_time += time_step;

        // log step
        if (mesh_arctic->GetMesh()->GetProcessorRank() == 0)
        {
            logger.Log("\n!!! Step " + std::to_string(stepnum) + " out of " + std::to_string(n_steps) + " !!!\n");
        }

        // update air velocty
        if ((forcing_counter%atm_freq) == 0)
        {
            cams_forcing.UpdateCAMS(mesh::meshVar::ua, {"u10", "v10"}, 100.0, 0.0, 0.0, atm_forcing_ind);
            ++atm_forcing_ind;
        }

        // update ocean velocty and ssh
        if ((forcing_counter%ocn_freq) == 0)
        {
            topaz_forcing.UpdateTOPAZ(mesh::meshVar::uw, {"u", "v"}, 10.0, 0.0, 0.0, ocn_forcing_ind);
            topaz_forcing.UpdateTOPAZ(mesh::meshVar::hw, {"ssh"}, 10.0, 0.0, 0.0, ocn_forcing_ind);
            ++ocn_forcing_ind;
        }
        ++forcing_counter;

        // transport scalars
        advection.TransportScalars();

        // update ice velocity
        momentum.ComputeVelocity();

        
        // write output to file
        if ((stepnum % output_frequency == 0) or (stepnum == (n_steps)))
        {       
            mesh_arctic->SaveVTU("CgridArctic", stepnum);
        }
        BARRIER
    }

    // delete mesh and slae solver
    delete mesh_arctic;

    BARRIER
}

int main(int argc, char* argv[])
{
// start MPI activity
    int rank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    
    run_model(300.0,
              50,
              adv::timeScheme::Euler,
              adv::spaceScheme::FVupwind,
              adv::advFilter::none,
              {}, 
              {2000.0, 2000.0},
              {200});

// end MPI activity
#ifdef USE_MPI
    MPI_Finalize();
#endif
}