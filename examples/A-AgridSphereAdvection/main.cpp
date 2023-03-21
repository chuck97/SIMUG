#include "test_forcing.hpp"
#include "parser.hpp"

#include <sstream>
#include <iomanip>

#ifdef USE_MPI
#include "mpi.h"
#endif


using namespace INMOST;
using namespace SIMUG;
using namespace std;

double ComputeRelativeL2Error(IceMesh* ice_mesh, INMOST::Tag mass_tag, auto scalar_func)
{
    INMOST::Mesh* mesh = ice_mesh->GetMesh();
    INMOST::Tag geo_coords_tag = ice_mesh->GetGridInfo(mesh::gridElemType::Node)->coords[coord::coordType::geo];
    INMOST::Tag trian_area_tag = ice_mesh->GetGridInfo(mesh::gridElemType::Trian)->GetCartesianSize();
    double local_sum_sqr_diff = 0.0;
    double local_sum_sqr_init = 0.0;

    for(auto trianit = mesh->BeginCell(); trianit != mesh->EndCell(); ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            double trian_area = trianit->Real(trian_area_tag);
            ElementArray<Node> adj_nodes = trianit->getNodes();
            for (int j = 0; j < 3; ++j)
            {
                double lon = adj_nodes[j]->RealArray(geo_coords_tag)[0];
                double lat = adj_nodes[j]->RealArray(geo_coords_tag)[1];
                double exact_value = scalar_func(std::pair<double, double>{lon, lat}, 0.0)[0];
                double computed_value = adj_nodes[j]->Real(mass_tag);
                local_sum_sqr_diff += ((computed_value - exact_value)*(computed_value - exact_value)/3.0)*trian_area;
                local_sum_sqr_init += ((exact_value*exact_value)/3.0)*trian_area;
            }
        }
    }
    BARRIER

    // compute sum of (exact - computed)**2 from all processes
    double all_sum_sqr_diff = local_sum_sqr_diff;

#if defined(USE_MPI)
   MPI_Allreduce(&local_sum_sqr_diff, &all_sum_sqr_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    BARRIER
    // compute sum of exact**2 from all processes
    double all_sum_sqr_init = local_sum_sqr_init;

#if defined(USE_MPI)
   MPI_Allreduce(&local_sum_sqr_init, &all_sum_sqr_init, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    BARRIER

    // compute relative L2_error
    double L2_error = std::sqrt(all_sum_sqr_diff) / std::sqrt(all_sum_sqr_init);

    return L2_error;
}

void run_model(double courant_number,                        // desired Courant number
               const std::string& mesh_path,                 // path to mesh .pmf file
               int output_frequency,                         // output frequency
               initScalar init_scalar,                       // initial scalar distridution
               velField vel_field,                           // analytical velocity field
               adv::timeScheme advection_time_scheme,        // advection time scheme
               adv::spaceScheme advection_space_scheme,      // advection space scheme
               adv::advFilter advection_filter,              // advection filter type
               const std::string& output_prefix,             // output prefix
               const std::string& output_folder)             // output folder
{
    // initialize logger
    SIMUG::Logger logger(std::cout);

    // mesh initialization
    IceMesh* mesh_sphere = new IceMesh(mesh_path, output_folder, mesh::surfType::sphere, mesh::gridType::Agrid); 
 
    // forcing initialization
    Forcing test_forcing(mesh_sphere);

    // figure out init scalar function and set m forcing
    auto scalar_func = init_mass_gaussian;
    if (init_scalar == initScalar::GH)
    {
        scalar_func = init_mass_gaussian;
    }
    else if (init_scalar == initScalar::SC)
    {
        scalar_func = init_mass_slotted_cylinders;
    }

    test_forcing.SetAnalytical(mesh::meshVar::mi, 0, scalar_func);
    test_forcing.Update(mesh::meshVar::mi, 0, coord::coordType::geo, 0.0);

    // figure out velocity field and set vel forcing
    auto vector_func = non_div_velocity_1;
    if (vel_field == velField::ND1)
    {
        vector_func = non_div_velocity_1;
    }
    else if (vel_field == velField::ND2)
    {
        vector_func = non_div_velocity_2;
    }
    else if (vel_field == velField::D)
    {
        vector_func = div_velocity;
    }

    test_forcing.SetAnalytical(mesh::meshVar::ui, vector_func);
    test_forcing.Update(mesh::meshVar::ui, coord::coordType::geo, 0.0);

    // get mass and velocity tag
    INMOST::Tag mass_tag = mesh_sphere->GetDataMulti(mesh::gridElemType::Node, 0)->Get(mesh::meshVar::mi);
    INMOST::Tag vel_tag = mesh_sphere->GetDataSingle(mesh::gridElemType::Node)->Get(mesh::meshVar::ui);

    // initialize SLAE solver
    INMOST::Solver* slae_solver = new INMOST::Solver("inner_ilu2"); 
    slae_solver->SetParameter("absolute_tolerance", "1e-9");

    // initialize advection solver
    AgridAdvectionSolver advection(mesh_sphere,
                                   1.0, 
                                   vel_tag,
                                   slae_solver,
                                   advection_time_scheme,
                                   advection_space_scheme,
                                   advection_filter,
                                   {0.5});
    
    // add mass scalar for advection
    advection.AddScalar(mass_tag);
    
    // compute time step for desirable Courant number
    double time_step = courant_number/advection.GetMaxUdivDx();

    // setup computed time step
    advection.SetTimeStep(time_step); 

    // compute new courant number to check
    double check_courant = advection.GetMaxCourant();

    // compute number of steps
    size_t n_steps = (size_t)(FINAL_TIME/time_step);

    // log Courant number, time step and number of steps
    if (mesh_sphere->GetMesh()->GetProcessorRank() == 0)
    {
        logger.Log("Courant number: " + std::to_string(check_courant) +"\n");
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

        // update forcing
        test_forcing.Update(mesh::meshVar::ui, coord::coordType::geo, current_time);

        // transport scalars
        advection.TransportScalars();
        
        // write output to file
        if ((stepnum % output_frequency == 0) or (stepnum == n_steps))
        {       
            mesh_sphere->SaveVTU(output_prefix, stepnum);
        }
        BARRIER
    }

    // compute error
    double error = ComputeRelativeL2Error(mesh_sphere, mass_tag, scalar_func);

    // write error
    if (mesh_sphere->GetMesh()->GetProcessorRank() == 0)
    {
        logger.Log("\n##### L2_ERROR = " + std::to_string(error) + " #####\n");
    }

    delete slae_solver;
    delete mesh_sphere;
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

    run_model(config.courant_number, 
              config.grid_file, 
              config.output_frequency, 
              config.init_scalar,
              config.vel_field, 
              config.time_scheme, 
              config.space_scheme, 
              config.adv_filter, 
              config.output_prefix,
              config.output_dir);

// end MPI activity
#ifdef USE_MPI
    MPI_Finalize();
#endif
}