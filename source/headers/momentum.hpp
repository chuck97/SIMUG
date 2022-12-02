#pragma once

#include "mesh.hpp"
#include "dynvar.hpp"
#include "vecmath.hpp"
#include "local_assembling.hpp"
#include "tools.hpp"

namespace SIMUG
{
    class MomentumSolver
    {
    // special type binding for scalar and vector tags
    public:
        typedef INMOST::Tag scalar_tag;
        typedef INMOST::Tag velocity_tag; 
        
    // constructor
    public:
        MomentumSolver(SIMUG::IceMesh* mesh_,
                       double time_step_,
                       velocity_tag vel_tag_,
                       scalar_tag mass_tag_,
                       scalar_tag conc_tag_,
                       scalar_tag thick_tag_,
                       SIMUG::dyn::timeScheme mom_time_scheme_,
                       SIMUG::dyn::spaceScheme mom_space_scheme_,
                       SIMUG::dyn::pressParam mom_press_param_=SIMUG::dyn::pressParam::clas,
                       SIMUG::dyn::bcType mom_bc_type_=SIMUG::dyn::bcType::noslip);
    
    // manually set time step
    public:
        inline void SetTimeStep(double time_step_)
        {time_step = time_step_;};
    
    // virtual function for velocity computation
    protected:
        virtual void ComputeVelocity() = 0;

    // virtual function for profiling
    protected:
        virtual void PrintProfiling() = 0;
    
    // common parameters
    protected:
        SIMUG::IceMesh* mesh;
        double time_step;
        velocity_tag vel_tag;
        scalar_tag mass_tag;
        scalar_tag conc_tag;
        scalar_tag thick_tag;
        SIMUG::dyn::timeScheme mom_time_scheme;
        SIMUG::dyn::spaceScheme mom_space_scheme;
        SIMUG::dyn::pressParam mom_press_param;
        SIMUG::dyn::bcType mom_bc_type;
    };

    class AgridMomentumSolver: public MomentumSolver
    {
    // constructor
    public:
        AgridMomentumSolver(SIMUG::IceMesh* mesh_,
                            double time_step_,
                            velocity_tag vel_tag_,
                            scalar_tag mass_tag_,
                            scalar_tag conc_tag_,
                            scalar_tag thick_tag_,
                            SIMUG::dyn::timeScheme mom_time_scheme_,
                            SIMUG::dyn::pressParam mom_press_param_,
                            SIMUG::dyn::bcType mom_bc_type_,
                            const std::vector<double>& real_params_,
                            const std::vector<int>& integer_params_);
    // main computation process
    protected:
        virtual void ComputeVelocity() = 0;
        virtual void PrintProfiling() = 0;
    
    // auxilary parameters
    protected:

    //  solver parameters
        std::vector<double> real_params;
        std::vector<int> integer_params;
    };

    class Agrid_mEVP_Solver : public AgridMomentumSolver
    {
    // constructor
    public:
        Agrid_mEVP_Solver(SIMUG::IceMesh* mesh_,
                          double time_step_,
                          velocity_tag vel_tag_,
                          scalar_tag mass_tag_,
                          scalar_tag conc_tag_,
                          scalar_tag thick_tag_,
                          SIMUG::dyn::pressParam mom_press_param_=SIMUG::dyn::pressParam::clas,
                          SIMUG::dyn::bcType mom_bc_type_=SIMUG::dyn::bcType::noslip,
                          const std::vector<double>& real_params_=std::vector<double>{500.0, 500.0},
                          const std::vector<int>& integer_params_=std::vector<int>{500});
    public:
        void ComputeVelocity() override;

    private:
        void PrintProfiling() override;

    };
}