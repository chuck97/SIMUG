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
                       scalar_tag scal_tag_,
                       SIMUG::dyn::scheme mom_time_scheme_,
                       SIMUG::dyn::press press_param_,
                       SIMUG::dyn::bc bc_type_);
    
    // manually set time step
    public:
        inline void SetTimeStep(double time_step_)
        {time_step = time_step_;};
    
    // virtual function for velocity computation
    public:
        virtual void ComputeVelocity() = 0;

    // virtual function for profiling
    protected:
        virtual void PrintProfiling() = 0;
    
    // common parameters
    protected:
        SIMUG::IceMesh* mesh;
        double time_step;
        velocity_tag vel_tag;
        scalar_tag scal_tag;
        SIMUG::dyn::scheme mom_time_scheme;
        SIMUG::dyn::press press_param;
        SIMUG::dyn::bc bc_type;
    };

    class AgridMomentumSolver: public MomentumSolver
    {
    // constructor
    public:
        AgridMomentumSolver(SIMUG::IceMesh* mesh_,
                            double time_step_,
                            velocity_tag vel_tag_,
                            scalar_tag scal_tag_,
                            SIMUG::dyn::scheme mom_time_scheme_,
                            SIMUG::dyn::press press_param_,
                            SIMUG::dyn::bc bc_type_,
                            const std::vector<double>& params_);
    // main computation process
    public:
        virtual void ComputeVelocity() override;
    
    // auxilary parameters
    private:
        // to do
    };
}