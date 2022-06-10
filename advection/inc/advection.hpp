#pragma once

#include "mesh.hpp"
#include "advvar.hpp"
#include "vecmath.hpp"
#include "local_assembling.hpp"

namespace SIMUG
{
    class AdvectionSolver
    {

    // special type binding for scalar and vector tags
    public:
        typedef INMOST::Tag scalar_tag;
        typedef INMOST::Tag velocity_tag; 

    // constructor
    public:
        inline AdvectionSolver(SIMUG::IceMesh* mesh_,
                               velocity_tag vel_tag_,
                               double time_step_):
            mesh(mesh_),
            vel_tag(vel_tag_),
            time_step(time_step_)
        {};
        
    // tags management and manual setting the time step
    public:
        inline void SetTimeStep(double time_step_)
        {time_step = time_step_;};

        inline void AddScalar(scalar_tag input_scal_tag)
        {scal_tags.push_back(input_scal_tag);};

    
        inline void AddScalar(const std::vector<scalar_tag>& input_scal_tags)
        {for (const auto& sc_tag: input_scal_tags){scal_tags.push_back(sc_tag);}};

    // virtual function for transporting single scalar
    protected:
        virtual void Evaluate(velocity_tag vel_tag, scalar_tag scal_tag) = 0;
    
    // virtual function for profiling
    protected:
        virtual void PrintProfiling() = 0;

    // main function for transporting all scalar
    public:
        void TransportScalars();

    // auxilary functions for Courant number calculation
    public:
        virtual double GetMaxUdivDx() = 0;
        virtual double GetMaxCourant() = 0;

    // configuration parameters
    protected:
        SIMUG::IceMesh* mesh;
        velocity_tag vel_tag;
        double time_step;
        std::vector<scalar_tag> scal_tags;
    };

    class AgridAdvectionSolver: public AdvectionSolver
    {
    // constructor
    public:
        AgridAdvectionSolver(SIMUG::IceMesh* mesh_,
                             double time_step_,
                             velocity_tag vel_tag_,
                             INMOST::Solver* slae_solver_,
                             SIMUG::adv::timeScheme adv_time_scheme_,
                             SIMUG::adv::spaceScheme adv_space_scheme_,
                             SIMUG::adv::advFilter adv_filter_,
                             const std::vector<double>& params_);
        
    // auxilary functions for Courant number calculation
    public:
        double GetMaxUdivDx() override;
        double GetMaxCourant() override;

    // main evaluate function
    private:
        void Evaluate(velocity_tag vel_tag, scalar_tag scal_tag) override;

    // auxiliary functions for assembling global LHS and RHS
    private:
        enum StepNumber {first, second};
        void AssembleLHS();
        void AssembleSingleStepRHS(velocity_tag vel_tag, scalar_tag scal_tag);
        void AssembleDoubleStepRHS(velocity_tag vel_tag, scalar_tag scal_tag, scalar_tag scal_half_tag,  StepNumber step_num);

    // configuration parameters
    private:
        INMOST::Solver* slae_solver;
        adv::timeScheme adv_time_scheme;
        adv::spaceScheme adv_space_scheme;
        adv::advFilter adv_filter;
        std::vector<double> params;
    
    // auxilary data for filter
    private:
        std::vector<std::vector<std::vector<double>>> M_C_minus_M_L;

    // sparse matricies and vectors 
    private:
        INMOST::Sparse::Matrix LHS;
        INMOST::Sparse::Matrix LHS_low;

        INMOST::Sparse::Vector RHS;
        INMOST::Sparse::Vector RHS_low;

    // time for profiling
    private:
        void PrintProfiling() override;
        double RHS_assembling_time = 0.0;
        double limiter_time = 0.0;
        double matrix_invertion_time = 0.0;

    };

/*
    class CgridAdvectionSolver : public AdvectionSolver
    {
    public:
        inline CgridAdvectionSolver(SIMUG::IceMesh* mesh_,
                                    double time_step_):
            AdvectionSolver(mesh_, time_step_)
        {};

    private:
        void Evaluate(velocity_tag vel_tag, scalar_tag scal_tag) override;
    };
*/
}