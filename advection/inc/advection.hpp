#pragma once

#include "mesh.hpp"
#include "advvar.hpp"
#include "vecmath.hpp"

namespace SIMUG
{
    class AdvectionSolver
    {
    public:
        typedef INMOST::Tag scalar_tag;
        typedef INMOST::Tag velocity_tag; 
    public:
        inline AdvectionSolver(SIMUG::IceMesh* mesh_,
                               velocity_tag vel_tag_,
                               double time_step_):
            mesh(mesh_),
            vel_tag(vel_tag_),
            time_step(time_step_)
        {};

        inline void SetTimeStep(double time_step_)
        {time_step = time_step_;};

        // add scalar for advection 
        inline void AddScalar(scalar_tag input_scal_tag)
        {scal_tags.push_back(input_scal_tag);};

        // add vector of scalars for advection 
        inline void AddScalar(const std::vector<scalar_tag>& input_scal_tags)
        {for (const auto& sc_tag: input_scal_tags){scal_tags.push_back(sc_tag);}};

    protected:
        // transport one scalar
        //virtual void Evaluate(velocity_tag vel_tag, scalar_tag scal_tag) = 0;

    public:
        // transport all scalars
        /*
        inline void TransportScalars()
        {
            for (auto& scal_tag: scal_tags)
            {    
                Evaluate(vel_tag, scal_tag);
            } 
            BARRIER
        };
        */

        // function for evaluating time step to assign for given Courant number 
        virtual double GetMaxUdivDx() = 0;

        // obtain max Courant number
        virtual double GetMaxCourant() = 0;

    protected:
        SIMUG::IceMesh* mesh;
        velocity_tag vel_tag;
        double time_step;
        std::vector<scalar_tag> scal_tags;
    };

    class AgridAdvectionSolver: public AdvectionSolver
    {
    public:
        AgridAdvectionSolver(SIMUG::IceMesh* mesh_,
                            double time_step_,
                            velocity_tag vel_tag_,
                            INMOST::Solver* slae_solver_,
                            SIMUG::adv::timeScheme adv_time_scheme_,
                            SIMUG::adv::spaceScheme adv_space_scheme_,
                            SIMUG::adv::advFilter adv_filter_,
                            const std::vector<double>& params_);
        
        // function for evaluating time step to assign for given Courant number 
        double GetMaxUdivDx() override;

        // obtain max Courant number
        double GetMaxCourant() override;

    private:
        //void Evaluate(velocity_tag vel_tag, scalar_tag scal_tag) override;
        void AssembleLHS();
        void AssembleRHS();

    private:
        INMOST::Solver* slae_solver;
        adv::timeScheme adv_time_scheme;
        adv::spaceScheme adv_space_scheme;
        adv::advFilter adv_filter;
        std::vector<double> params;

    private:
        INMOST::Sparse::Matrix LHS;
        INMOST::Sparse::Matrix LHS_low;

        INMOST::Sparse::Vector RHS;
        INMOST::Sparse::Vector RHS_low;
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