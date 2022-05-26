#pragma once

#include "mesh.hpp"

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

        // add scalar for advection 
        inline void AddScalar(scalar_tag scal_tag_)
        {scal_tags.push_back(scal_tag);};

        // add vector of scalars for advection 
        void AddScalar(const std::vector<scalar_tag>& scal_tags_)
        {for (auto& sc_tag: scal_tags_){scal_tags.push_back(sc_tag);}};

    protected:
        // transport one scalar
        virtual void Evaluate(scalar_tag scal_tag) = 0;

    public:
        // transport all scalars
        inline void TransportScalars()
        {
            for (auto& scal_tag: scal_tags)
            {    
                Evaluate(vel_tag, scal_tag);
            } 
            BARRIER
        };

    protected:
        SIMUG::IceMesh* mesh;
        double time_step;
        velocity_tag vel_tag;
        std::vector<scalar_tag> scal_tags;
    };

    class AgridAdvectionSolver: public AdvectionSolver
    {
    public:
        AgridAdvectionSolver(SIMUG::IceMesh* mesh_,
                             double time_step_,
                             bool is_fct,
                             double fct_cd);
    private:
        void Evaluate(scalar_tag scal_tag) override;

        double ComputeMaxCourant();

        void AssembleLHS();
        void AssembleRHS();
        void Evaluate();

    private:
        INMOST::Solver* slae_solver;
        adv::scheme adv_scheme;
        adv::filter adv_filter;
        std::vector<double> params;

    private:
        INMOST::Sparse::Matrix LHS;
        std::optional<INMOST::Sparse::Matrix> LHS_low;

        INMOST::Sparse::Vector RHS;
        std::optional<INMOST::Sparse::Vector> RHS_low;



    private:
        void AssembleSingleStepRHS();
        void AssembleDoubleStepRHS(StepNumber step_num);
        void LogInit();
        void LogStep();

        IceMesh& ice_mesh;
        AdvectionParams& advection_params;
        ModelParams& model_params;
        MeshParams& mesh_params;
        ModelVariableNotation transported_scalar;
        ModelVariableNotation transporting_velocity;
        INMOST::Tag m_tag;
        INMOST::Tag u_tag;
        INMOST::Tag m_high_tag;
        INMOST::Tag m_low_tag;
        INMOST::Tag m_half_tag;
        INMOST::Sparse::Matrix LHS;
        INMOST::Sparse::Matrix LHS_low;
        INMOST::Sparse::Vector RHS;
        INMOST::Sparse::Vector RHS_low;
        INMOST::Solver& Sol;
        std::vector<std::vector<std::vector<double>>> M_C_minus_M_L;
        bool is_verbose_advection = true;
        double init_mass_integral;
    };

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
}