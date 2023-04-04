#pragma once

#include "mesh.hpp"
#include "dynvar.hpp"
#include "vecmath.hpp"
#include "local_assembling.hpp"
#include "tools.hpp"
#include "constants.hpp"

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
    
    // variables
    protected:

        // gradients of basis functions (in trian basis)
        INMOST::Tag grad_basis_func_tags;

        // lumped mass matrix tag
        INMOST::Tag lumped_mass_matrix_tag;

    private:

        // compute gradients of basis functions in trian basis
        void ComputeGradientsBasisFuncs(INMOST::Tag grad_bas_tags);

        // compute lumped mass matrix
        void ComputeMassMatrix(INMOST::Tag mm_tag);
        
    };

    class CgridMomentumSolver: public MomentumSolver
    {
        // constructor
    public:
        CgridMomentumSolver(SIMUG::IceMesh* mesh_,
                            double time_step_,
                            velocity_tag vel_tag_,
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
    
    // variables
    protected:

        //  solver parameters
        std::vector<double> real_params;
        std::vector<int> integer_params;
    
        // opposite edge num for every node on every triangle
        INMOST::Tag opposite_edge_for_node_tags;

        // opposite node num for every edge on every triangle
        INMOST::Tag opposite_node_for_edge_tags;

        // gradients of basis functions (in trian basis)
        INMOST::Tag grad_basis_func_tags; 

        // mass matrix tag
        INMOST::Tag mass_matrix_tag;
    
    // private functions
    private:
        // computation of opposite node number for every edge of triangle
        void ComputeOppositeEdges(INMOST::Tag op_edge_tags);

        // computation of opposite node number for every edge of triangle
        void ComputeOppositeNodes(INMOST::Tag op_node_tags);

        // compute gradients of basis functions in trian basis
        void ComputeGradientsBasisFuncs(INMOST::Tag grad_bas_tags);

        // compute lumped mass matrix
        void ComputeMassMatrix(INMOST::Tag mm_tag);

    };


    class Agrid_mEVP_Solver : public AgridMomentumSolver
    {
    // constructor
    public:
        Agrid_mEVP_Solver(SIMUG::IceMesh* mesh_,
                          double time_step_,
                          velocity_tag vel_tag_,
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

        void UpdateScalars();
        void ComputeP();
        void ComputeVarepsilonDelta();
        void AssembleForceVector();
        void UpdateSigmaMevp();
        void ComputeLevelVector();
        void UpdateVelocityMevp();
        void UpdateVelocity();
        void UpdateSigma();
        void MoveVectors(bool to_local);
        void Initialize();
        void Finalize();
        void LogError(int presudostep);
        void ComputeShearDeformation();
        double GetMaxVelocity();
    
    private:
        INMOST::Tag ua_tags;
        INMOST::Tag uw_tags;

        INMOST::Tag old_vel_tag;
        INMOST::Tag prev_vel_tag;
        INMOST::Tag new_vel_tag;

        INMOST::Tag prev_sigma_tag;
        INMOST::Tag new_sigma_tag;
        INMOST::Tag old_sigma_tag;

        INMOST::Tag delta_tag;
        INMOST::Tag P_tag;
        INMOST::Tag vareps_tag;
        INMOST::Tag force_tags;
        INMOST::Tag disc_level_tags;
        INMOST::Tag level_tag;
        
        INMOST::Tag shear_deformation_tag;
        INMOST::Tag z_component_tag;

        // profiling
    private:
        double strain_rate_computation_time = 0.0;
        double force_assembling_time = 0.0;
        double velocity_computation_time = 0.0;
    };


    class Cgrid_mEVP_Solver : public CgridMomentumSolver
    {
    // constructor
    public:
        Cgrid_mEVP_Solver(SIMUG::IceMesh* mesh_,
                          double time_step_,
                          velocity_tag vel_tag_,
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
        
        void LogError(int presudostep);
        void Initialize();
        void Finalize();
        void UpdateSigma();
        void UpdateVelocity();
        void UpdateScalars();
        void MoveVectors(bool to_local);
        void ComputeP();
        void ComputeVarepsilonDelta();
        void AssembleForceVector();
        void UpdateSigmaMevp();
        void UpdateVelocityMevp();
        void ComputeEdgeStabilizationSum();
        void ComputeEdgeStabilization(double alpha);
        void ComputeLevelVector();
        void ComputeShearDeformation();
        double GetMaxVelocity();
        
    
    private:

        INMOST::Tag edge_thick_tag;
        INMOST::Tag edge_conc_tag;

        INMOST::Tag ua_tags;
        INMOST::Tag uw_tags;

        INMOST::Tag old_vel_tag;
        INMOST::Tag prev_vel_tag;
        INMOST::Tag new_vel_tag;

        INMOST::Tag old_sigma_tag;
        INMOST::Tag prev_sigma_tag;
        INMOST::Tag new_sigma_tag;

        INMOST::Tag delta_tag;
        INMOST::Tag P_tag;
        INMOST::Tag vareps_tag;
        INMOST::Tag xi_edge_tag;
        INMOST::Tag force_tags;
        INMOST::Tag edge_stab_tags;
        INMOST::Tag edge_stab_sum_tags;
        INMOST::Tag disc_level_tags;
        INMOST::Tag level_tag;

        INMOST::Tag shear_deformation_tag;
        INMOST::Tag check_vel_tag;
        INMOST::Tag z_component_tag;
        INMOST::Tag check_ghost_edge_tag;
    
    // profiling
    private:
        double strain_rate_computation_time = 0.0;
        double force_assembling_time = 0.0;
        double velocity_computation_time = 0.0;
        double stabilization_assembling_time = 0.0;
    };
}