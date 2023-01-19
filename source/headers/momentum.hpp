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

    class CgridMomentumSolver: public MomentumSolver
    {
        // constructor
    public:
        CgridMomentumSolver(SIMUG::IceMesh* mesh_,
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
    
        // edge based mass matrix entries
        INMOST::Tag mass_matrix_entry_tag;
    
        // opposite edge num for every node on every triangle
        INMOST::Tag opposite_edge_for_node_tags;

        // opposite node num for every edge on every triangle
        INMOST::Tag opposite_node_for_edge_tags

        // edge basis in trian coords tags
        INMOST::Tag outward_edge_basis_in_trian_coords_tags;

        // height of triangle to edge
        INMOST::Tag trian_height_to_edge_tags;

        // gradients of basis functions (in trian basis)
        INMOST::Tag grad_basis_func_tags; 
    
    // private functions
    private:
        // procedure of assembling edge-based mass matrix
        void AssembleMassMatrix(INMOST::Tag mass_matrix_tag);

        // computation of opposite node number for every edge of triangle
        void ComputeOppositeEdges(INMOST::Tag op_edge_tags);

        // computation of opposite node number for every edge of triangle
        void ComputeOppositeNodes(INMOST::Tag op_node_tags);

        // compute local edge basis in trian coords
        void ComputeEdgeBasisInTrianCoords(INMOST::Tag edge_basis_in_trian_coords_tags);

        // compute trian heights
        void ComputeTrianHeights(INMOST::Tag trian_hight_tags);

        // compute gradients of basis functions in trian basis
        void ComputeGradientsBasisFuncs(INMOST::Tag grad_bas_tags);

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
    
    private:
    };


    class Cgrid_mEVP_Solver : public CgridMomentumSolver
    {
    // constructor
    public:
        Cgrid_mEVP_Solver(SIMUG::IceMesh* mesh_,
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
        
        void ComputeP();
        void ComputeVarepsilonDelta(INMOST::Tag vel_tag);
        void AssembleForceVector(INMOST::Tag sig_tag);

    
    protected:
        INMOST::Tag mass_matrix_entry;
    
    private:
        INMOST::Tag prev_vel_tag;
        INMOST::Tag new_vel_tag;
        INMOST::Tag sigma_tag;
        INMOST::Tag delta_tag;
        INMOST::Tag P_tag;
        INMOST::Tag vareps_tag;
        INMOST::Tag force_tags;
    };

    
}