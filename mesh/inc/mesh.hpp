#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include "inmost.h"

#include "defines.hpp"
#include "modvar.hpp"
#include "coordvar.hpp"
#include "basis.hpp"
#include "model_config.hpp"
#include "mesh_config.hpp"

namespace SIMUG::mesh
{
    struct IdInterval
    {
        int id_min;
        int id_max;
    };


    struct MeshInfo
    {
        int num_nodes;
        int num_bnd_edges;
        int num_elems;

        int num_nodes_bnd;

        int num_nodes_not_bnd;

        IdInterval id_interval_nodes;
        IdInterval id_interval_edges;
        IdInterval id_interval_elems;

        IdInterval id_interval_nodes_no_bnd;
    };

    template<typename RT>
    struct MeshData
    {
        std::map<var, INMOST::Tag> data;
        int id;
    };

    template<typename RT>
    struct NodeData : public MeshData<RT>
    {
        std::map<NodeCoordsNotation, INMOST::Tag> coords;
        int id_no_bnd;
    };

    template<typename RT>
    struct EdgeData : public MeshData<RT>
    {
        RT size;
    };

    template<typename RT>
    struct ElemData : public MeshData<RT>
    {
        RT size;
        RT area;
    };

    template<typename RT>
    class IceMesh
    {
    public:
        IceMesh(const MeshConfig& mesh_config_,
                const ModelConfig& model_config_,
                std::ostream& os_);  
        ~IceMesh();                                                            
        INMOST::Mesh* GetMesh();
        NodeData& GetNodeData();
        EdgeData& GetEdgeData();
        ElemData& GetElemData();
        MeshInfo& GetMeshInfo();

        ElemBasisData& GetBasisData();

        void SavePVTU(const std::string& filename) const;    

    private:

        void Partition();
        void AssignVariables();
        void AssignCoords();
        void SetBoundaryNodes(bool display);
        void AssignIntervals();
        void ComputeMeshInfo();
        void AssembleBasisData(bool display);


        INMOST::Mesh* ice_mesh;
        
        MeshInfo mesh_info;
        
        MeshConfig mesh_config;
        ModelConfig model_config;
        
        NodeData<RT> node_data;
        EdgeData<RT> edge_data;
        ElemData<RT> elem_data;

        ElemBasisData<RT> basis_data;

        std::ostream os;
    };
}