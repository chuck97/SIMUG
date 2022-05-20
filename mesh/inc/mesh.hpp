#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <filesystem>
#include <cmath>
#include "inmost.h"

#include "defines.hpp"
#include "data.hpp"
#include "timer.hpp"
#include "logger.hpp"
#include "mesh_info.hpp"
#include "stringtrim.hpp"
#include "coordvar.hpp"
#include "gridvars.hpp"
#include "vecmath.hpp"

#ifdef USE_MPI
#include "mpi.h"
#endif

//#include "modvar.hpp"
//#include "coordvar.hpp"
//#include "basis.hpp"
//#include "model_config.hpp"
//#include "mesh_config.hpp"

namespace SIMUG
{
    // general mesh information on processor
    struct MeshInfo
    {   
        struct IdInterval
        {
            int id_min;
            int id_max;
        };

        mesh::surfType surface_type;
        mesh::gridType grid_type;
        std::string output_folder;
        int num_ice_layers;
        std::map<mesh::meshVar, mesh::gridElemType> prog_elems;
        std::map<mesh::meshVar, mesh::gridElemType> forc_elems;

        int num_nodes;
        int num_edges;
        int num_trians;

        int num_bnd_nodes;
        int num_bnd_edges;
        int num_bnd_trians;

        IdInterval id_interval_nodes;
        IdInterval id_interval_edges;
        IdInterval id_interval_trians;

        IdInterval id_interval_nodes_no_bnd;
        IdInterval id_interval_edges_no_bnd;
        IdInterval id_interval_trians_no_bnd;
    };

    // grid information for particular element
    class GridInfo
    {
    public:
        GridInfo(INMOST::Mesh* ice_mesh_): ice_mesh(ice_mesh_) {}; 

    public:
        INMOST::Tag id;                                         // global element id 
        INMOST::Tag id_no_bnd;                                  // global element id without bnd elements
        INMOST::Tag is_bnd;                                     // is node bnd (1 = true, 0 = false)
        
        std::map<coord::coordType, INMOST::Tag> coords;         // coordinates of element centroid (model, geographical, cartesian)
        
        std::vector<INMOST::Tag> geo_basis;                     // geographical basis vectors coordinates
        std::vector<INMOST::Tag> cart_basis;                    // local cartesian basis vectors coordinates
        
        INMOST::Tag trans_matr_from_geo_to_elem;                // transition 2x2 matrix for vector/tensor coordinates while switching from geo to element basis
        INMOST::Tag trans_matr_from_elem_to_geo;                // transition 2x2 matrix for vector/tensor coordinates while switching from element to geo basis ( = trans_matr_from_geo_to_elem^-1)

        virtual std::vector<INMOST::Tag>&  GetTransMatrToNode() {SIMUG_ERR("Current elent type doesn't have transition matrix to node!");};
        virtual std::vector<INMOST::Tag>&  GetTransMatrToEdge() {SIMUG_ERR("Current elent type doesn't have transition matrix to edge!");};
        virtual std::vector<INMOST::Tag>&  GetTransMatrToTrian(){SIMUG_ERR("Current elent type doesn't have transition matrix to triangle!");};

    public:

        virtual void Mute() = 0;
        virtual void UnMute() = 0;

    protected:
        INMOST::Mesh* ice_mesh;

    };

    class NodeInfo : public GridInfo
    {
    public:
        NodeInfo(INMOST::Mesh* ice_mesh_): GridInfo(ice_mesh_) {};
        std::vector<INMOST::Tag>&  GetTransMatrToEdge() {return trans_matr_to_edge;};
        std::vector<INMOST::Tag>&  GetTransMatrToTrian() {return trans_matr_to_trian;};

    public:

        void Mute();
        void UnMute();
    
    private:
        std::vector<INMOST::Tag> trans_matr_to_edge;
        std::vector<INMOST::Tag> trans_matr_to_trian;
    };

    class EdgeInfo : public GridInfo
    {
    public:
        EdgeInfo(INMOST::Mesh* ice_mesh_): GridInfo(ice_mesh_) {};
        std::vector<INMOST::Tag>&  GetTransMatrToNode() {return trans_matr_to_node;};
        std::vector<INMOST::Tag>&  GetTransMatrToTrian() {return trans_matr_to_trian;};

    public:
        void Mute();
        void UnMute();
    
    private:
        std::vector<INMOST::Tag> trans_matr_to_node;
        std::vector<INMOST::Tag> trans_matr_to_trian;
    };

    class TrianInfo : public GridInfo
    {
    public:
        TrianInfo(INMOST::Mesh* ice_mesh_): GridInfo(ice_mesh_) {};
        std::vector<INMOST::Tag>&  GetTransMatrToNode() {return trans_matr_to_node;};
        std::vector<INMOST::Tag>&  GetTransMatrToEdge() {return trans_matr_to_edge;};

    public:
        void Mute();
        void UnMute();
    private:
        std::vector<INMOST::Tag> trans_matr_to_node;
        std::vector<INMOST::Tag> trans_matr_to_edge;
    };

    class IceMesh
    {
        typedef std::map<int, std::shared_ptr<GridData>> LayersDataMap;

    public:

        // construct manually with number of ice layers (with output folder)
        IceMesh(const std::string& path_to_file_,
                const std::string& output_folder_,
                const mesh::surfType& surf_type_,
                const mesh::gridType& grid_type_,
                const int& n_ice_layers_);
        
        // construct manually with number of ice layers (without output folder)
        IceMesh(const std::string& path_to_file_,
                const mesh::surfType& surf_type_,
                const mesh::gridType& grid_type_,
                const int& n_ice_layers_);

        // construct manually with 1 ice layer (with output folder)
        IceMesh(const std::string& path_to_file_,
                const std::string& output_folder_,
                const mesh::surfType& surf_type_,
                const mesh::gridType& grid_type_);
        
        // construct manually with 1 ice layer (without output folder)
        IceMesh(const std::string& path_to_file_,
                const mesh::surfType& surf_type_,
                const mesh::gridType& grid_type_);

            

        // construct with mesh config class and number of ice layers (to do)
        //IceMesh(const MeshConfig& mesh_config_);

        // Get INMOST::Mesh pointer
        inline INMOST::Mesh* GetMesh() {return ice_mesh.get();};

        // Get and Mute prognostic data
        inline std::shared_ptr<GridData>& GetProgData(const mesh::gridElemType& gdtype, int layer) {return prognostic_data[gdtype][layer];};
        inline const std::shared_ptr<GridData>& GetProgData(const mesh::gridElemType& gdtype, int layer) const {return prognostic_data.at(gdtype).at(layer);};
        inline void MuteProgData(int layer) {for (auto& [var, elem]: mesh_info.prog_elems){prognostic_data[elem][layer]->Mute(var);};};
        inline void MuteProgData() {for (int layer = 1; layer < mesh_info.num_ice_layers; ++layer){MuteProgData(layer);};};

        // Get and Mute forcing data
        inline std::shared_ptr<GridData>& GetForcData(const mesh::gridElemType& gdtype) {return forcing_data[gdtype];};
        inline const std::shared_ptr<GridData>& GetForcData(const mesh::gridElemType& gdtype) const {return forcing_data.at(gdtype);};
        void MuteForcData() {for (auto& [var, elem]: mesh_info.forc_elems){forcing_data[elem]->Mute(var);};};

        // Get mesh information (number of elements and local processor ids) 
        inline MeshInfo& GetMeshInfo() {return mesh_info;}; 
        inline const MeshInfo& GetMeshInfo() const {return mesh_info;};

        // Get mesh information (number of elements and local processor ids) 
        inline std::shared_ptr<GridInfo>& GetGridInfo(const mesh::gridElemType& gdtype) {return grid_info[gdtype];}; 
        inline const std::shared_ptr<GridInfo>& GetGridInfo(const mesh::gridElemType& gdtype) const {return grid_info.at(gdtype);};
    
        // save mesh with data to a file
        void SaveVTU(const std::string& meshname) const;

        // get vector of bnd nodes on current processor
        inline INMOST::ElementArray<INMOST::Node>& GetBndNodes() {return bnd_nodes;};
        const INMOST::ElementArray<INMOST::Node>& GetBndNodes() const {return bnd_nodes;};

        // get array of bnd edges on current processor
        INMOST::ElementArray<INMOST::Face>& GetBndEdges() {return bnd_edges;};
        const INMOST::ElementArray<INMOST::Face>& GetBndEdges() const {return bnd_edges;};

        // get array of bnd triangles on current processor
        INMOST::ElementArray<INMOST::Cell>& GetBndTrians() {return bnd_trians;};
        const INMOST::ElementArray<INMOST::Cell>& GetBndTrians() const {return bnd_trians;};

        // Destructor that frees up the grid data
        //~IceMesh();   

    private:

        void Partition();
        void ComputeMeshInfo();
        void SelectBndNodes();
        void SelectBndEdges();
        void SelectBndEdgesNoGhost();
        void SelectBndTrians();
        void AssignPrognosticVariables();
        void AssignForcingVariables();
        void AssignCoords();
        void AssignIds();
        void AssignIdIntervals();

        void AssembleBasisData();
        void AssembleGeoElementBasis();
        void AssembleCartesianElementBasis();
        void AssembleGeoToElementTransitionMatricies();
        void AssembleElementToElementTransitionMatricies();

    private:

        std::shared_ptr<INMOST::Mesh> ice_mesh;
        std::map<mesh::gridElemType, LayersDataMap> prognostic_data;
        std::map<mesh::gridElemType, std::shared_ptr<GridData>> forcing_data;
        std::map<mesh::gridElemType, std::shared_ptr<GridInfo>> grid_info;
        INMOST::ElementArray<INMOST::Node> bnd_nodes;
        INMOST::ElementArray<INMOST::Face> bnd_edges;
        INMOST::ElementArray<INMOST::Cell> bnd_trians;
        MeshInfo mesh_info;
    };
}