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

#ifdef USE_MPI
#include "mpi.h"
#endif

//#include "modvar.hpp"
//#include "coordvar.hpp"
//#include "basis.hpp"
//#include "model_config.hpp"
//#include "mesh_config.hpp"

namespace SIMUG::mesh
{
    // general mesh information on processor
    struct MeshInfo
    {   
        struct IdInterval
        {
            int id_min;
            int id_max;
        };

        surfType surface_type;
        gridType grid_type;

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
        INMOST::Tag id;
        INMOST::Tag id_no_bnd;
        INMOST::Tag is_bnd;
        std::map<SIMUG::coord::coordType, INMOST::Tag> coords;

    public:
        virtual void Exchange() = 0;
        virtual void Mute() = 0;
        virtual void UnMute() = 0;

    protected:
        INMOST::Mesh* ice_mesh;
        void ExchangeAll(const unsigned char& GridElem);

    };

    class NodeInfo : public GridInfo
    {
    public:
        NodeInfo(INMOST::Mesh* ice_mesh_): GridInfo(ice_mesh_) {};
    
    public:
        void Exchange();
        void Mute();
        void UnMute();
    };

    class EdgeInfo : public GridInfo
    {
    public:
        EdgeInfo(INMOST::Mesh* ice_mesh_): GridInfo(ice_mesh_) {};
    
    public:
        void Exchange();
        void Mute();
        void UnMute();

    };

    class TrianInfo : public GridInfo
    {
    public:
        TrianInfo(INMOST::Mesh* ice_mesh_): GridInfo(ice_mesh_) {};
    
    public:
        void Exchange();
        void Mute();
        void UnMute();
    };

    class IceMesh
    {
        typedef std::map<int, std::shared_ptr<GridData>> LayersDataMap;

    public:

        // construct manually with number of ice layers
        IceMesh(const std::string& path_to_file_,
                const surfType& surf_type_,
                const gridType& grid_type_,
                const int& n_ice_layers_);

        // construct manually with 1 ice layer
        IceMesh(const std::string& path_to_file_,
                const surfType& surf_type_,
                const gridType& grid_type_);

        // construct with mesh config class and number of ice layers (to do)
        //IceMesh(const MeshConfig& mesh_config_);

        // Get INMOST::Mesh pointer
        inline INMOST::Mesh* GetMesh() {return ice_mesh.get();};

        // Get map <layer number -> pointer to GridData> (could be (bnd)node, (bnd)edge, (bnd)trian)
        inline std::shared_ptr<GridData>& GetData(const gridElemType& gdtype, int layer) {return grid_data[gdtype][layer];};
        inline const std::shared_ptr<GridData>& GetData(const gridElemType& gdtype, int layer) const {return grid_data.at(gdtype).at(layer);};

        // Get mesh information (number of elements and local processor ids) 
        inline MeshInfo& GetMeshInfo() {return mesh_info;}; 
        inline const MeshInfo& GetMeshInfo() const {return mesh_info;};

        // Get mesh information (number of elements and local processor ids) 
        inline std::map<gridElemType, std::shared_ptr<GridInfo>>& GetGridInfo() {return grid_info;}; 
        inline const std::map<gridElemType, std::shared_ptr<GridInfo>>& GetGridInfo() const {return grid_info;};
        
        //Get number of ice layers 
        inline const int GetNumLayers() const
        {return num_ice_layers;};
    
        // save mesh with data to a file
        void SaveVTU(const std::string& filename) const;

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
        void AssignGridVariables();
        void AssignCoords();
        void AssignIds();
        void AssignIdIntervals();
        void AssembleBasisData();

    private:

        std::shared_ptr<INMOST::Mesh> ice_mesh;
        std::map<gridElemType, LayersDataMap> grid_data;
        std::map<gridElemType, std::shared_ptr<GridInfo>> grid_info;
        INMOST::ElementArray<INMOST::Node> bnd_nodes;
        INMOST::ElementArray<INMOST::Face> bnd_edges;
        INMOST::ElementArray<INMOST::Cell> bnd_trians;
        MeshInfo mesh_info;
        int num_ice_layers;
    };
}