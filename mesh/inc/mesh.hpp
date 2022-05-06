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
        // id interval for set of elements on processor   
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
    struct GridInfo
    {
        INMOST::Tag id;
        INMOST::Tag id_no_bnd;
        INMOST::Tag is_bnd;
        std::map<SIMUG::coord::coordType, INMOST::Tag> coords;

        virtual void Exchange(INMOST::Mesh* ice_mesh) = 0;
        //virtual void Mute(INMOST::Mesh* ice_mesh) = 0;
        //virtual void UnMute(INMOST::Mesh* ice_mesh) = 0;

        void ExchangeAll(INMOST::Mesh* ice_mesh, const unsigned char& GridElem)
        {
            ice_mesh->ExchangeData(id, GridElem, 0);
            ice_mesh->ExchangeData(id_no_bnd, GridElem, 0);
            ice_mesh->ExchangeData(is_bnd, GridElem, 0);

            for (auto [key, value]: coords)
                ice_mesh->ExchangeData(value, GridElem, 0);
            BARRIER
        };

        void MuteAll(INMOST::Mesh* ice_mesh)
        {

        };

        void UnMuteAll(INMOST::Mesh* ice_mesh)
        {

        };
    };

    struct NodeInfo : public GridInfo
    {
        inline void Exchange(INMOST::Mesh* ice_mesh)
        {ExchangeAll(ice_mesh, INMOST::NODE);};

        //inline void Mute(INMOST::Mesh* ice_mesh)
        //{};

        // additional local basis data
    };

    struct EdgeInfo : public GridInfo
    {
        inline void Exchange(INMOST::Mesh* ice_mesh)
        {ExchangeAll(ice_mesh, INMOST::FACE);};

        // additional local basis data
    };

    struct TrianInfo : public GridInfo
    {
        inline void Exchange(INMOST::Mesh* ice_mesh)
        {ExchangeAll(ice_mesh, INMOST::CELL);};

        // additional local basis data
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
        inline LayersDataMap& GetData(const gridElemType& gdtype) {return grid_data[gdtype];};
        inline const LayersDataMap& GetData(const gridElemType& gdtype) const {return grid_data.at(gdtype);};

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
        void SelectBndTrians();
        void AssignVariables();
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