#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <filesystem>
#include "inmost.h"

#include "defines.hpp"
#include "data.hpp"
#include "timer.hpp"
#include "logger.hpp"
#include "mesh_info.hpp"

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
    // id interval for set of elements on processor   
    struct IdInterval
    {
        int id_min;
        int id_max;
    };

    // general mesh information on processor
    struct MeshInfo
    {
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

    class IceMesh
    {
        typedef std::map<int, std::unique_ptr<GridData>> LayersDataMap;

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
        inline INMOST::Mesh* GetMesh()
        {return ice_mesh.get();};

        // Get map <layer number -> pointer to GridData> (could be (bnd)node, (bnd)edge, (bnd)trian)
        inline LayersDataMap& GetData(const gridElemType& gdtype)
        {return grid_data[gdtype];};

        // const version of GetData()
        inline const LayersDataMap& GetData(const gridElemType& gdtype) const
        {return grid_data.at(gdtype);};

        // Get mesh information (number of elements and local processor ids) 
        inline MeshInfo& GetInfo()
        {return mesh_info;}; 

        // const version of GetInfo()
        inline const MeshInfo& GetInfo() const
        {return mesh_info;};

        // get number of ice layers (always const)
        inline int GetNumLayers() const
        {return num_ice_layers;};
    
        // save mesh with data to a file
        void SaveVTU(const std::string& filename) const;

        // Destructor that frees up the grid data
        //~IceMesh();   

    private:

        void Partition();
        void SetBoundaryElements(bool display);
        void AssignVariables();
        void AssignCoords();
        void AssignIntervals();
        void ComputeMeshInfo();
        void AssembleBasisData(bool display);

    private:
        std::shared_ptr<INMOST::Mesh> ice_mesh;
        std::map<gridElemType, LayersDataMap> grid_data;
        MeshInfo mesh_info;
        int num_ice_layers;
    };
}