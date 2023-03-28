#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <memory>
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
#include "coords_rotation.hpp"

#ifdef INTEL_COMPILER
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

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
        std::map<mesh::meshVar, mesh::gridElemType> multi_elems;
        std::map<mesh::meshVar, mesh::gridElemType> single_elems;

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
        
        INMOST::Tag trans_matr_from_geo_to_elem;                // tranfile_infsition 2x2 matrix for vector/tensor coordinates while switching from geo to element basis
        INMOST::Tag trans_matr_from_elem_to_geo;                // transition 2x2 matrix for vector/tensor coordinates while switching from element to geo basis ( = trans_matr_from_geo_to_elem^-1)

        std::vector<double> VectorFromGeoToElem(const std::vector<double>& vec_geo_coords);
        std::vector<double> VectorFromElemToGeo(const std::vector<double>& vec_elem_coords);

        std::vector<double> TensorFromGeoToElem(const std::vector<double>& tens_geo_coords);
        std::vector<double> TensorFromElemToGeo(const std::vector<double>& tens_elem_coords);

        virtual std::vector<INMOST::Tag>&  GetTransMatrToNode() {SIMUG_ERR("Current elent type doesn't have transition matrix to node!");};
        virtual std::vector<INMOST::Tag>&  GetTransMatrToEdge() {SIMUG_ERR("Current elent type doesn't have transition matrix to edge!");};
        virtual std::vector<INMOST::Tag>&  GetTransMatrToTrian(){SIMUG_ERR("Current elent type doesn't have transition matrix to triangle!");};
        virtual std::vector<INMOST::Tag>&  GetNodeCoordsInTrianBasis() {SIMUG_ERR("Current elent type doesn't have node coords in triangular basis!");};
        virtual std::vector<INMOST::Tag>& GetIsXedgeBasisIsNormal() {SIMUG_ERR("Current elent type doesn't contain information about normal!");};

    public:

        virtual void Mute() = 0;
        virtual void UnMute() = 0;
        virtual INMOST::Tag& GetCartesianSize() = 0;

    protected:
        INMOST::Mesh* ice_mesh;

    };

    class NodeInfo : public GridInfo
    {
    public:
        NodeInfo(INMOST::Mesh* ice_mesh_): GridInfo(ice_mesh_) {};
        std::vector<INMOST::Tag>&  GetTransMatrToEdge() {return trans_matr_to_edge;};
        std::vector<INMOST::Tag>&  GetTransMatrToTrian() {return trans_matr_to_trian;};
        INMOST::Tag&  GetCartesianSize() {return node_cart_size;};

    public:

        void Mute();
        void UnMute();
    
    private:
        std::vector<INMOST::Tag> trans_matr_to_edge;
        std::vector<INMOST::Tag> trans_matr_to_trian;
        INMOST::Tag node_cart_size;
    };

    class EdgeInfo : public GridInfo
    {
    public:
        EdgeInfo(INMOST::Mesh* ice_mesh_): GridInfo(ice_mesh_) {};
        std::vector<INMOST::Tag>&  GetTransMatrToNode() {return trans_matr_to_node;};
        std::vector<INMOST::Tag>&  GetTransMatrToTrian() {return trans_matr_to_trian;};
        INMOST::Tag&  GetCartesianSize() {return edge_cart_size;};

    public:
        void Mute();
        void UnMute();
    
    private:
        std::vector<INMOST::Tag> trans_matr_to_node;
        std::vector<INMOST::Tag> trans_matr_to_trian;
        INMOST::Tag edge_cart_size;
    };

    class TrianInfo : public GridInfo
    {
    public:
        TrianInfo(INMOST::Mesh* ice_mesh_): GridInfo(ice_mesh_) {};
        std::vector<INMOST::Tag>&  GetTransMatrToNode() {return trans_matr_to_node;};
        std::vector<INMOST::Tag>&  GetTransMatrToEdge() {return trans_matr_to_edge;};
        std::vector<INMOST::Tag>&  GetNodeCoordsInTrianBasis() {return node_coords_in_trian_basis;};
        std::vector<INMOST::Tag>& GetIsXedgeBasisIsNormal() {return is_x_edge_basis_normal_vec;};
        INMOST::Tag&  GetCartesianSize() {return trian_cart_size;};

    public:
        void Mute();
        void UnMute();
    private:
        std::vector<INMOST::Tag> trans_matr_to_node;
        std::vector<INMOST::Tag> trans_matr_to_edge;

        std::vector<INMOST::Tag> node_coords_in_trian_basis;
        std::vector<INMOST::Tag> is_x_edge_basis_normal_vec;
        INMOST::Tag trian_cart_size;
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
        inline std::shared_ptr<GridData>& GetDataMulti(const mesh::gridElemType& gdtype, int layer) {return multi_data[gdtype][layer];};
        inline const std::shared_ptr<GridData>& GetDataMulti(const mesh::gridElemType& gdtype, int layer) const {return multi_data.at(gdtype).at(layer);};
        inline void MuteDataMulti(int layer) {for (auto& [var, elem]: mesh_info.multi_elems){multi_data[elem][layer]->Mute(var);};};
        inline void MuteDataMulti() {for (int layer = 1; layer < mesh_info.num_ice_layers; ++layer){MuteDataMulti(layer);};};

        // Get and Mute forcing data
        inline std::shared_ptr<GridData>& GetDataSingle(const mesh::gridElemType& gdtype) {return single_data[gdtype];};
        inline const std::shared_ptr<GridData>& GetDataSingle(const mesh::gridElemType& gdtype) const {return single_data.at(gdtype);};
        void MuteDataSingle() {for (auto& [var, elem]: mesh_info.single_elems){single_data[elem]->Mute(var);};};

        // Get mesh information (number of elements and local processor ids) 
        inline MeshInfo& GetMeshInfo() {return mesh_info;}; 
        inline const MeshInfo& GetMeshInfo() const {return mesh_info;};

        // Get mesh information (number of elements and local processor ids) 
        inline std::shared_ptr<GridInfo>& GetGridInfo(const mesh::gridElemType& gdtype) {return grid_info[gdtype];}; 
        inline const std::shared_ptr<GridInfo>& GetGridInfo(const mesh::gridElemType& gdtype) const {return grid_info.at(gdtype);};
    
        // save mesh with data to a file
        void SaveVTU(const std::string& meshname) const;
        void SaveVTU(const std::string& meshname, int postscript) const;

        // get vector of bnd nodes on current processor
        inline INMOST::ElementArray<INMOST::Node>& GetBndNodes() {return bnd_nodes;};
        const INMOST::ElementArray<INMOST::Node>& GetBndNodes() const {return bnd_nodes;};

        // get array of bnd edges on current processor
        INMOST::ElementArray<INMOST::Face>& GetBndEdges() {return bnd_edges;};
        const INMOST::ElementArray<INMOST::Face>& GetBndEdges() const {return bnd_edges;};

        // get array of bnd triangles on current processor
        INMOST::ElementArray<INMOST::Cell>& GetBndTrians() {return bnd_trians;};
        const INMOST::ElementArray<INMOST::Cell>& GetBndTrians() const {return bnd_trians;};

        // vector transition functions between element bases
        std::vector<double> VecTransition(const std::vector<double>& vec_coords, const INMOST::Node& node, const INMOST::Face& edge);
        std::vector<double> VecTransition(const std::vector<double>& vec_coords, const INMOST::Node& node, const INMOST::Cell& trian);
        std::vector<double> VecTransition(const std::vector<double>& vec_coords, const INMOST::Face& edge, const INMOST::Node& node);
        std::vector<double> VecTransition(const std::vector<double>& vec_coords, const INMOST::Face& edge, const INMOST::Cell& trian);
        std::vector<double> VecTransition(const std::vector<double>& vec_coords, const INMOST::Cell& trian, const INMOST::Node& node);
        std::vector<double> VecTransition(const std::vector<double>& vec_coords, const INMOST::Cell& trian, const INMOST::Face& edge);

        // tensor transition functions between element bases
        std::vector<std::vector<double>> TensTransition(const std::vector<std::vector<double>>& tens_coords, const INMOST::Node& node, const INMOST::Face& edge);
        std::vector<std::vector<double>> TensTransition(const std::vector<std::vector<double>>& tens_coords, const INMOST::Node& node, const INMOST::Cell& trian);
        std::vector<std::vector<double>> TensTransition(const std::vector<std::vector<double>>& tens_coords, const INMOST::Face& edge, const INMOST::Node& node);
        std::vector<std::vector<double>> TensTransition(const std::vector<std::vector<double>>& tens_coords, const INMOST::Face& edge, const INMOST::Cell& trian);
        std::vector<std::vector<double>> TensTransition(const std::vector<std::vector<double>>& tens_coords, const INMOST::Cell& trian, const INMOST::Node& node);
        std::vector<std::vector<double>> TensTransition(const std::vector<std::vector<double>>& tens_coords, const INMOST::Cell& trian, const INMOST::Face& edge);

        // vector transition functions between geo and element basis
        std::vector<double> VecTransitionToElemBasis(const std::vector<double>& vec_coords, const INMOST::Node& node);
        std::vector<double> VecTransitionToElemBasis(const std::vector<double>& vec_coords, const INMOST::Face& edge);
        std::vector<double> VecTransitionToElemBasis(const std::vector<double>& vec_coords, const INMOST::Cell& trian);

        std::vector<double> VecTransitionToGeoBasis(const std::vector<double>& vec_coords, const INMOST::Node& node);
        std::vector<double> VecTransitionToGeoBasis(const std::vector<double>& vec_coords, const INMOST::Face& edge);
        std::vector<double> VecTransitionToGeoBasis(const std::vector<double>& vec_coords, const INMOST::Cell& trian);

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
        void ComputeNodeCoordsInTrianBasis();
        void ComputeIfXedgeBasisIsNormalToTrian();
        void ComputeElementsCartesianSize();

    private:

        std::shared_ptr<INMOST::Mesh> ice_mesh;
        std::map<mesh::gridElemType, LayersDataMap> multi_data;
        std::map<mesh::gridElemType, std::shared_ptr<GridData>> single_data;
        std::map<mesh::gridElemType, std::shared_ptr<GridInfo>> grid_info;
        INMOST::ElementArray<INMOST::Node> bnd_nodes;
        INMOST::ElementArray<INMOST::Face> bnd_edges;
        INMOST::ElementArray<INMOST::Cell> bnd_trians;
        MeshInfo mesh_info;
    };
}