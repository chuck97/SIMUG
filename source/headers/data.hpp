#pragma once

#include <map> 
#include <vector>
#include <string>
#include <memory>
#include <tuple>
#include <optional>
#include "inmost.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "mesh_data.hpp"
#include "defines.hpp"
#include "gridvars.hpp"


namespace SIMUG
{
    // Grid Data base class
    class GridData
    {
    public:
        typedef std::vector<double>(*FuncPtr)(std::pair<double, double>, double);

        typedef std::map<mesh::meshVar, std::pair<INMOST::Tag, std::optional<FuncPtr>>> ProgData;
        typedef std::map<std::string, std::pair<INMOST::Tag, std::optional<FuncPtr>>> TempData;

        typedef std::map<mesh::meshVar, size_t> ProgDataSize;
        typedef std::map<std::string, size_t> TempDataSize;

    protected:  
        INMOST::Mesh* ice_mesh;    

        ProgData prog_data;
        TempData temp_data;

        ProgDataSize prog_data_size;
        TempDataSize temp_data_size;

        std::optional<int> layer;

    public:
        inline GridData(INMOST::Mesh* ice_mesh_, int layer_):
            ice_mesh(ice_mesh_),
            layer(layer_)
        {};

        inline GridData(INMOST::Mesh* ice_mesh_):
            ice_mesh(ice_mesh_)
        {};

        // Create mesh data
        virtual void Create(const mesh::meshVar& pNot,  const mesh::meshDim& pDim, const INMOST::DataType& dtype) = 0;
        virtual void Create(const std::string& tVar, const mesh::meshDim& tDim, const INMOST::DataType& dtype) = 0;
        virtual void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) = 0;

        // Obtain mesh data
        inline INMOST::Tag& Get(const mesh::meshVar& pNot)
        {return prog_data.at(pNot).first;};

        inline INMOST::Tag& Get(const std::string& tName)
        {return temp_data.at(tName).first;};

        inline const INMOST::Tag& Get(const mesh::meshVar& pNot) const
        {return prog_data.at(pNot).first;};

        inline const INMOST::Tag& Get(const std::string& tName) const 
        {return temp_data.at(tName).first;};

        // Obtain analytical expressions
        std::optional<FuncPtr>& GetExpr(const mesh::meshVar& pNot)
        {return prog_data.at(pNot).second;};

        std::optional<FuncPtr>& GetExpr(const std::string& tName)
        {return temp_data.at(tName).second;};

        const std::optional<FuncPtr>& GetExpr(const mesh::meshVar& pNot) const 
        {return prog_data.at(pNot).second;};

        const std::optional<FuncPtr>& GetExpr(const std::string& tName) const
        {return temp_data.at(tName).second;};

        // Obtain size of mesh variable
        inline size_t GetSize(const mesh::meshVar& pNot) const
        {return prog_data_size.at(pNot);};

        inline size_t GetSize(const std::string& tName) const
        {return temp_data_size.at(tName);};


        // Exchange mesh data
        virtual void Exchange(const mesh::meshVar& pNot) = 0;
        virtual void Exchange(const std::string& tName) = 0; 

        // Mute data
        inline void Mute(const mesh::meshVar& pNot)
        {
            if (layer.has_value())
                ice_mesh->SetFileOption("Tag:" + mesh::meshVarName.at(pNot) + " " + std::to_string(layer.value()), "nosave");
            else
                ice_mesh->SetFileOption("Tag:" + mesh::meshVarName.at(pNot), "nosave");
        };

        inline void Mute(const std::string& tName)
        {
            if (layer.has_value())
                ice_mesh->SetFileOption("Tag:" + tName + " " + std::to_string(layer.value()), "nosave");
            else
                ice_mesh->SetFileOption("Tag:" + tName, "nosave");
        };

        // Unmute data
        inline void Unmute(const mesh::meshVar& pNot)
        {
            if (layer.has_value())
                ice_mesh->SetFileOption("Tag:" + mesh::meshVarName.at(pNot) + " " + std::to_string(layer.value()), "save");
            else
                ice_mesh->SetFileOption("Tag:" + mesh::meshVarName.at(pNot), "save");
        };
        
        inline void Unmute(const std::string& tName)
        {
            if (layer.has_value())
                ice_mesh->SetFileOption("Tag:" + tName + " " + std::to_string(layer.value()), "save");
            else
                ice_mesh->SetFileOption("Tag:" + tName, "save");
        };


        // Delete mesh data
        virtual void Delete(const mesh::meshVar& pNot) = 0;
        virtual void Delete(const std::string& tVar) = 0;
        virtual void Delete(const std::vector<mesh::meshVar>& pNots) = 0;
        virtual void Delete(const std::vector<std::string>& tVars) = 0;

    protected:
        void GridCreateData(const mesh::meshVar& pNot,  const mesh::meshDim& pDim, 
                            const INMOST::DataType& InmostDataType,
                            const unsigned char& GridElem,
                            const unsigned char& GridSparse);

        void GridCreateData(const std::string& tVar, const mesh::meshDim& tDim,
                            const INMOST::DataType& InmostDataType,
                            const unsigned char& GridElem,
                            const unsigned char& GridSparse);

        void GridCreateData(const std::string& tVar, const int& vSize,
                            const INMOST::DataType& InmostDataType,
                            const unsigned char& GridElem,
                            const unsigned char& GridSparse);

        inline void GridDeleteData(const mesh::meshVar& pNot, unsigned char InmostGridElem)
        {if (prog_data.count(pNot)==0) {SIMUG_ERR("variable "+mesh::meshVarName.at(pNot)+" doesn't exist")};
        ice_mesh->DeleteTag(prog_data.at(pNot).first, InmostGridElem);};

        inline void GridDeleteData(const std::string& tVar, unsigned char InmostGridElem)
        {if (temp_data.count(tVar)==0) {SIMUG_ERR("variable "+tVar+" doesn't exist")};
        ice_mesh->DeleteTag(temp_data.at(tVar).first, InmostGridElem);};

        inline void GridDeleteData(const std::vector<mesh::meshVar>& pNots, unsigned char InmostGridElem)
        {for (const auto& item: pNots){if (prog_data.count(item)==0) {SIMUG_ERR("variable "+mesh::meshVarName.at(item)+" doesn't exist")};
        ice_mesh->DeleteTag(prog_data.at(item).first, InmostGridElem);};};

        inline void GridDeleteData(const std::vector<std::string>& tVars, unsigned char InmostGridElem)
        {for (const auto& item: tVars){if (temp_data.count(item)==0) {SIMUG_ERR("variable "+item+" doesn't exist")};
        ice_mesh->DeleteTag(temp_data.at(item).first, InmostGridElem);};};
   
    };


    // Node data inharited class
    class NodeData: public GridData
    {
    public:
        inline NodeData(INMOST::Mesh* ice_mesh_):
            GridData(ice_mesh_)
        {};

        inline NodeData(INMOST::Mesh* ice_mesh_, int layer_):
            GridData(ice_mesh_, layer_)
        {};

        // Create data on all nodes
        inline virtual void Create(const mesh::meshVar& pNot,  const mesh::meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::NODE, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const mesh::meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::NODE, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::NODE, INMOST::NONE);};

        // Delete data on all nodes
        inline void Delete(const mesh::meshVar& pNot) override
        {GridDeleteData(pNot, INMOST::NODE);};

        inline void Delete(const std::string& tVar) override
        {GridDeleteData(tVar, INMOST::NODE);};

        inline void Delete(const std::vector<mesh::meshVar>& pNots) override
        {GridDeleteData(pNots, INMOST::NODE);};

        inline void Delete(const std::vector<std::string>& tVars) override
        {GridDeleteData(tVars, INMOST::NODE);};

        // Exchane data
        inline void Exchange(const mesh::meshVar& pNot)
        {ice_mesh->ExchangeData(prog_data.at(pNot).first, INMOST::NODE, 0);};

        inline void Exchange(const std::string& tName)
        {ice_mesh->ExchangeData(temp_data.at(tName).first, INMOST::NODE, 0);};
    };

    // Edge data inharited class
    class EdgeData: public GridData
    {
    public:
        inline EdgeData(INMOST::Mesh* ice_mesh_):
            GridData(ice_mesh_)
        {};

        inline EdgeData(INMOST::Mesh* ice_mesh_, int layer_):
            GridData(ice_mesh_, layer_)
        {};

        // Create data on all edges
        inline virtual void Create(const mesh::meshVar& pNot,  const mesh::meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::FACE, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const mesh::meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::FACE, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::FACE, INMOST::NONE);};

        // Delete data on all edges
        inline void Delete(const mesh::meshVar& pNot) override
        {GridDeleteData(pNot, INMOST::FACE);};

        inline void Delete(const std::string& tVar) override
        {GridDeleteData(tVar, INMOST::FACE);};

        inline void Delete(const std::vector<mesh::meshVar>& pNots) override
        {GridDeleteData(pNots, INMOST::FACE);};

        inline void Delete(const std::vector<std::string>& tVars) override
        {GridDeleteData(tVars, INMOST::FACE);};

        // Exchane data
        inline void Exchange(const mesh::meshVar& pNot)
        {ice_mesh->ExchangeData(prog_data.at(pNot).first, INMOST::FACE, 0);};

        inline void Exchange(const std::string& tName)
        {ice_mesh->ExchangeData(temp_data.at(tName).first, INMOST::FACE, 0);};
    };

    // Triangle data inharited class
    class TrianData: public GridData
    {
    public:
        inline TrianData(INMOST::Mesh* ice_mesh_):
            GridData(ice_mesh_)
        {};

        inline TrianData(INMOST::Mesh* ice_mesh_, int layer_):
            GridData(ice_mesh_, layer_)
        {};

        // Create data on all triangles
        inline virtual void Create(const mesh::meshVar& pNot,  const mesh::meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::CELL, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const mesh::meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::CELL, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::CELL, INMOST::NONE);};

        // Delete data on all triangles
        inline void Delete(const mesh::meshVar& pNot) override
        {GridDeleteData(pNot, INMOST::CELL);};

        inline void Delete(const std::string& tVar) override
        {GridDeleteData(tVar, INMOST::CELL);};

        inline void Delete(const std::vector<mesh::meshVar>& pNots) override
        {GridDeleteData(pNots, INMOST::CELL);};

        inline void Delete(const std::vector<std::string>& tVars) override
        {GridDeleteData(tVars, INMOST::CELL);};

        // Exchane data
        inline void Exchange(const mesh::meshVar& pNot)
        {ice_mesh->ExchangeData(prog_data.at(pNot).first, INMOST::CELL, 0);};

        inline void Exchange(const std::string& tName)
        {ice_mesh->ExchangeData(temp_data.at(tName).first, INMOST::CELL, 0);};
    };

    // Boundary node data inharited class
    class BndNodeData: public NodeData
    {
    public:
        inline BndNodeData(INMOST::Mesh* ice_mesh_):
            NodeData(ice_mesh_)
        {};

        inline BndNodeData(INMOST::Mesh* ice_mesh_, int layer_):
            NodeData(ice_mesh_, layer_)
        {};

        // Create on all boundary nodes
        inline void Create(const mesh::meshVar& pNot,  const mesh::meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::NODE, INMOST::NODE);};

        inline void Create(const std::string& tVar, const mesh::meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::NODE, INMOST::NODE);};

        inline void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::NODE, INMOST::NODE);};
    };

    // Boundary edge data inharited class
    class BndEdgeData: public EdgeData
    {
    public:
        inline BndEdgeData(INMOST::Mesh* ice_mesh_):
            EdgeData(ice_mesh_)
        {};

        inline BndEdgeData(INMOST::Mesh* ice_mesh_, int layer_):
            EdgeData(ice_mesh_, layer_)
        {};

        // Create data on all boundary edges
        inline void Create(const mesh::meshVar& pNot,  const mesh::meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::FACE, INMOST::FACE);};

        inline void Create(const std::string& tVar, const mesh::meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::FACE, INMOST::FACE);};

        inline void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::FACE, INMOST::FACE);};
    };

    // Boundary edge data inharited class
    class BndTrianData: public TrianData
    {
    public:
        inline BndTrianData(INMOST::Mesh* ice_mesh_):
            TrianData(ice_mesh_)
        {};

        inline BndTrianData(INMOST::Mesh* ice_mesh_, int layer_):
            TrianData(ice_mesh_, layer_)
        {};

        // Create data on all boundary triangles
        inline void Create(const mesh::meshVar& pNot,  const mesh::meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::CELL, INMOST::CELL);};

        inline void Create(const std::string& tVar, const mesh::meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::CELL, INMOST::CELL);};

        inline void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::CELL, INMOST::CELL);};
    };
}