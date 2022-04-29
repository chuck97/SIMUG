#pragma once

#include <map> 
#include <vector>
#include <string>
#include <memory>
#include "inmost.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "mesh_data.hpp"
#include "defines.hpp"


namespace SIMUG::mesh
{
    // Grid Data base class
    class GridData
    {
        typedef std::map<meshVar, INMOST::Tag> ProgData;
        typedef std::map<std::string, INMOST::Tag> TempData;

    protected:  
        INMOST::Mesh* ice_mesh;    

        ProgData ProgData_;
        TempData TempData_;

    public:
        inline GridData(INMOST::Mesh* ice_mesh_):
            ice_mesh(ice_mesh_)
        {};

        // Create mesh data
        virtual void CreateData(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) = 0;
        virtual void CreateData(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) = 0;
        virtual void CreateData(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) = 0;

        // Obtain mesh data
        inline INMOST::Tag& GetData(const meshVar& pNot)
        {return ProgData_.at(pNot);};

        inline INMOST::Tag& GetData(const std::string& tName)
        {return TempData_.at(tName);};

        inline const INMOST::Tag& GetData(const meshVar& pNot) const
        {return ProgData_.at(pNot);};

        inline const INMOST::Tag& GetData(const std::string& tName) const 
        {return TempData_.at(tName);};

        // Delete mesh data
        virtual void DeleteData(const meshVar& pNot) = 0;
        virtual void DeleteData(const std::string& tVar) = 0;
        virtual void DeleteData(const std::vector<meshVar>& pNots) = 0;
        virtual void DeleteData(const std::vector<std::string>& tVars) = 0;

    protected:
        void GridCreateData(const meshVar& pNot,  const meshDim& pDim, 
                            const INMOST::DataType& InmostDataType,
                            const unsigned char& GridElem,
                            const unsigned char& GridSparse);

        void GridCreateData(const std::string& tVar, const meshDim& tDim,
                            const INMOST::DataType& InmostDataType,
                            const unsigned char& GridElem,
                            const unsigned char& GridSparse);

        void GridCreateData(const std::string& tVar, const int& vSize,
                            const INMOST::DataType& InmostDataType,
                            const unsigned char& GridElem,
                            const unsigned char& GridSparse);

        inline void GridDeleteData(const meshVar& pNot, unsigned char InmostGridElem)
        {if (ProgData_.count(pNot)==0) {SIMUG_ERR("variable "+meshVarName.at(pNot)+" doesn't exist")};
        ice_mesh->DeleteTag(ProgData_.at(pNot), InmostGridElem);};

        inline void GridDeleteData(const std::string& tVar, unsigned char InmostGridElem)
        {if (TempData_.count(tVar)==0) {SIMUG_ERR("variable "+tVar+" doesn't exist")};
        ice_mesh->DeleteTag(TempData_.at(tVar), InmostGridElem);};

        inline void GridDeleteData(const std::vector<meshVar>& pNots, unsigned char InmostGridElem)
        {for (const auto& item: pNots){if (ProgData_.count(item)==0) {SIMUG_ERR("variable "+meshVarName.at(item)+" doesn't exist")};
        ice_mesh->DeleteTag(ProgData_.at(item), InmostGridElem);};};

        inline void GridDeleteData(const std::vector<std::string>& tVars, unsigned char InmostGridElem)
        {for (const auto& item: tVars){if (TempData_.count(item)==0) {SIMUG_ERR("variable "+item+" doesn't exist")};
        ice_mesh->DeleteTag(TempData_.at(item), InmostGridElem);};};
    };


    // Node data inharited class
    class NodeData: public GridData
    {
    public:
        inline NodeData(INMOST::Mesh* ice_mesh_):
            GridData(ice_mesh_)
        {};

        // Create data on all nodes
        inline virtual void CreateData(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::NODE, INMOST::NONE);};

        inline virtual void CreateData(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::NODE, INMOST::NONE);};

        inline virtual void CreateData(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::NODE, INMOST::NONE);};

        // Delete data on all nodes
        inline void DeleteData(const meshVar& pNot) override
        {GridDeleteData(pNot, INMOST::NODE);};

        inline void DeleteData(const std::string& tVar) override
        {GridDeleteData(tVar, INMOST::NODE);};

        inline void DeleteData(const std::vector<meshVar>& pNots) override
        {GridDeleteData(pNots, INMOST::NODE);};

        inline void DeleteData(const std::vector<std::string>& tVars) override
        {GridDeleteData(tVars, INMOST::NODE);};
    };

    // Edge data inharited class
    class EdgeData: public GridData
    {
    public:
        inline EdgeData(INMOST::Mesh* ice_mesh_):
            GridData(ice_mesh_)
        {};

        // Create data on all edges
        inline virtual void CreateData(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::FACE, INMOST::NONE);};

        inline virtual void CreateData(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::FACE, INMOST::NONE);};

        inline virtual void CreateData(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::FACE, INMOST::NONE);};

        // Delete data on all edges
        inline void DeleteData(const meshVar& pNot) override
        {GridDeleteData(pNot, INMOST::FACE);};

        inline void DeleteData(const std::string& tVar) override
        {GridDeleteData(tVar, INMOST::FACE);};

        inline void DeleteData(const std::vector<meshVar>& pNots) override
        {GridDeleteData(pNots, INMOST::FACE);};

        inline void DeleteData(const std::vector<std::string>& tVars) override
        {GridDeleteData(tVars, INMOST::FACE);};
    };

    // Triangle data inharited class
    class TriangleData: public GridData
    {
    public:
        inline TriangleData(INMOST::Mesh* ice_mesh_):
            GridData(ice_mesh_)
        {};

        // Create data on all triangles
        inline virtual void CreateData(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::CELL, INMOST::NONE);};

        inline virtual void CreateData(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::CELL, INMOST::NONE);};

        inline virtual void CreateData(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::CELL, INMOST::NONE);};

        // Delete data on all triangles
        inline void DeleteData(const meshVar& pNot) override
        {GridDeleteData(pNot, INMOST::CELL);};

        inline void DeleteData(const std::string& tVar) override
        {GridDeleteData(tVar, INMOST::CELL);};

        inline void DeleteData(const std::vector<meshVar>& pNots) override
        {GridDeleteData(pNots, INMOST::CELL);};

        inline void DeleteData(const std::vector<std::string>& tVars) override
        {GridDeleteData(tVars, INMOST::CELL);};
    };

    // Boundary node data inharited class
    class BndNodeData: public NodeData
    {
    public:
        inline BndNodeData(INMOST::Mesh* ice_mesh_):
            NodeData(ice_mesh_)
        {};

        // Create on all boundary nodes
        inline void CreateData(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::NODE, INMOST::NODE);};

        inline void CreateData(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::NODE, INMOST::NODE);};

        inline void CreateData(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::NODE, INMOST::NODE);};
    };

    // Boundary edge data inharited class
    class BndEdgeData: public EdgeData
    {
    public:
        inline BndEdgeData(INMOST::Mesh* ice_mesh_):
            EdgeData(ice_mesh_)
        {};

        // Create data on all boundary edges
        inline void CreateData(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::FACE, INMOST::FACE);};

        inline void CreateData(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::FACE, INMOST::FACE);};

        inline void CreateData(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::FACE, INMOST::FACE);};
    };

    // Boundary edge data inharited class
    class BndTriangleData: public TriangleData
    {
    public:
        inline BndTriangleData(INMOST::Mesh* ice_mesh_):
            TriangleData(ice_mesh_)
        {};

        // Create data on all boundary triangles
        inline void CreateData(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::CELL, INMOST::CELL);};

        inline void CreateData(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::CELL, INMOST::CELL);};

        inline void CreateData(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::CELL, INMOST::CELL);};
    };
}