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
    enum gridElemType
    {
        Node,
        Edge,
        Trian,
        bndNode,
        bndEdge,
        bndTrian
    };


    // Grid Data base class
    class GridData
    {
    public:
        typedef std::map<meshVar, INMOST::Tag> ProgData;
        typedef std::map<std::string, INMOST::Tag> TempData;

    protected:  
        INMOST::Mesh* ice_mesh;    

        ProgData prog_data;
        TempData temp_data;

    public:
        inline GridData(INMOST::Mesh* ice_mesh_):
            ice_mesh(ice_mesh_)
        {};

        // Create mesh data
        virtual void Create(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) = 0;
        virtual void Create(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) = 0;
        virtual void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) = 0;

        // Obtain mesh data
        inline INMOST::Tag& Get(const meshVar& pNot)
        {return prog_data.at(pNot);};

        inline INMOST::Tag& Get(const std::string& tName)
        {return temp_data.at(tName);};

        inline const INMOST::Tag& Get(const meshVar& pNot) const
        {return prog_data.at(pNot);};

        inline const INMOST::Tag& Get(const std::string& tName) const 
        {return temp_data.at(tName);};

        // Delete mesh data
        virtual void Delete(const meshVar& pNot) = 0;
        virtual void Delete(const std::string& tVar) = 0;
        virtual void Delete(const std::vector<meshVar>& pNots) = 0;
        virtual void Delete(const std::vector<std::string>& tVars) = 0;

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
        {if (prog_data.count(pNot)==0) {SIMUG_ERR("variable "+meshVarName.at(pNot)+" doesn't exist")};
        ice_mesh->DeleteTag(prog_data.at(pNot), InmostGridElem);};

        inline void GridDeleteData(const std::string& tVar, unsigned char InmostGridElem)
        {if (temp_data.count(tVar)==0) {SIMUG_ERR("variable "+tVar+" doesn't exist")};
        ice_mesh->DeleteTag(temp_data.at(tVar), InmostGridElem);};

        inline void GridDeleteData(const std::vector<meshVar>& pNots, unsigned char InmostGridElem)
        {for (const auto& item: pNots){if (prog_data.count(item)==0) {SIMUG_ERR("variable "+meshVarName.at(item)+" doesn't exist")};
        ice_mesh->DeleteTag(prog_data.at(item), InmostGridElem);};};

        inline void GridDeleteData(const std::vector<std::string>& tVars, unsigned char InmostGridElem)
        {for (const auto& item: tVars){if (temp_data.count(item)==0) {SIMUG_ERR("variable "+item+" doesn't exist")};
        ice_mesh->DeleteTag(temp_data.at(item), InmostGridElem);};};
   
    };


    // Node data inharited class
    class NodeData: public GridData
    {
    public:
        inline NodeData(INMOST::Mesh* ice_mesh_):
            GridData(ice_mesh_)
        {};

        // Create data on all nodes
        inline virtual void Create(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::NODE, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::NODE, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::NODE, INMOST::NONE);};

        // Delete data on all nodes
        inline void Delete(const meshVar& pNot) override
        {GridDeleteData(pNot, INMOST::NODE);};

        inline void Delete(const std::string& tVar) override
        {GridDeleteData(tVar, INMOST::NODE);};

        inline void Delete(const std::vector<meshVar>& pNots) override
        {GridDeleteData(pNots, INMOST::NODE);};

        inline void Delete(const std::vector<std::string>& tVars) override
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
        inline virtual void Create(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::FACE, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::FACE, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::FACE, INMOST::NONE);};

        // Delete data on all edges
        inline void Delete(const meshVar& pNot) override
        {GridDeleteData(pNot, INMOST::FACE);};

        inline void Delete(const std::string& tVar) override
        {GridDeleteData(tVar, INMOST::FACE);};

        inline void Delete(const std::vector<meshVar>& pNots) override
        {GridDeleteData(pNots, INMOST::FACE);};

        inline void Delete(const std::vector<std::string>& tVars) override
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
        inline virtual void Create(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::CELL, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::CELL, INMOST::NONE);};

        inline virtual void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::CELL, INMOST::NONE);};

        // Delete data on all triangles
        inline void Delete(const meshVar& pNot) override
        {GridDeleteData(pNot, INMOST::CELL);};

        inline void Delete(const std::string& tVar) override
        {GridDeleteData(tVar, INMOST::CELL);};

        inline void Delete(const std::vector<meshVar>& pNots) override
        {GridDeleteData(pNots, INMOST::CELL);};

        inline void Delete(const std::vector<std::string>& tVars) override
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
        inline void Create(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::NODE, INMOST::NODE);};

        inline void Create(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
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

        // Create data on all boundary edges
        inline void Create(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::FACE, INMOST::FACE);};

        inline void Create(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::FACE, INMOST::FACE);};

        inline void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
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
        inline void Create(const meshVar& pNot,  const meshDim& pDim, const INMOST::DataType& dtype) override
        {GridCreateData(pNot, pDim, dtype, INMOST::CELL, INMOST::CELL);};

        inline void Create(const std::string& tVar, const meshDim& tDim, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, tDim, dtype, INMOST::CELL, INMOST::CELL);};

        inline void Create(const std::string& tVar, const int& vSize, const INMOST::DataType& dtype) override
        {GridCreateData(tVar, vSize, dtype, INMOST::CELL, INMOST::CELL);};
    };
}