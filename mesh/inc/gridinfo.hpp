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
    // grid information for particular element
    class GridInfo
    {
    public:
        inline GridInfo(INMOST::Mesh* ice_mesh_):ice_mesh(ice_mesh_) {};

    protected:
        INMOST::Mesh* ice_mesh;
        INMOST::Tag id;
        INMOST::Tag id_no_bnd;
        INMOST::Tag is_bnd;
        std::map<SIMUG::coord::coordType, INMOST::Tag> coords;
    
    public:
        virtual void Exchange() = 0;
        virtual void Mute() = 0;
        virtual void UnMute(INMOST::Mesh* ice_mesh) = 0;
        virtual void AssignCoords() = 0;
        //virtual void AssempleBasisData() = 0;
    
    protected:
        void ExchangeAll(INMOST::Mesh* ice_mesh, const unsigned char& GridElem)
        {
            ice_mesh->ExchangeData(id, GridElem, 0);
            ice_mesh->ExchangeData(id_no_bnd, GridElem, 0);
            ice_mesh->ExchangeData(is_bnd, GridElem, 0);

            for (auto [key, value]: coords)
                ice_mesh->ExchangeData(value, GridElem, 0);
            BARRIER
        };

        void MuteAll(INMOST::Mesh* ice_mesh);
        void UnMuteAll(INMOST::Mesh* ice_mesh);

        void SelectBndNodes();
        void SelectBndEdges();
        void SelectBndTrians();
        
    };

    class NodeInfo : public GridInfo
    {
    public:
        NodeInfo(INMOST::Mesh* ice_mesh_):
            GridInfo(ice_mesh_)
        {
            SIMUG::Logger mesh_log(cout);
            SIMUG::Timer mesh_timer;
            double duration;

            mesh_timer.Launch();
            SelectBndEdges();
            mesh_timer.Stop();
            duration = mesh_timer.GetMaxTime();
            mesh_timer.Reset();
            if (ice_mesh->GetProcessorRank()==0)
            mesh_log.Log("Boundary edges selected succeessfully! (" + to_string(duration) + " ms)\n");
            BARRIER
        };

        inline void Exchange(INMOST::Mesh* ice_mesh)
        {ExchangeAll(ice_mesh, INMOST::NODE);};

        void Mute(INMOST::Mesh* ice_mesh);
        void UnMute(INMOST::Mesh* ice_mesh);

        inline INMOST::ElementArray<INMOST::Node>& GetBndNodes() {return bnd_nodes;};
        inline const INMOST::ElementArray<INMOST::Node>& GetBndNodes() const {return bnd_nodes;};

    private:
        INMOST::ElementArray<INMOST::Node> bnd_nodes;
    };

    class EdgeInfo : public GridInfo
    {
    public:
        EdgeInfo(INMOST::Mesh* ice_mesh);

        inline void Exchange(INMOST::Mesh* ice_mesh)
        {ExchangeAll(ice_mesh, INMOST::FACE);};

        inline INMOST::ElementArray<INMOST::Face>& GetBndEdges() {return bnd_edges;};
        inline const INMOST::ElementArray<INMOST::Face>& GetBndEdges() const {return bnd_edges;};

    private:
        INMOST::ElementArray<INMOST::Face> bnd_edges;
        
    };

    class TrianInfo : public GridInfo
    {
    public:
        TrianInfo(INMOST::Mesh* ice_mesh);

        inline void Exchange(INMOST::Mesh* ice_mesh)
        {ExchangeAll(ice_mesh, INMOST::CELL);};

        inline INMOST::ElementArray<INMOST::Cell>& GetBndTrians() {return bnd_trians;};
        inline const INMOST::ElementArray<INMOST::Cell>& GetBndTrians() const {return bnd_trians;};
    
    private:
        INMOST::ElementArray<INMOST::Cell> bnd_trians;
    };
}