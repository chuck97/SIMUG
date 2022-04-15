#pragma once

#include <iostream>
#include <string>

#include "meshdata.hpp"
#include "meshinfo.hpp"

namespace SIMUG
{
    enum varType
    {
        varNoType, 
        varPhysConst,   
        varModel,  
        varAdvection,    
        varDynamics,    
        varMeshData,
        varMeshInfo,   
        varForcing    
    };

    // variable base class
    class Var
    {
    protected:
        std::string varName_;
        varType varType_;

    public:
        Var();

        Var(const std::string& _varName_,
            const varType& _varType_);

        inline void SetName(const std::string& _varName_);
        inline const std::string& GetName() const;

        virtual inline const varType& GetType() const = 0;
    };

    // mesh data variable class
    class MeshDataVar : public Var
    {
    private:
        std::string varUnit_;
        mesh::meshDim varDim_;

    public:
        MeshDataVar();
        MeshDataVar(const std::string&    _varName_, 
                    const std::string&    _varUnit_,
                    const mesh::meshDim&  _varDim_
                    );

        inline void SetUnit(const std::string&  _varUnit_);
        inline void SetDim(const mesh::meshDim& _varDim_);

        inline const varType&     GetType() const;
        inline const std::string& GetUnit() const;
        inline const mesh::meshDim&  GetDim() const;
    };

    std::ostream& operator<< (std::ostream& out, const MeshDataVar& var);

    /*
    // mesh information variable class
    template <typename RT>
    class MeshInfoVar: public Var<RT>
    {
    private:
        surfType surfType_;
        gridType gridType_;
    
    public:
        MeshInfoVar();
        MeshInfoVar(const std::string& _varName_,
                    const surfType&    _surfType_,
                    const gridType&    _gridType_);

        inline void SetSurfType(const surfType& _surfType_);
        inline void SetGridType(const gridType& _gridType_);

        inline varType&  GetType() const;
        inline surfType& GetSurfType() const;
        inline gridType& GetGridType() const;
    };
    */
}