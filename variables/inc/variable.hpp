#pragma once

#include <iostream>
#include <string>

#include "mesh_data.hpp"
#include "mesh_info.hpp"
#include "mesh_info.hpp"
#include "model_var.hpp"

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
        inline Var():
            varName_(""),
            varType_(varNoType)
        {};

        inline Var(const std::string& _varName_,
                   const varType&     _varType_):
            varName_(_varName_),
            varType_(_varType_)
        {};

        inline void SetName(const std::string& _varName_)
        {varName_ = _varName_;};

        inline const std::string& GetName() const
        {return varName_;};

        virtual inline const varType& GetType() const = 0;
    };

    // mesh data variable class
    class MeshDataVar : public Var
    {
    private:
        std::string varUnit_;
        mesh::meshDim varDim_;

    public:
        inline MeshDataVar():
            Var("", varType::varMeshData)
        {};

        inline MeshDataVar(const std::string&    _varName_, 
                           const std::string&    _varUnit_,
                           const mesh::meshDim&  _varDim_ ):
            Var(_varName_, varType::varMeshData),
            varUnit_(_varUnit_),
            varDim_(_varDim_)
        {};

        inline void SetUnit(const std::string&  _varUnit_)
        {varUnit_ = _varUnit_;};

        inline void SetDim(const mesh::meshDim& _varDim_)
        {varDim_ = _varDim_;};

        inline const varType&     GetType() const
        {return varType_;};

        inline const std::string& GetUnit() const
        {return varUnit_;};

        inline const mesh::meshDim&  GetDim() const
        {return varDim_;};
    };

    std::ostream& operator<< (std::ostream& out, const MeshDataVar& var)
    {
        out << "Mesh variable \'" << var.GetName()
        << "\': unit = \'" << var.GetUnit()
        << "\', dim = \'" << mesh::meshDimName.at(var.GetDim()) << "\';\n";
        return out;
    }
    
    // mesh information variable class
    class MeshInfoVar: public Var
    {
    private:
        mesh::surfType surfType_;
        mesh::gridType gridType_;
    
    public:
        inline MeshInfoVar():
            Var("", varType::varMeshInfo)
        {};

        inline MeshInfoVar(const std::string& _varName_,
                           const mesh::surfType&    _surfType_,
                           const mesh::gridType&    _gridType_):
            Var(_varName_, varType::varMeshInfo),
            surfType_(_surfType_),
            gridType_(_gridType_)
        {};

        inline void SetSurfType(const mesh::surfType& _surfType_)
        {surfType_ = _surfType_;};

        inline void SetGridType(const mesh::gridType& _gridType_)
        {gridType_ = _gridType_;};

        inline const varType&  GetType() const
        {return varType_;};

        inline const mesh::surfType& GetSurfType() const
        {return surfType_;};

        inline const mesh::gridType& GetGridType() const
        {return gridType_;};
    };

    std::ostream& operator<< (std::ostream& out, const MeshInfoVar& var)
    {
        out << "Mesh type variable \'" << var.GetName()
        << "\': surface type = \'" << mesh::surfTypeName.at(var.GetSurfType())
        << "\', grid type = \'" << mesh::gridTypeName.at(var.GetGridType()) << "\';\n";
        return out;
    };

    // model variable class
    template<typename T>
    class ModelVar: public Var
    {
    private:
        std::string varUnit_;
        T           varValue_;
    public:
        inline ModelVar():
            Var("", varType::varModel)
        {};

        inline ModelVar(const std::string& _varName_,
                        const std::string& _varUnit_,
                        const T&           _varValue_):
            Var(_varName_, varType::varModel),
            varUnit_(_varUnit_),
            varValue_(_varValue_)
        {};

        inline void SetUnit(const std::string& _varUnit_)
        {varUnit_ = _varUnit_;};

        inline void SetValue(const T& _varValue_)
        {varValue_ = _varValue_;};

        inline const varType& GetType() const
        {return varType_;};

        inline const std::string& GetUnit() const
        {return varUnit_;};

        inline const T& GetValue() const
        {return varValue_;};
    };

    template class ModelVar<int>;
    template class ModelVar<float>;
    template class ModelVar<double>;
    template class ModelVar<std::string>;

    template <typename T>
    std::ostream& operator<< (std::ostream& out, const ModelVar<T>& var)
    {
        out << "Model variable \'" << var.GetName()
        << "\': value = \'" << var.GetValue()
        << "\', unit = \'" << var.GetUnit() << "\';\n";
        return out;
    }

    
}