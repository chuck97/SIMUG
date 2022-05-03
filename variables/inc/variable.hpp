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
        varScal,
        varVec,
        varTens,
        varNoType, 
        varValue,
        varNoUnit   
    };

    // variable base class
    class VarBase
    {
    protected:
        std::string varName_;
        varType varType_;

    public:
        inline VarBase():
            varName_(""),
            varType_(varNoType)
        {};

        inline VarBase(const std::string& _varName_,
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
    template <typename T>
    class VarNoUnit : public VarBase
    {
    protected:
        std::string varValue_;

    public:
        inline VarNoUnit():
            VarBase("", varType::varValueNoUnit)
        {};

        inline VarNoUnit(const std::string&    _varName_, 
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
    class ValueVar: public Var
    {
    private:
        std::string varUnit_;
        T           varValue_;
    public:
        inline ValueVar():
            Var("", varType::varModel)
        {};

        inline ValueVar(const std::string& _varName_,
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

    template class ValueVar<int>;
    template class ValueVar<float>;
    template class ValueVar<double>;
    template class ValueVar<std::string>;

    template <typename T>
    std::ostream& operator<< (std::ostream& out, const ValueVar<T>& var)
    {
        out << "Model variable \'" << var.GetName()
        << "\': value = \'" << var.GetValue()
        << "\', unit = \'" << var.GetUnit() << "\';\n";
        return out;
    }

    // physics consts class
    template <typename T>
    class PhysConstVar: public Var
    {
    private:
        std::string varUnit_;
        T           varValue_;
    public:
        inline PhysConstVar():
            Var("", varType::varPhysConst)
        {};

        inline PhysConstVar(const std::string& _varName_,
                            const std::string& _varName_)
    }
    
}