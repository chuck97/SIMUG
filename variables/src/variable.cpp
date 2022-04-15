#include "variable.hpp"

using namespace std;
using namespace SIMUG;

Var::Var()
{
}

Var::Var(const std::string& _varName_,
         const varType& _varType_):
    varName_(_varName_),
    varType_(_varType_)
{
}

inline void Var::SetName(const std::string& _varName_)
{
    varName_ = _varName_;
}

inline const std::string& Var::GetName() const
{
    return varName_;
}

MeshDataVar::MeshDataVar():
    Var("", varType::varMeshData)
{
}

MeshDataVar::MeshDataVar(const std::string&       _varName_, 
                         const std::string&       _varUnit_,
                         const mesh::meshDim&     _varDim_):
    Var(_varName_, varType::varMeshData),
    varUnit_(_varUnit_),
    varDim_(_varDim_)
{
}

inline void MeshDataVar::SetUnit(const std::string& _varUnit_)
{
    varUnit_ = _varUnit_;
}

inline void MeshDataVar::SetDim(const mesh::meshDim& _varDim_)
{
    varDim_ = _varDim_;
}

inline const varType& MeshDataVar::GetType() const
{
    return varType_;
}

inline const std::string& MeshDataVar::GetUnit() const
{
    return varUnit_;
}

inline const mesh::meshDim&  MeshDataVar::GetDim() const
{
    return varDim_;
}

ostream& operator<< (ostream& out, const MeshDataVar& var)
{
    out << "Mesh variable \'" << var.GetName()
        << "\': unit = \'" << var.GetUnit()
        << "\', dim = \'" << mesh::meshDimName[var.GetDim()] << "\';\n";
    return out;
}