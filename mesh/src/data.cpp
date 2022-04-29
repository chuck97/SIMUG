#include "data.hpp"

using namespace SIMUG::mesh;
using namespace INMOST;
using namespace std;

// Create prognostic data on mesh (scalar, vector or tensor)
void GridData::GridCreateData(const meshVar& pNot, const meshDim& pDim,
                              const INMOST::DataType& InmostDataType,
                              const unsigned char& GridElem,
                              const unsigned char& GridSparse)
{
    INMOST::Tag tag;
    switch(pDim)
    {
        case scalar:
            tag = ice_mesh->CreateTag(meshVarName.at(pNot), InmostDataType, GridElem, GridSparse, 1);
            break;
                
        case vector:
            tag = ice_mesh->CreateTag(meshVarName.at(pNot), InmostDataType, GridElem, GridSparse, 3);
            break;
                
        case tensor:
            tag = ice_mesh->CreateTag(meshVarName.at(pNot), InmostDataType, GridElem, GridSparse, 4);
    }
    ProgData_[pNot] = tag;
};

// Create temporal data on mesh (scalar, vector or tensor)
void GridData::GridCreateData(const std::string& tVar, const meshDim& tDim,
                              const INMOST::DataType& InmostDataType,
                              const unsigned char& GridElem,
                              const unsigned char& GridSparse)
{
    INMOST::Tag tag;
    switch(tDim)
    {
        case scalar:
            tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, 1);
            break;
                
        case vector:
            tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, 3);
            break;
                
        case tensor:
            tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, 4);
    }
    TempData_[tVar] = tag;
};

void GridData::GridCreateData(const std::string& tVar, const int& vSize,
                              const INMOST::DataType& InmostDataType,
                              const unsigned char& GridElem,
                              const unsigned char& GridSparse)
{
    INMOST::Tag tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, vSize);
    TempData_[tVar] = tag;
};