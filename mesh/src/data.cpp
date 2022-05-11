#include "data.hpp"

using namespace SIMUG::mesh;
using namespace INMOST;
using namespace std;

// Create prognostic grid data (scalar, vector or tensor)
void GridData::GridCreateData(const meshVar& pNot, const meshDim& pDim,
                              const INMOST::DataType& InmostDataType,
                              const unsigned char& GridElem,
                              const unsigned char& GridSparse)
{
    INMOST::Tag tag;
    switch(pDim)
    {
        case scalar:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(meshVarName.at(pNot) + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 1);
            else
                tag = ice_mesh->CreateTag(meshVarName.at(pNot), InmostDataType, GridElem, GridSparse, 1);
            break;
                
        case vector:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(meshVarName.at(pNot) + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 3);
            else
                tag = ice_mesh->CreateTag(meshVarName.at(pNot), InmostDataType, GridElem, GridSparse, 3);
            break;
                
        case tensor:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(meshVarName.at(pNot) + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 4);
            else
                tag = ice_mesh->CreateTag(meshVarName.at(pNot), InmostDataType, GridElem, GridSparse, 4);
    }
    prog_data[pNot] = tag;
};

// Create temporal grid data (scalar, vector or tensor)
void GridData::GridCreateData(const std::string& tVar, const meshDim& tDim,
                              const INMOST::DataType& InmostDataType,
                              const unsigned char& GridElem,
                              const unsigned char& GridSparse)
{
    INMOST::Tag tag;
    switch(tDim)
    {
        case scalar:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(tVar + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 1);
            else
                tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, 1);
            break;
                
        case vector:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(tVar + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 3);
            else
                tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, 3);
            break;
                
        case tensor:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(tVar + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 4);
            else
                tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, 4);
    }
    temp_data[tVar] = tag;
};

// Create temporal grid data with given size
void GridData::GridCreateData(const std::string& tVar, const int& vSize,
                              const INMOST::DataType& InmostDataType,
                              const unsigned char& GridElem,
                              const unsigned char& GridSparse)
{
    INMOST::Tag tag;
    if (layer.has_value())
        tag = ice_mesh->CreateTag(tVar + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, vSize);
    else
        tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, vSize);

    temp_data[tVar] = tag;
};