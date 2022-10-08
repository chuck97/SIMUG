#include "data.hpp"

using namespace INMOST;
using namespace std;

namespace SIMUG
{

// Create prognostic grid data (scalar, vector or tensor)
void GridData::GridCreateData(const mesh::meshVar& pNot, const mesh::meshDim& pDim,
                              const INMOST::DataType& InmostDataType,
                              const unsigned char& GridElem,
                              const unsigned char& GridSparse)
{
    INMOST::Tag tag;
    switch(pDim)
    {
        case mesh::scalar:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(mesh::meshVarName.at(pNot) + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 1);
            else
                tag = ice_mesh->CreateTag(mesh::meshVarName.at(pNot), InmostDataType, GridElem, GridSparse, 1);
            prog_data_size[pNot] = 1;
            break;
                
        case mesh::vector:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(mesh::meshVarName.at(pNot) + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 3);
            else
                tag = ice_mesh->CreateTag(mesh::meshVarName.at(pNot), InmostDataType, GridElem, GridSparse, 3);
            prog_data_size[pNot] = 3;
            break;
                
        case mesh::tensor:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(mesh::meshVarName.at(pNot) + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 4);
            else
                tag = ice_mesh->CreateTag(mesh::meshVarName.at(pNot), InmostDataType, GridElem, GridSparse, 4);
            prog_data_size[pNot] = 4;
    }
    prog_data[pNot].first = tag;
};

// Create temporal grid data (scalar, vector or tensor)
void GridData::GridCreateData(const std::string& tVar, const mesh::meshDim& tDim,
                              const INMOST::DataType& InmostDataType,
                              const unsigned char& GridElem,
                              const unsigned char& GridSparse)
{
    INMOST::Tag tag;
    switch(tDim)
    {
        case mesh::scalar:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(tVar + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 1);
            else
                tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, 1);
            temp_data_size[tVar] = 1;
            break;
                
        case mesh::vector:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(tVar + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 3);
            else
                tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, 3);
            temp_data_size[tVar] = 3;
            break;
                
        case mesh::tensor:
            if (layer.has_value())
                tag = ice_mesh->CreateTag(tVar + " " + std::to_string(layer.value()), InmostDataType, GridElem, GridSparse, 4);
            else
                tag = ice_mesh->CreateTag(tVar, InmostDataType, GridElem, GridSparse, 4);
            temp_data_size[tVar] = 4;
    }
    temp_data[tVar].first = tag;
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
    
    temp_data_size[tVar] = vSize;
    temp_data[tVar].first = tag;
};

}