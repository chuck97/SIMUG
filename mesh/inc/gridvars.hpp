#pragma once
#include "mesh_data.hpp"
#include "defines.hpp"

#include <map>

namespace SIMUG::mesh
{
    // grid element types
    enum gridElemType
    {
        Node,
        Edge,
        Trian,
        bndNode,
        bndEdge,
        bndTrian
    };

    // prognostic and forcing variable dimensions
    static const std::map<meshVar, meshDim> progDims = 
    {
        {meshVar::mi, meshDim::scalar},
        {meshVar::hi, meshDim::scalar},
        {meshVar::ai, meshDim::scalar},
        {meshVar::P0, meshDim::scalar},
        {meshVar::del, meshDim::scalar},
        {meshVar::ui, meshDim::vector},
        {meshVar::sig, meshDim::tensor},
        {meshVar::eps, meshDim::tensor}
    };

    static const std::map<meshVar, meshDim> forcDims = 
    {
        {meshVar::hw, meshDim::scalar},
        {meshVar::ua, meshDim::vector},
        {meshVar::uw, meshDim::vector}
    };

    // Agrid elements for prognostic and forcing variables
    static const std::map<meshVar, gridElemType> gridA_progElems = 
    {
        {meshVar::mi, gridElemType::Node},
        {meshVar::hi, gridElemType::Node},
        {meshVar::ai, gridElemType::Node},
        {meshVar::P0, gridElemType::Node},
        {meshVar::del, gridElemType::Trian},
        {meshVar::ui, gridElemType::Node},
        {meshVar::sig, gridElemType::Trian},
        {meshVar::eps, gridElemType::Trian}
    };

    static const std::map<meshVar, gridElemType> gridA_forcElems
    {
        {meshVar::hw, gridElemType::Node},
        {meshVar::ua, gridElemType::Node},
        {meshVar::uw, gridElemType::Node}
    };
}