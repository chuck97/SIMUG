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
    static const std::map<meshVar, meshDim> multiDims = 
    {
        {meshVar::mi, meshDim::scalar},
        {meshVar::hi, meshDim::scalar},
        {meshVar::ai, meshDim::scalar}
    };

    static const std::map<meshVar, meshDim> singleDims = 
    {
        {meshVar::P0, meshDim::scalar},
        {meshVar::del, meshDim::scalar},
        {meshVar::sig, meshDim::tensor},
        {meshVar::eps, meshDim::tensor},
        {meshVar::ui, meshDim::vector},
        {meshVar::hw, meshDim::scalar},
        {meshVar::ua, meshDim::vector},
        {meshVar::uw, meshDim::vector}
    };

    // Agrid elements for prognostic and forcing variables
    static const std::map<meshVar, gridElemType> gridA_multiElems = 
    {
        {meshVar::mi, gridElemType::Node},
        {meshVar::hi, gridElemType::Node},
        {meshVar::ai, gridElemType::Node}
    };

    static const std::map<meshVar, gridElemType> gridA_singleElems
    {
        {meshVar::P0, gridElemType::Node},
        {meshVar::del, gridElemType::Trian},
        {meshVar::ui, gridElemType::Node},
        {meshVar::sig, gridElemType::Trian},
        {meshVar::eps, gridElemType::Trian},
        {meshVar::hw, gridElemType::Node},
        {meshVar::ua, gridElemType::Node},
        {meshVar::uw, gridElemType::Node}
    };

    // Cgrid elements for prognostic and forcing variables
    static const std::map<meshVar, gridElemType> gridC_multiElems = 
    {
        {meshVar::mi, gridElemType::Trian},
        {meshVar::hi, gridElemType::Trian},
        {meshVar::ai, gridElemType::Trian}
    };

    static const std::map<meshVar, gridElemType> gridC_singleElems
    {
        {meshVar::P0, gridElemType::Edge},
        {meshVar::del, gridElemType::Trian},
        {meshVar::ui, gridElemType::Edge},
        {meshVar::sig, gridElemType::Trian},
        {meshVar::eps, gridElemType::Trian},
        {meshVar::hw, gridElemType::Trian},
        {meshVar::ua, gridElemType::Edge},
        {meshVar::uw, gridElemType::Edge}
    };
}