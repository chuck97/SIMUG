#pragma once
#include "inmost.h"

namespace SIMUG::mesh
{
    template <typename RT>
    struct ElemBasisData
    {
        // transition matrix nodal -> element (three 2x2 matricies written in array with size 12)
        INMOST::Tag trans_nodal_to_elem;

        // transition matrix element -> nodal (three 2x2 matricies written in array with size 12)
        INMOST::Tag trans_elem_to_nodal;

        // nodal coords in element basis (three 2x1 vectors written in array with size 6)
        INMOST::Tag node_coords_elem_basis;

        // nodal coords in nodal bases 
        // v_ijk: i={0,1,2} - node number, j={0(x),1(y)}, k={0,1,2} - nodal basis number
        // result is written in array with size 18
        INMOST::Tag node_coords_nodal_basis;
    };
}