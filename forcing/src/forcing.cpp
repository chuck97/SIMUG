#include "forcing.hpp"

using namespace INMOST;
using namespace SIMUG;

void Forcing::SetAnalytical(mesh::meshVar mesh_var,
                            int layer,
                            FuncPtr func_ptr
                           )
{
    // figure out if it scalar of forcing variable 
    if (mesh::progDims.count(mesh_var) == 1)
    {
        // setup expression
        mesh->GetProgData(mesh->GetMeshInfo().prog_elems[mesh_var], layer)->GetExpr(mesh_var) = func_ptr;
    }
    else if (mesh::forcDims.count(mesh_var) == 1)
    {
        if (layer > 0)
        {
            SIMUG_ERR("forcing variables hawe only one layer!");
        }
        // setup expression
        mesh->GetForcData(mesh->GetMeshInfo().forc_elems[mesh_var])->GetExpr(mesh_var) = func_ptr;
    }
    else
    {
        SIMUG_ERR("unknown mesh variable!");
    }
    BARRIER
}

void Forcing::SetAnalytical(mesh::meshVar mesh_var,
                            FuncPtr func_ptr
                           )
{
    SetAnalytical(mesh_var, 0, func_ptr);
}

void Forcing::SetAnalytical(const std::string& varname, 
                            mesh::gridElemType elem_type,
                            int layer,
                            FuncPtr func_ptr
                            )
{   
    // check that expression output size is equal to variable size
    if(func_ptr({1.0, 1.0}, 1.0).size() != mesh->GetProgData(elem_type, layer)->GetSize(varname))
    {
        SIMUG_ERR("the output for expression should have the same size as variable!");
    }

    mesh->GetProgData(elem_type, layer)->GetExpr(varname) = func_ptr;
    
    BARRIER
}

void Forcing::SetAnalytical(const std::string& varname, 
                            mesh::gridElemType elem_type,
                            FuncPtr func_ptr
                            )
{
    SetAnalytical(varname, elem_type, 0, func_ptr);
    BARRIER
}

void Forcing::UpdateVariable(INMOST::Tag var_tag,
                             mesh::gridElemType elem_type,
                             coord::coordType coord_type,
                             FuncPtr func_ptr,
                             double time)
{
    // output size of expression
    size_t output_size = func_ptr({1.0, 1.0}, 1.0).size();

    INMOST::Mesh* ice_mesh = mesh->GetMesh();

    // update grid variable for every element
    if (elem_type == SIMUG::mesh::gridElemType::Node)
    {
        INMOST::Tag coords_tag = mesh->GetGridInfo(mesh::gridElemType::Node)->coords.at(coord_type);

        for (auto nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit)
        {
            if (nodeit->GetStatus() != Element::Ghost)
            {
                double coord_x = nodeit->RealArray(coords_tag)[0];
                double coord_y = nodeit->RealArray(coords_tag)[1];

                if (output_size == 1)
                {
                    nodeit->Real(var_tag) = func_ptr({coord_x, coord_y}, time)[0];
                }
                else
                {
                    for (size_t i = 0; i < output_size; ++i)
                    {
                        nodeit->RealArray(var_tag)[i] = func_ptr({coord_x, coord_y}, time)[i];
                    }
                }
            }
        }
        ice_mesh->ExchangeData(var_tag, INMOST::NODE, 0);
    }
    else if (elem_type == SIMUG::mesh::gridElemType::Edge)
    {
        INMOST::Tag coords_tag = mesh->GetGridInfo(mesh::gridElemType::Edge)->coords.at(coord_type);

        for (auto edgeit = ice_mesh->BeginFace(); edgeit != ice_mesh->EndFace(); ++edgeit)
        {
            if (edgeit->GetStatus() != Element::Ghost)
            {
                double coord_x = edgeit->RealArray(coords_tag)[0];
                double coord_y = edgeit->RealArray(coords_tag)[1];

                if (output_size == 1)
                {
                    edgeit->Real(var_tag) = func_ptr({coord_x, coord_y}, time)[0];
                }
                else
                {
                    for (size_t i = 0; i < output_size; ++i)
                    {
                        edgeit->RealArray(var_tag)[i] = func_ptr({coord_x, coord_y}, time)[i];
                    }
                }
            }
        }
        ice_mesh->ExchangeData(var_tag, INMOST::FACE, 0);
    }
    else if (elem_type == SIMUG::mesh::gridElemType::Trian)
    {
        INMOST::Tag coords_tag = mesh->GetGridInfo(mesh::gridElemType::Trian)->coords.at(coord_type);

        for (auto trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit)
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                double coord_x = trianit->RealArray(coords_tag)[0];
                double coord_y = trianit->RealArray(coords_tag)[1];

                if (output_size == 1)
                {
                    trianit->Real(var_tag) = func_ptr({coord_x, coord_y}, time)[0];
                }
                else
                {
                    for (size_t i = 0; i < output_size; ++i)
                    {
                        trianit->RealArray(var_tag)[i] = func_ptr({coord_x, coord_y}, time)[i];
                    }
                }
            }
        }
        ice_mesh->ExchangeData(var_tag, INMOST::CELL, 0);
    }
    else if (elem_type == SIMUG::mesh::gridElemType::bndNode)
    {
        INMOST::Tag coords_tag = mesh->GetGridInfo(SIMUG::mesh::gridElemType::Node)->coords.at(coord_type);

        for (size_t bnd_node_num = 0; bnd_node_num < mesh->GetBndNodes().size(); ++bnd_node_num)
        {
            if (mesh->GetBndNodes()[bnd_node_num].GetStatus() != Element::Ghost)
            {
                double coord_x = mesh->GetBndNodes()[bnd_node_num].RealArray(coords_tag)[0];
                double coord_y = mesh->GetBndNodes()[bnd_node_num].RealArray(coords_tag)[1];

                if (output_size == 1)
                {
                    mesh->GetBndNodes()[bnd_node_num].Real(var_tag) = func_ptr({coord_x, coord_y}, time)[0];
                }
                else
                {
                    for (size_t i = 0; i < output_size; ++i)
                    {
                        mesh->GetBndNodes()[bnd_node_num].RealArray(var_tag)[i] = func_ptr({coord_x, coord_y}, time)[i];
                    }
                }
            }
        }
        ice_mesh->ExchangeData(var_tag, INMOST::NODE, 0);
    }
    else if (elem_type == SIMUG::mesh::gridElemType::bndEdge)
    {
        INMOST::Tag coords_tag = mesh->GetGridInfo(SIMUG::mesh::gridElemType::Edge)->coords.at(coord_type);

        for (size_t bnd_edge_num = 0; bnd_edge_num < mesh->GetBndEdges().size(); ++bnd_edge_num)
        {
            if (mesh->GetBndEdges()[bnd_edge_num].GetStatus() != Element::Ghost)
            {
                double coord_x = mesh->GetBndEdges()[bnd_edge_num].RealArray(coords_tag)[0];
                double coord_y = mesh->GetBndEdges()[bnd_edge_num].RealArray(coords_tag)[1];

                if (output_size == 1)
                {
                    mesh->GetBndEdges()[bnd_edge_num].Real(var_tag) = func_ptr({coord_x, coord_y}, time)[0];
                }
                else
                {
                    for (size_t i = 0; i < output_size; ++i)
                    {
                        mesh->GetBndEdges()[bnd_edge_num].RealArray(var_tag)[i] = func_ptr({coord_x, coord_y}, time)[i];
                    }
                }
            }
        }
        ice_mesh->ExchangeData(var_tag, INMOST::FACE, 0);
    }
    else if (elem_type == SIMUG::mesh::gridElemType::bndTrian)
    {
        INMOST::Tag coords_tag = mesh->GetGridInfo(SIMUG::mesh::gridElemType::Trian)->coords.at(coord_type);

        for (size_t bnd_trian_num = 0; bnd_trian_num < mesh->GetBndTrians().size(); ++bnd_trian_num)
        {
            if (mesh->GetBndTrians()[bnd_trian_num].GetStatus() != Element::Ghost)
            {
                double coord_x = mesh->GetBndTrians()[bnd_trian_num].RealArray(coords_tag)[0];
                double coord_y = mesh->GetBndTrians()[bnd_trian_num].RealArray(coords_tag)[1];

                if (output_size == 1)
                {
                    mesh->GetBndTrians()[bnd_trian_num].Real(var_tag) = func_ptr({coord_x, coord_y}, time)[0];
                }
                else
                {
                    for (size_t i = 0; i < output_size; ++i)
                    {
                        mesh->GetBndTrians()[bnd_trian_num].RealArray(var_tag)[i] = func_ptr({coord_x, coord_y}, time)[i];
                    }
                }
            }
        }
        ice_mesh->ExchangeData(var_tag, INMOST::CELL, 0);
    }
    else 
    {
        SIMUG_ERR("unknown type of element!");
    }
    BARRIER
}

void Forcing::Update(mesh::meshVar mesh_var,
                     int layer,
                     coord::coordType coord_type,
                     double time)
{
    // figure out if it prognostic of forcing variable 
    if (mesh::progDims.count(mesh_var) == 1)
    {
        // get element type
        mesh::gridElemType elem_type = mesh->GetMeshInfo().prog_elems[mesh_var];

        // get expression
        std::optional<FuncPtr> expr_ptr = mesh->GetProgData(elem_type, layer)->GetExpr(mesh_var);

        // check if expression assigned
        if (!expr_ptr)
        {
            SIMUG_ERR("expression is not assigned - can not update mesh variable!");
        }

        // get variable tag
        INMOST::Tag var_tag = mesh->GetProgData(elem_type, layer)->Get(mesh_var);

        // update prognostic value 
        UpdateVariable(var_tag, elem_type, coord_type, expr_ptr.value(), time);
    }
    else if (mesh::forcDims.count(mesh_var) == 1)
    {
        // get element type
        mesh::gridElemType elem_type = mesh->GetMeshInfo().forc_elems[mesh_var];

        // get expression
        std::optional<FuncPtr> expr_ptr = mesh->GetForcData(elem_type)->GetExpr(mesh_var);

        // check if expression assigned
        if (!expr_ptr)
        {
            SIMUG_ERR("expression is not assign - can not update mesh variable!");
        }

        // get variable tag
        INMOST::Tag var_tag = mesh->GetForcData(elem_type)->Get(mesh_var);

        // update forcing value 
        UpdateVariable(var_tag, elem_type, coord_type, expr_ptr.value(), time);
    }
    else
    {
        SIMUG_ERR("unknown mesh variable!");
    }
    BARRIER
}

void Forcing::Update(mesh::meshVar mesh_var,
                     coord::coordType coord_type,
                     double time)
{
    Update(mesh_var, 0, coord_type, time);
    BARRIER
}

void Forcing::Update(const std::string& varname,
                     mesh::gridElemType elem_type,
                     int layer,
                     coord::coordType coord_type,
                     double time)
{
    // get expression
    std::optional<FuncPtr> expr_ptr = mesh->GetProgData(elem_type, layer)->GetExpr(varname);

    // check if expression assigned
    if (!expr_ptr)
    {
        SIMUG_ERR("expression is not assigned - can not update mesh variable!");
    }

    // get variable tag
    INMOST::Tag var_tag = mesh->GetProgData(elem_type, layer)->Get(varname);

    // update prognostic value 
    UpdateVariable(var_tag, elem_type, coord_type, expr_ptr.value(), time);
    BARRIER
}

void Forcing::Update(const std::string& varname,
                     mesh::gridElemType elem_type,
                     coord::coordType coord_type,
                     double time)
{
    Update(varname, elem_type, 0, coord_type, time);
    BARRIER
}