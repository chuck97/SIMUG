#include "mesh.hpp"

using namespace SIMUG::mesh;

void GridInfo::ExchangeAll(const unsigned char& GridElem)
{
    ice_mesh->ExchangeData(id, GridElem, 0);
    ice_mesh->ExchangeData(id_no_bnd, GridElem, 0);
    ice_mesh->ExchangeData(is_bnd, GridElem, 0);

    for (auto [key, value]: coords)
        ice_mesh->ExchangeData(value, GridElem, 0);
    BARRIER
}

void NodeInfo::Exchange()
{
    ExchangeAll(INMOST::NODE);
}

void EdgeInfo::Exchange()
{
    ExchangeAll(INMOST::FACE);
}

void TrianInfo::Exchange()
{
    ExchangeAll(INMOST::CELL);
}

void NodeInfo::Mute()
{
    ice_mesh->SetFileOption("Tag:id node", "nosave");
    ice_mesh->SetFileOption("Tag:id node no bnd", "nosave");
    ice_mesh->SetFileOption("Tag:is node bnd", "nosave");
    ice_mesh->SetFileOption("Tag:model coords node", "nosave");
    ice_mesh->SetFileOption("Tag:cart coords node", "nosave");
    ice_mesh->SetFileOption("Tag:geo coords node", "nosave");
    BARRIER
}

void EdgeInfo::Mute()
{
    ice_mesh->SetFileOption("Tag:id edge", "nosave");
    ice_mesh->SetFileOption("Tag:id edge no bnd", "nosave");
    ice_mesh->SetFileOption("Tag:is edge bnd", "nosave");
    ice_mesh->SetFileOption("Tag:model coords edge", "nosave");
    ice_mesh->SetFileOption("Tag:cart coords edge", "nosave");
    ice_mesh->SetFileOption("Tag:geo coords edge", "nosave");
    BARRIER
}

void TrianInfo::Mute()
{
    ice_mesh->SetFileOption("Tag:id trian", "nosave");
    ice_mesh->SetFileOption("Tag:id trian no bnd", "nosave");
    ice_mesh->SetFileOption("Tag:is trian bnd", "nosave");
    ice_mesh->SetFileOption("Tag:model coords trian", "nosave");
    ice_mesh->SetFileOption("Tag:cart coords trian", "nosave");
    ice_mesh->SetFileOption("Tag:geo coords trian", "nosave");
    BARRIER
}

