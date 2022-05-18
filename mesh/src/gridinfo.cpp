#include "mesh.hpp"

using namespace SIMUG::mesh;

void NodeInfo::Mute()
{
    ice_mesh->SetFileOption("Tag:id node", "nosave");
    ice_mesh->SetFileOption("Tag:id node no bnd", "nosave");
    ice_mesh->SetFileOption("Tag:is node bnd", "nosave");
    ice_mesh->SetFileOption("Tag:model coords node", "nosave");
    ice_mesh->SetFileOption("Tag:cart coords node", "nosave");
    ice_mesh->SetFileOption("Tag:geo coords node", "nosave");

    ice_mesh->SetFileOption("Tag:geo basis x node", "nosave");
    ice_mesh->SetFileOption("Tag:geo basis y node", "nosave");
    ice_mesh->SetFileOption("Tag:geo basis z node", "nosave");

    ice_mesh->SetFileOption("Tag:cart basis x node", "nosave");
    ice_mesh->SetFileOption("Tag:cart basis y node", "nosave");
    ice_mesh->SetFileOption("Tag:cart basis z node", "nosave");

    ice_mesh->SetFileOption("Tag:geo to elem trans matr node", "nosave");
    ice_mesh->SetFileOption("Tag:elem to geo trans matr node", "nosave");

    BARRIER
}

void NodeInfo::UnMute()
{
    ice_mesh->SetFileOption("Tag:id node", "save");
    ice_mesh->SetFileOption("Tag:id node no bnd", "save");
    ice_mesh->SetFileOption("Tag:is node bnd", "save");
    ice_mesh->SetFileOption("Tag:model coords node", "save");
    ice_mesh->SetFileOption("Tag:cart coords node", "save");
    ice_mesh->SetFileOption("Tag:geo coords node", "save");

    ice_mesh->SetFileOption("Tag:geo basis x node", "save");
    ice_mesh->SetFileOption("Tag:geo basis y node", "save");
    ice_mesh->SetFileOption("Tag:geo basis z node", "save");

    ice_mesh->SetFileOption("Tag:cart basis x node", "save");
    ice_mesh->SetFileOption("Tag:cart basis y node", "save");
    ice_mesh->SetFileOption("Tag:cart basis z node", "save");

    ice_mesh->SetFileOption("Tag:geo to elem trans matr node", "save");
    ice_mesh->SetFileOption("Tag:elem to geo trans matr node", "save");

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

    ice_mesh->SetFileOption("Tag:geo basis x edge", "nosave");
    ice_mesh->SetFileOption("Tag:geo basis y edge", "nosave");
    ice_mesh->SetFileOption("Tag:geo basis z edge", "nosave");

    ice_mesh->SetFileOption("Tag:cart basis x edge", "nosave");
    ice_mesh->SetFileOption("Tag:cart basis y edge", "nosave");
    ice_mesh->SetFileOption("Tag:cart basis z edge", "nosave");

    ice_mesh->SetFileOption("Tag:geo to elem trans matr edge", "nosave");
    ice_mesh->SetFileOption("Tag:elem to geo trans matr edge", "nosave");

    BARRIER
}

void EdgeInfo::UnMute()
{
    ice_mesh->SetFileOption("Tag:id edge", "save");
    ice_mesh->SetFileOption("Tag:id edge no bnd", "save");
    ice_mesh->SetFileOption("Tag:is edge bnd", "save");
    ice_mesh->SetFileOption("Tag:model coords edge", "save");
    ice_mesh->SetFileOption("Tag:cart coords edge", "save");
    ice_mesh->SetFileOption("Tag:geo coords edge", "save");

    ice_mesh->SetFileOption("Tag:geo basis x edge", "save");
    ice_mesh->SetFileOption("Tag:geo basis y edge", "save");
    ice_mesh->SetFileOption("Tag:geo basis z edge", "save");

    ice_mesh->SetFileOption("Tag:cart basis x edge", "save");
    ice_mesh->SetFileOption("Tag:cart basis y edge", "save");
    ice_mesh->SetFileOption("Tag:cart basis z edge", "save");

    ice_mesh->SetFileOption("Tag:geo to elem trans matr edge", "save");
    ice_mesh->SetFileOption("Tag:elem to geo trans matr edge", "save");

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

    ice_mesh->SetFileOption("Tag:geo basis x trian", "nosave");
    ice_mesh->SetFileOption("Tag:geo basis y trian", "nosave");
    ice_mesh->SetFileOption("Tag:geo basis z trian", "nosave");

    ice_mesh->SetFileOption("Tag:cart basis x trian", "nosave");
    ice_mesh->SetFileOption("Tag:cart basis y trian", "nosave");
    ice_mesh->SetFileOption("Tag:cart basis z trian", "nosave");

    ice_mesh->SetFileOption("Tag:geo to elem trans matr trian", "nosave");
    ice_mesh->SetFileOption("Tag:elem to geo trans matr trian", "nosave");

    BARRIER
}

void TrianInfo::UnMute()
{
    ice_mesh->SetFileOption("Tag:id trian", "save");
    ice_mesh->SetFileOption("Tag:id trian no bnd", "save");
    ice_mesh->SetFileOption("Tag:is trian bnd", "save");
    ice_mesh->SetFileOption("Tag:model coords trian", "save");
    ice_mesh->SetFileOption("Tag:cart coords trian", "save");
    ice_mesh->SetFileOption("Tag:geo coords trian", "save");

    ice_mesh->SetFileOption("Tag:geo basis x trian", "save");
    ice_mesh->SetFileOption("Tag:geo basis y trian", "save");
    ice_mesh->SetFileOption("Tag:geo basis z trian", "save");

    ice_mesh->SetFileOption("Tag:cart basis x trian", "save");
    ice_mesh->SetFileOption("Tag:cart basis y trian", "save");
    ice_mesh->SetFileOption("Tag:cart basis z trian", "save");

    ice_mesh->SetFileOption("Tag:geo to elem trans matr trian", "save");
    ice_mesh->SetFileOption("Tag:elem to geo trans matr trian", "save");

    BARRIER
}


