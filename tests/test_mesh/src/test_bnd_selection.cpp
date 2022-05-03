#include "tests.hpp"

using namespace INMOST;
using namespace SIMUG::mesh;
using namespace std;

bool test_bnd_selection()
{
    IceMesh imesh(MESH_PATH,
                 surfType::plane,
                 gridType::Agrid);

    INMOST::Tag tag_bnd_nodes = imesh.GetMesh()->CreateTag("bnd_nodes", DATA_INTEGER, INMOST::NODE, INMOST::NODE, 1);
    INMOST::Tag tag_bnd_trians = imesh.GetMesh()->CreateTag("bnd_trians", DATA_INTEGER, INMOST::CELL, INMOST::CELL, 1);


    for (auto bnodeit = imesh.GetBndNodes().begin(); bnodeit != imesh.GetBndNodes().end(); ++bnodeit)
        bnodeit->Integer(tag_bnd_nodes) = 1.0;

    for (auto btrianit = imesh.GetBndTrians().begin(); btrianit != imesh.GetBndTrians().end(); ++btrianit)
        btrianit->Integer(tag_bnd_trians) = 1.0;


    if (imesh.GetMesh()->GetProcessorsNumber() > 1)
        imesh.SaveVTU("./bndtest.pvtu");
    else
        imesh.SaveVTU("./bndtest.vtu");

    return true;
}