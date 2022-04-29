#include "tests.hpp"

using namespace INMOST;
using namespace SIMUG::mesh;
using namespace std;

bool test_mesh_load()
{
    IceMesh(MESH_PATH,
            surfType::plane,
            gridType::Agrid);
    
    return true;
}