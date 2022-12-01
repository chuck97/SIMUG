#include <iostream>

#include "simug.hpp"


using namespace std;
using namespace SIMUG;

/*
bool advvar_test()
{
    for (SIMUG::adv::scheme item: SIMUG::adv::schemes)
    {
        std::cout << SIMUG::adv::advschname[item] << ", ";
    }
    std::cout << std::endl;
    return true;
}

bool dynvar_test()
{
    for (SIMUG::dyn::scheme item: SIMUG::dyn::schemes)
    {
        std::cout << SIMUG::dyn::momschname[item] << ", ";
    }
    std::cout << std::endl;
    return true;
}
*/

bool meshvar_test()
{
    MeshDataVar var("test variable", "m", mesh::meshDim::scalar);
    
    cout << var << endl;
    var.SetName("new variable");
    cout << var << endl;
    var.SetUnit("s");
    cout << var << endl;
    var.SetDim(mesh::meshDim::tensor);
    cout << var << endl;

    if ( var.GetType() != varType::varMeshData or
         var.GetName() != "new variable" or
         var.GetUnit() != "s" or
         var.GetDim()  != mesh::meshDim::tensor)
    {
        return false;
    }
    return true;
}

bool meshinfo_test()
{
    MeshInfoVar var("mesh information", mesh::surfType::plane, mesh::gridType::Agrid);
    
    cout << var << endl;
    var.SetName("another mesh info");
    cout << var << endl;
    var.SetSurfType(mesh::surfType::basin);
    cout << var << endl;
    var.SetGridType(mesh::gridType::Cgrid);
    cout << var << endl;

    if ( var.GetType() != varType::varMeshInfo or
         var.GetName() != "another mesh info" or
         var.GetSurfType() != mesh::surfType::basin or
         var.GetGridType()  != mesh::gridType::Cgrid)
    {
        return false;
    }
    return true;
}

bool modvar_test()
{
    ModelVar<float> varf("float variable", "m s-1", 1e-5);
    std::cout << varf << std::endl;
    ModelVar<float> vari("int variable", "", 30);
    std::cout << vari << std::endl;
    ModelVar<std::string> vars("string variable", "", "meshes/mesh.txt");
    std::cout << vars << std::endl;
    vars.SetName("another string var");
    std::cout << vars << std::endl;
    vars.SetUnit("mesh path");
    std::cout << vars << std::endl;
    vars.SetValue("meshes/another_mesh.txt");
    std::cout << vars << std::endl;

    if ( vars.GetType() != varType::varModel or
         vars.GetName() != "another string var" or
         vars.GetUnit() != "mesh path" or
         vars.GetValue()  != "meshes/another_mesh.txt")
    {
        return false;
    }
    return true;
}

/*
bool physvar_test()
{
    for (SIMUG::phys::var item: SIMUG::phys::vars)
    {
        std::cout << SIMUG::phys::name[item] << ", ";
    }
    std::cout << std::endl;
    return true;
}

bool confvar_test()
{
    for (SIMUG::conf::type item: SIMUG::conf::types)
    {
        std::cout << SIMUG::conf::name[item] << ", ";
    }
    std::cout << std::endl;
    return true;
}
*/

int main(int argc, char *argv[])
{
#if defined (USE_MPI)
    MPI_Init()
#endif

    if (meshvar_test())
        std::cout << "MeshDataVar test: OK!\n";
    else
        SIMUG_ERR("MeshDataVar test: FAILED!\n");

    if (meshinfo_test())
        std::cout << "MeshInfoVar test: OK!\n";
    else
        SIMUG_ERR("MeshInfoVar test: FAILED!\n");

    if (modvar_test())
        std::cout << "ModelVar test: OK!\n";
    else
        SIMUG_ERR("ModelVar test: FAILED!\n");

    return 0;

#if defined (USE_MPI)
    MPI_Finalize()
#endif
}