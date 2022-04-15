#include <iostream>

#include "defines.hpp"
#include "variable.hpp"


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
    MeshDataVar var1("test variable", "m", mesh::meshDim::scalar);
    
    cout << var1 << endl;
    var1.SetName("new variable");
    cout << var1 << endl;
    var1.SetUnit("s");
    cout << var1 << endl;
    var1.SetDim(mesh::meshDim::tensor);
    cout << var1 << endl;
    if ( var1.GetType() != varType::varMeshData or
         var1.GetName() != "new variable" or
         var1.GetUnit() != "s" or
         var1.GetDim()  != mesh::meshDim::tensor)
    {
        return false;
    }
    return true;
}
/*
bool modvar_test()
{
    for (SIMUG::model::var item: SIMUG::model::vars)
    {
        std::cout << SIMUG::model::name[item] << ", ";
    }
    std::cout << std::endl;
    return true;
}

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
    /*
    if (advvar_test())
        std::cout << "Advvar test: OK!\n";
    else
        SIMUG_ERR("Advvar test: FAILED!\n");

    if (dynvar_test())
        std::cout << "Dynvar test: OK!\n";
    else
        SIMUG_ERR("Dynvar test: FAILED!\n");
    */

    if (meshvar_test())
        std::cout << "Meshvar test: OK!\n";
    else
        SIMUG_ERR("Meshvar test: FAILED!\n");

/*
    if (modvar_test())
        std::cout << "Modvar test: OK!\n";
    else
        SIMUG_ERR("Modvar test: FAILED!\n");
    
    if (physvar_test())
        std::cout << "Physvar test: OK!\n";
    else
        SIMUG_ERR("Physvar test: FAILED!\n");

    if (confvar_test())
        std::cout << "Confvar test: OK!\n";
    else
        SIMUG_ERR("Confvar test: FAILED!\n");

*/
    return 0;

#if defined (USE_MPI)
    MPI_Finalize()
#endif
}