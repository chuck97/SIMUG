#include <iostream>
#include "defines.hpp"
#include "allvars.hpp"
#include "model_config.hpp"
#include "phys_config.hpp"
#include "mesh_config.hpp"
#include "adv_config.hpp"
#include "dyn_config.hpp"
#include "forc_config.hpp"


const std::string test_filename = "/data90t/geosci/spetrov/SIMUG/Configs/test_config.json";

using namespace SIMUG;

bool ModelConfig_tests()
{
    ModelConfig<double> MC1;
    std::cout << "It is " << conf::name[MC1.GetConfigType()] << "\n";
    MC1.Log(std::cout);
    ModelConfig<float> MC2(test_filename);
    MC2.Log(std::cout);
    Var<float> test_var = {30.0, "m"};
    MC2.SetVar(model::var::dt, test_var);
    MC2.Log(std::cout);
    return true;
};

bool PhysConfig_tests()
{
    PhysConfig<double> PC1;
    std::cout << "It is " << conf::name[PC1.GetConfigType()] << "\n";
    PC1.Log(std::cout);
    PhysConfig<float> PC2(test_filename);
    PC2.Log(std::cout);
    Var<float> test_var = {20.0, "m s-2"};
    PC2.SetVar(phys::var::g, test_var);
    PC2.Log(std::cout);
    return true;
};

bool MeshConfig_tests()
{
    MeshConfig<double> MC1;
    std::cout << "It is " << conf::name[MC1.GetConfigType()] << "\n";
    MC1.Log(std::cout);
    MeshConfig<float> MC2(test_filename);
    MC2.Log(std::cout);
    Var<std::string> test_var = {"ice mass", "kg"};
    MC2.SetVar(mesh::var::mi, test_var);
    MC2.Log(std::cout);
    return true;
};

bool AdvConfig_tests()
{
    AdvConfig<double> AC1;
    std::cout << "It is " << conf::name[AC1.GetConfigType()] << "\n";
    AC1.Log(std::cout);
    AdvConfig<float> AC2(test_filename);
    AC2.Log(std::cout);
    AC2.SetScheme(adv::scheme::TTG4);
    AC2.SetFilter(adv::filter::none);
    AC2.Log(std::cout);
    return true;
};

bool DynConfig_tests()
{
    DynConfig<double> DC1;
    std::cout << "It is " << conf::name[DC1.GetConfigType()] << "\n";
    DC1.Log(std::cout);
    DynConfig<float> DC2(test_filename);
    DC2.Log(std::cout);
    DC2.SetScheme(dyn::scheme::aEVP);
    DC2.SetPress(dyn::press::repl);
    DC2.SetBC(dyn::bc::fric);
    DC2.Log(std::cout);
    return true;
};

bool ForcConfig_tests()
{
    ForcConfig<double, forc::ice> FC1;
    std::cout << "It is " << conf::name[FC1.GetConfigType()] << "\n";
    FC1.Log(std::cout);
    Var<std::string> test_var = {"mass", "kg m-2"};
    FC1.SetVar(mesh::var::mi, {test_var, forc::setup::file});
    FC1.Log(std::cout);
    ForcConfig<float, forc::water> FC2(test_filename);
    FC2.Log(std::cout);
    ForcConfig<float, forc::ice> FC3(test_filename);
    FC3.Log(std::cout);
    return true;
};


int main()
{
    std::cout << "\n";

    if (!ModelConfig_tests())
        SIMUG_ERR("Model config test: NOT PASSED!");
    
    if (!PhysConfig_tests())
       SIMUG_ERR("Phys config test: NOT PASSED!");
    
    if (!MeshConfig_tests())
       SIMUG_ERR("Mesh config test: NOT PASSED!");

    if (!AdvConfig_tests())
       SIMUG_ERR("Adv config test: NOT PASSED!");
    
    if (!DynConfig_tests())
       SIMUG_ERR("Dyn config test: NOT PASSED!");

    if (!ForcConfig_tests())
       SIMUG_ERR("Forc config test: NOT PASSED!");

    return 0;
}