#include "defines.hpp"
#include "stringtrim.hpp"
#include <iostream>

bool trim_test()
{
    std::string a = " abcd ";
    std::string b = " abcd ";
    std::string c = " abcd ";

    SIMUG::tools::ltrim(a);
    SIMUG::tools::rtrim(b);
    SIMUG::tools::trim(c);
    
    if (a != "abcd ")
        return false;

    if (b != " abcd")
        return false;

    if (c != "abcd")
        return false;
    
    return true;
}

int main(int argc, char *argv[])
{
    if (trim_test())
        std::cout << "Trim test: OK!\n";
    else
        SIMUG_ERR("Trim test: FAILED!\n")
    return 0;
}