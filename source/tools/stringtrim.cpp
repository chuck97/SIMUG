#include "stringtrim.hpp"

void SIMUG::tools::ltrim(std::string& s)
{
    s.erase(s.begin(), find_if(s.begin(), s.end(), 
    [](unsigned char ch) 
    {
        return !isspace(ch);
    }));
};

void SIMUG::tools::rtrim(std::string &s) 
{
    s.erase(find_if(s.rbegin(), s.rend(), 
    [](unsigned char ch) 
    {
        return !isspace(ch);
    }).base(), s.end());
};

void SIMUG::tools::trim(std::string &s) 
{
    ltrim(s);
    rtrim(s);
};