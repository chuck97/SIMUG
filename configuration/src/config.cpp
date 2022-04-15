#include "config.hpp"

using namespace SIMUG;

// empty default constructor
template <typename RT>
Configuration<RT>::Configuration()
{
};

// json-input constructor
template <typename RT>
Configuration<RT>::Configuration(const std::string& json_path_)
{
    std::string trimed_json = json_path_;
    tools::trim(trimed_json);
    json_filename = trimed_json;
    if(json_filename.substr(json_filename.size()-5, 5) != ".json")
        SIMUG_ERR("input configuration file shoud be ended by .json!");
    
    std::ifstream infile(json_filename);
    if (!infile.good())
        SIMUG_ERR("config file \"" + json_filename + "\" doesn't exists!");
};

// output: json-filename
template <typename RT>
std::string Configuration<RT>::GetJsonFilename() const
{
    if (json_filename.empty())
        return "";
    else
        return json_filename;
};

// templates instantiation
template struct SIMUG::Var<float>;
template struct SIMUG::Var<double>;
template struct SIMUG::Var<std::string>;

template class SIMUG::Configuration<float>;
template class SIMUG::Configuration<double>;