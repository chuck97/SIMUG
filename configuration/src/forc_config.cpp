#include "forc_config.hpp"

using namespace SIMUG;

template<typename RT, forc::forc ftype>
ForcConfig<RT, ftype>::ForcConfig()
{
    this->config_type = conf::type::forc;
    forctype = ftype; 
};

template<typename RT, forc::forc ftype>
ForcConfig<RT, ftype>::ForcConfig(const std::string& json_path_):
    Configuration<RT>(json_path_)
{
    this->config_type = conf::type::forc;
    forctype = ftype;
    std::ifstream ifs(this->json_filename);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    nlohmann::json j_forc = j_input[conf::jsonkey[conf::type::forc]];

    // Parse forcing
    if (!j_forc[forc::forcname[forctype]].empty())
    {      
        nlohmann::json j_forc_curr = j_forc[forc::forcname[forctype]];

        // parse forcing variables
        for (mesh::var item: mesh::vars)
        {
            if (!j_forc_curr[mesh::varname[item]].empty())
            {
                Var<std::string> cur_var = {mesh::varname[item],
                                            j_forc_curr[mesh::varname[item]][0]};
            
                ForcVars[item] = {cur_var, forc::forcsetupnotation.at(j_forc_curr[mesh::varname[item]][1])};
            }
        }  
    }
};

template<typename RT, forc::forc ftype>
void ForcConfig<RT, ftype>::SetVar(mesh::var mesh_var_notation_,
                            const std::pair<Var<std::string>, forc::setup>& mesh_var_)
{
    ForcVars[mesh_var_notation_] = mesh_var_;
};

template<typename RT, forc::forc ftype>
std::pair<Var<std::string>, forc::setup> 
ForcConfig<RT, ftype>::GetVar(mesh::var mesh_var_notation_)
{
    return ForcVars.at(mesh_var_notation_);
};

template<typename RT, forc::forc ftype>
conf::type ForcConfig<RT, ftype>::GetConfigType() const
{
    return this->config_type;
};

template<typename RT, forc::forc ftype>
forc::forc ForcConfig<RT, ftype>::GetForcingType() const
{
    return forctype;
};

template<typename RT, forc::forc ftype>
void ForcConfig<RT, ftype>::Log(std::ostream& os) const
{
    os << "===================================================\n";
    os << forc::forcname.at(forctype) <<" FORCING VARIABLES:\n";
    int i = 1;
    for (const auto& [key, val]: ForcVars)
    {
        os << i << ". " << mesh::varname[key] << ": " 
           << val.first.unit << ", " 
           << forc::forcsetupname.at(val.second) << ";\n";
        ++i;
    }
    os << "===================================================\n";
};

// template class instantiation
template class ForcConfig<float, forc::ice>;
template class ForcConfig<double, forc::ice>;

template class ForcConfig<float, forc::air>;
template class ForcConfig<double, forc::air>;

template class ForcConfig<float, forc::water>;
template class ForcConfig<double, forc::water>;