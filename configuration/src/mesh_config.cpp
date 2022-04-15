#include "mesh_config.hpp"

using namespace SIMUG;

template<typename RT>
MeshConfig<RT>::MeshConfig()
{
    this->config_type = conf::type::mesh;
};

template<typename RT>
MeshConfig<RT>::MeshConfig(const std::string& json_path_):
    Configuration<RT>(json_path_)
{
    this->config_type = conf::type::mesh;
    std::ifstream ifs(this->json_filename);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    nlohmann::json j_mesh = j_input[conf::jsonkey[conf::type::mesh]];

    // Parse all mesh data
    for (mesh::var item: mesh::vars)
    {
        // Parse all mesh variables
        if (!j_mesh["variables"][mesh::varname[item]].empty())
        {
            Var<std::string> cur_var = {mesh::varname[item],
                                        j_mesh["variables"][mesh::varname[item]]};
            
            MeshVars[item] = cur_var;
        }

        // Parse surface type
        if (!j_mesh["surface"].empty())
        {
            surface = mesh::surfnotation[j_mesh["surface"]];
        }
    }  
};

template<typename RT>
void MeshConfig<RT>::SetSurf(mesh::surf surf_)
{
    surface = surf_;
};

template<typename RT>
mesh::surf MeshConfig<RT>::GetSurf() const
{
    return surface;
};


template<typename RT>
void MeshConfig<RT>::SetVar(mesh::var var_notation_,
                            const Var<std::string>& var_)
{
    MeshVars[var_notation_] = var_;
};

template<typename RT>
void MeshConfig<RT>::SetVars(const std::map<mesh::var, Var<std::string>>& vars_)
{
    MeshVars = vars_;
};

template<typename RT>
Var<std::string> MeshConfig<RT>::GetVar(mesh::var var_notation_) const
{
    return MeshVars.at(var_notation_);
};

template<typename RT>
conf::type MeshConfig<RT>::GetConfigType() const
{
    return this->config_type;
};

template<typename RT>
void MeshConfig<RT>::Log(std::ostream& os) const
{
    os << "===================================================\n";
    os << "MESH VARIABLES:\n";
    os << "surface type: " << mesh::surfname[surface] << ";\n";
    int i = 1;
    for (const auto& [key, val]: MeshVars)
    {
        os << i << ". " << mesh::varname[key] << ": " << val.unit << ";\n";
        ++i;
    }
    os << "===================================================\n";
};

template class MeshConfig<float>;
template class MeshConfig<double>;